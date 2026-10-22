"""Fit phenology models to iNaturalist observations and validate every model.

    python pipeline/fit.py prepare     # obs CSVs -> obs table + cells per year
    (python pipeline/engine.py samples 2008-2026)   # driver curves at those cells
    python pipeline/fit.py fit         # -> registry/fitted.json
    (python pipeline/engine.py grids ...)           # event grids incl. fitted events
    python pipeline/fit.py validate    # -> site/data/validation.json

Fitting (per species x stage). Each annotated observation gives the value of
a driver (e.g. degree-days above 5 °C since Jan 1) on the day it was taken.
Onset / peak / end thresholds are the 10th / 50th / 90th percentiles of those
values. The driver is chosen among CANDIDATES by year-grouped cross-validation
on the training years (<= TEST_START - 1); the chosen model is then scored on
the held-out years against a calendar-only model that predicts the same date
everywhere. A model is published only if its error is at most 90% of the calendar
model's and its onset-end window holds at least 70% of test records.
"""
import csv
import datetime as dt
import gzip
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import registry  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "cache")
GRID = json.load(open(os.path.join(HERE, "grid.json")))
TEST_START = 2024
MIN_TRAIN, MIN_TEST = 100, 30
SKILL = 0.9          # publish only if error is at most 90% of the calendar model's
MIN_COVERAGE = 0.7   # and the onset-end window holds at least 70% of test records
MAX_ERROR = 28       # and its typical error is within four weeks
Q = {"onset": 0.10, "peak": 0.50, "end": 0.90}

_g = lambda b: {"type": "gdd", "base": b, "upper": 35, "method": "triangle", "start": "01-01"}
CANDIDATES = [
    {"id": "gdd0", "family": "gdd", "label": "degree-days above 0 °C since Jan 1", "model": _g(0)},
    {"id": "gdd5", "family": "gdd", "label": "degree-days above 5 °C since Jan 1", "model": _g(5)},
    {"id": "gdd10", "family": "gdd", "label": "degree-days above 10 °C since Jan 1", "model": _g(10)},
    {"id": "shortday", "family": "day", "label": "day length shortening after midsummer",
     "model": {"type": "shortday"}},
    {"id": "cold15_aug1", "family": "cold", "label": "chilling below 15 °C since Aug 1",
     "model": {"type": "cold", "base": 15, "start": "08-01"}},
    {"id": "cold15_p13", "family": "cold", "label": "chilling below 15 °C once days are shorter than 13 h",
     "model": {"type": "cold", "base": 15, "start": {"daylength": 13, "decreasing": True}}},
    {"id": "cold10_aug1", "family": "cold", "label": "chilling below 10 °C since Aug 1",
     "model": {"type": "cold", "base": 10, "start": "08-01"}},
]


def in_season(stage, doy):
    """Drop records outside the stage's season (e.g. last year's dried fruit in January)."""
    if stage == "Colored Leaves":
        return doy >= 182
    if stage == "Fruiting":
        return doy >= 91
    return doy >= 32


def cell_of(lat, lon):
    i = int((GRID["north"] - lat) // GRID["dlat"])
    j = int((lon - GRID["west"]) // GRID["dlon"])
    if 0 <= i < GRID["nrows"] and 0 <= j < GRID["ncols"]:
        return i * GRID["ncols"] + j
    return -1


# ------------------------------------------------------------------ prepare
def prepare():
    reg, _ = registry.load()
    rows = []
    for k, sp in enumerate(reg["species"]):
        path = os.path.join(CACHE, "obs", f"{sp['id']}.csv.gz")
        if not os.path.exists(path):
            continue
        with gzip.open(path, "rt") as f:
            for r in csv.DictReader(f):
                if not r["stage"]:
                    continue
                d = dt.date.fromisoformat(r["date"][:10])
                cell = cell_of(float(r["lat"]), float(r["lon"]))
                if cell < 0:
                    continue
                rows.append((k, r["stage"], d.year, d.timetuple().tm_yday, cell, int(r["id"])))
    if not rows:
        raise SystemExit("no observations found under cache/obs")
    sp_i, stage, year, doy, cell, oid = zip(*rows)
    np.savez_compressed(os.path.join(CACHE, "obs_table.npz"), species=np.array(sp_i), stage=np.array(stage),
                        year=np.array(year), doy=np.array(doy), cell=np.array(cell), id=np.array(oid))
    os.makedirs(os.path.join(CACHE, "obs_cells"), exist_ok=True)
    year, cell = np.array(year), np.array(cell)
    for y in np.unique(year):
        np.save(os.path.join(CACHE, "obs_cells", f"{y}.npy"), np.unique(cell[year == y]))
    print(f"prepared {len(rows)} annotated observations, years {year.min()}-{year.max()}")


def load_obs():
    z = np.load(os.path.join(CACHE, "obs_table.npz"))
    return {k: z[k] for k in z.files}


# ------------------------------------------------------------------ fitting
def curves_for(obs_idx, obs):
    """Daily candidate-driver curves (n, 366, n_candidates) for the given observations."""
    out = np.full((len(obs_idx), 366, len(CANDIDATES)), np.nan, np.float32)
    for y in np.unique(obs["year"][obs_idx]):
        path = os.path.join(CACHE, "samples", f"{y}.npz")
        if not os.path.exists(path):
            continue
        z = np.load(path)
        lookup = {c: i for i, c in enumerate(z["cells"])}
        sel = np.where(obs["year"][obs_idx] == y)[0]
        rows = [lookup.get(c, -1) for c in obs["cell"][obs_idx][sel]]
        ok = np.array(rows) >= 0
        out[sel[ok]] = z["curves"][np.array(rows)[ok]]
    return out


def pred_doy(curve, thr):
    """First day (1-based) each curve reaches thr; 366 if it never does."""
    hit = np.nan_to_num(curve, nan=-1) >= thr
    return np.where(hit.any(1), hit.argmax(1) + 1, 366)


def thresholds(x):
    return {q: float(np.quantile(x, p)) for q, p in Q.items()}


def score(curve, doy, th):
    on, pk, en = (pred_doy(curve, th[q]) for q in ("onset", "peak", "end"))
    return {"mae": float(np.mean(np.abs(doy - pk))), "coverage": float(np.mean((doy >= on) & (doy <= en))),
            "width": float(np.mean(en - on)), "early": float(np.mean(doy < on))}


def null_score(train_doy, doy):
    q = {k: np.quantile(train_doy, p) for k, p in Q.items()}
    return {"mae": float(np.mean(np.abs(doy - q["peak"]))),
            "coverage": float(np.mean((doy >= q["onset"]) & (doy <= q["end"]))),
            "width": float(q["end"] - q["onset"]), "early": float(np.mean(doy < q["onset"]))}


def region_keys(year, cell):
    """Group observations by year and 2° block, the scale managers plan at."""
    if cell is None:
        return None
    i, j = cell // GRID["ncols"], cell % GRID["ncols"]
    return year.astype(np.int64) * 100000 + (i // 24) * 100 + (j // 24)


def group_error(doy, pred_on, pred_pk, keys, null_on, null_pk, min_n=8):
    """Mean |error| of onset (10th pct) and peak (median) per region-year group.

    Returns (model error, calendar-model error, number of groups)."""
    em, en = [], []
    for g in np.unique(keys):
        sel = keys == g
        if sel.sum() < min_n:
            continue
        o_on, o_pk = np.quantile(doy[sel], Q["onset"]), np.median(doy[sel])
        em.append((abs(o_on - np.median(pred_on[sel])) + abs(o_pk - np.median(pred_pk[sel]))) / 2)
        en.append((abs(o_on - null_on) + abs(o_pk - null_pk)) / 2)
    return (float(np.mean(em)) if em else np.nan, float(np.mean(en)) if en else np.nan, len(em))


def candidates_for(stage, doy):
    """Heat drivers for everything; day-length and chilling drivers too for
    stages that peak after midsummer."""
    late = stage == "Colored Leaves" or np.median(doy) > 182
    return [k for k, c in enumerate(CANDIDATES) if c["family"] == "gdd" or (late and c["family"] in ("cold", "day"))]


def fit_stage(curves, doy, year, stage, cell=None):
    """Choose a driver, fit onset/peak/end thresholds and validate them.

    Error is measured per region-year group (2° block x year): the observed
    10th-percentile and median dates against the predicted onset and peak.
    That is the signal a manager plans on, and it is not swamped by the
    weeks-long spread of single records within a stage."""
    cands = candidates_for(stage, doy)
    ok = ~np.isnan(curves[:, 0, cands[0]])
    curves, doy, year = curves[ok], doy[ok], year[ok]
    cell = cell[ok] if cell is not None else None
    keys = region_keys(year, cell)
    x_all = curves[np.arange(len(doy)), doy - 1]  # (n, n_candidates)
    train = year < TEST_START
    if (train.sum() < MIN_TRAIN or len(np.unique(year[train])) < 3) and len(doy) >= MIN_TRAIN + MIN_TEST \
            and len(np.unique(year)) >= 3:
        train = np.ones_like(train)  # recent invader: validate by cross-validation over all years
    test = ~train
    res = {"n_train": int(train.sum()), "n_test": int(test.sum()),
           "years": [int(year.min()), int(year.max())] if len(year) else None}
    if train.sum() < MIN_TRAIN or len(np.unique(year[train])) < 3:
        return {**res, "publish": False, "reason": "too few observations"}

    def evaluate(k, fit_sel, eval_sel):
        th = thresholds(x_all[fit_sel, k])
        on = pred_doy(curves[eval_sel, :, k], th["onset"])
        pk = pred_doy(curves[eval_sel, :, k], th["peak"])
        nul_on, nul_pk = np.quantile(doy[fit_sel], Q["onset"]), np.median(doy[fit_sel])
        g = group_error(doy[eval_sel], on, pk, keys[eval_sel], nul_on, nul_pk) if keys is not None else (np.nan, np.nan, 0)
        obs_mae = float(np.mean(np.abs(doy[eval_sel] - pk)))
        return g, obs_mae, th

    # Driver choice: year-grouped CV on training years.
    tr_years = np.unique(year[train])
    folds = np.array_split(np.random.default_rng(0).permutation(tr_years), min(5, len(tr_years)))
    cv, cv_null = {}, []
    for k in cands:
        gm, gn, om = [], [], []
        for f in folds:
            va = train & np.isin(year, f)
            fi = train & ~va
            if va.sum() == 0 or fi.sum() < 20:
                continue
            (m, n0, ng), o, _ = evaluate(k, fi, va)
            om.append(o)
            if ng:
                gm.append(m)
                gn.append(n0)
        use_groups = len(gm) >= 2
        cv[CANDIDATES[k]["id"]] = float(np.mean(gm)) if use_groups else (float(np.mean(om)) if om else np.inf)
        if k == cands[0]:
            cv_null = gn if use_groups else []
    best = min(cv, key=cv.get)
    k = next(i for i, c in enumerate(CANDIDATES) if c["id"] == best)
    res.update({"driver": best, "cv_error": cv, "cv_error_null": float(np.mean(cv_null)) if cv_null else None})
    if test.sum() >= MIN_TEST:
        (gm, gn, ng), _, th_train = evaluate(k, train, test)
        res["test"] = score(curves[test, :, k], doy[test], th_train)
        res["test"]["region_error"], res["test"]["regions"] = gm, ng
        res["test_null"] = null_score(doy[train], doy[test])
        res["test_null"]["region_error"] = gn
        res["validated_by"] = f"held-out years {TEST_START}+"
        if ng >= 5:
            skill = gm <= SKILL * gn and gm <= MAX_ERROR
        else:  # few region-years: per-record test, and no worse where regions exist
            skill = res["test"]["mae"] <= SKILL * res["test_null"]["mae"] and res["test"]["mae"] <= MAX_ERROR
            if ng:
                skill = skill and gm <= gn
        skill = skill and res["test"]["coverage"] >= MIN_COVERAGE
    else:
        res["validated_by"] = "year-grouped cross-validation (too few records for a held-out test)"
        skill = res["cv_error_null"] is not None and cv[best] <= SKILL * res["cv_error_null"] and cv[best] <= MAX_ERROR
    all_th = thresholds(x_all[:, k])
    res.update({"publish": bool(skill), "thresholds": all_th, "driver_model": CANDIDATES[k]["model"],
                "driver_label": CANDIDATES[k]["label"],
                "reason": None if skill else "not clearly better than a calendar-only model"})
    return res


def fit():
    reg, _ = registry.load()
    obs = load_obs()
    out = {}
    for k, sp in enumerate(reg["species"]):
        stages = sp.get("fit_stages", [])
        for st in stages:
            idx = np.where((obs["species"] == k) & (obs["stage"] == st))[0]
            idx = idx[[in_season(st, d) for d in obs["doy"][idx]]] if len(idx) else idx
            if len(idx) == 0:
                out.setdefault(sp["id"], {})[st] = {"n_train": 0, "n_test": 0, "publish": False,
                                                     "reason": "no annotated observations"}
                continue
            cur = curves_for(idx, obs)
            r = fit_stage(cur, obs["doy"][idx], obs["year"][idx], st, obs["cell"][idx])
            out.setdefault(sp["id"], {})[st] = r
            t = r.get("test", {})
            print(f"{sp['id']:>18} {st:<18} n={r['n_train']}+{r['n_test']} driver={r.get('driver')} "
                  f"publish={r['publish']} region_err={t.get('region_error')} "
                  f"null={r.get('test_null', {}).get('region_error')} groups={t.get('regions')}", flush=True)
    with open(registry.FITTED, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, default=float)


# ------------------------------------------------------------------ validation
def validate():
    """Check literature onset events against observations of the matching stage,
    and collect fitted-model metrics, into site/data/validation.json."""
    reg, fitted = registry.load()
    obs = load_obs()
    out = {}
    years = {}
    for k, sp in enumerate(reg["species"]):
        for st, ev_id in sp.get("validate", {}).items():
            key = f"{sp['id']}.{ev_id}"
            idx = np.where((obs["species"] == k) & (obs["stage"] == st))[0]
            diffs, unreached = [], 0
            for y in np.unique(obs["year"][idx]):
                if y not in years:
                    p = os.path.join(CACHE, "years", f"{y}.npz")
                    years[y] = np.load(p) if os.path.exists(p) else None
                z = years[y]
                if z is None:
                    continue
                e = list(z["event_ids"]).index(key)
                sel = idx[obs["year"][idx] == y]
                pred = z["doy"][e].ravel()[obs["cell"][sel]].astype(int)
                unreached += int((pred == 0).sum())
                diffs.append(obs["doy"][sel][pred > 0] - pred[pred > 0])
            d = np.concatenate(diffs) if diffs else np.array([])
            out.setdefault(sp["id"], {})[st] = {
                "kind": "literature", "event": ev_id, "n": int(len(d) + unreached),
                "early": float(np.mean(d < 0)) if len(d) else None,
                "early_7d": float(np.mean(d < -7)) if len(d) else None,
                "median_after_onset": float(np.median(d)) if len(d) else None,
                "unreached": unreached}
        for st, r in fitted.get(sp["id"], {}).items():
            out.setdefault(sp["id"], {})[st] = {"kind": "fitted", **{k: r.get(k) for k in (
                "n_train", "n_test", "driver", "driver_label", "publish", "reason", "test", "test_null",
                "validated_by", "cv_error_null", "years")}, "cv_error": (r.get("cv_error") or {}).get(r.get("driver"))}
    path = os.path.join(HERE, "..", "site", "data", "validation.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, default=float)
    print(f"wrote {path}")


if __name__ == "__main__":
    {"prepare": prepare, "fit": fit, "validate": validate}[sys.argv[1]]()
