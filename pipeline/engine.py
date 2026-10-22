"""Run every registry event over nClimGrid, one year at a time.

    python pipeline/engine.py grids 1991-2026 --jobs 3   # event DOY grids,
                                                         # + climatology sums
    python pipeline/engine.py normals                    # 1991-2020 summary
    python pipeline/engine.py samples 2008-2026          # driver curves at
                                                         # observation cells
Outputs go to pipeline/cache/:
  years/<y>.npz    day of year each event first happens per 1/12° cell
                   (0 = not reached), and end-of-data value of each general
                   driver (e.g. degree-days so far)
  state/<y>.npz    every accumulator's state at the last observed day, which
                   nowcast.py continues with forecasts and normals
  clim.npz         1991-2020 daily climatology of every accumulator's
                   increment, as 4-harmonic coefficients per cell
  normals.npz      mean / SD / % years reached of each event, 1991-2020
  samples/<y>.npz  daily curves of the candidate drivers (fit.CANDIDATES) at
                   cells with iNaturalist observations, for model fitting
"""
import argparse
import concurrent.futures as cf
import os
import sys
import tempfile
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import drivers  # noqa: E402
import nclimgrid  # noqa: E402
import registry  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "cache")
NORMAL_YEARS = (1991, 2020)
N_HARM = 4
N_PENT = 74  # 5-day periods in a year (the last takes leftover days)


def parse_years(tokens):
    out = []
    for t in tokens:
        a, _, b = t.partition("-")
        out += list(range(int(a), int(b or a) + 1))
    return out


def harm_basis(doy):
    """[1, cos(2πkd/365.25), sin(2πkd/365.25)] for k = 1..N_HARM."""
    w = 2 * np.pi * np.asarray(doy, dtype=np.float64) / 365.25
    cols = [np.ones_like(w)]
    for k in range(1, N_HARM + 1):
        cols += [np.cos(k * w), np.sin(k * w)]
    return np.stack(cols, -1)


class DayStepper:
    """Shares state between events whose driver settings are identical."""

    def __init__(self, models, lat, shape, year):
        self.accs = {}
        for m in models:
            k = drivers.accumulator_key(m)
            if k not in self.accs:
                self.accs[k] = drivers.Accumulator(m, lat, shape, year)

    def step(self, doy, tmax, tmin, prcp):
        for a in self.accs.values():
            a.step(doy, tmax, tmin, prcp)

    def value(self, model):
        return self.accs[drivers.accumulator_key(model)].value


def month_urls(year):
    return [u for u in (nclimgrid.month_url(year, m) for m in range(1, 13)) if u]


def year_months(urls, year, workdir):
    with cf.ThreadPoolExecutor(3) as pool:
        futs = [pool.submit(nclimgrid.download, u, os.path.join(workdir, f"{year}-{os.path.basename(u)}"))
                for u in urls]
        for fut in futs:
            path = fut.result()
            out = nclimgrid.read_month(path)
            os.remove(path)
            yield out


def climatology_keys(events, general):
    """Accumulators that get a daily climatology (everything but rainfall)."""
    keys = {}
    for m in [e["model"] for e in events] + [g["model"] for g in general]:
        if m["type"] in ("gdd", "cold"):
            keys.setdefault(drivers.accumulator_key(m), m)
    return keys


def run_year(year, mode, workdir, harm=None):
    """Process one year. `harm` (dict) collects climatology sums when given."""
    reg, _ = registry.load()
    general = reg.get("drivers", [])
    events = registry.grid_events() if mode == "grids" else []
    sample_cells, cand = None, []
    if mode == "samples":
        import fit
        cand = fit.CANDIDATES
        path = os.path.join(CACHE, "obs_cells", f"{year}.npy")
        if not os.path.exists(path):
            print(f"{year}: no observation cells, skipped", flush=True)
            return
        sample_cells = np.load(path)
    models = [e["model"] for e in events] + [c["model"] for c in cand] + [g["model"] for g in general]
    clim = climatology_keys(events, general) if harm is not None else {}

    stepper = doy = valid = lat = lon = curves = None
    day, last_date = 0, None
    urls = month_urls(year)
    for lat_, lon_, times, data in year_months(urls, year, workdir):
        if stepper is None:
            lat, lon = lat_, lon_
            shape = (len(lat), len(lon))
            stepper = DayStepper(models, lat, shape, year)
            doy = np.zeros((len(events),) + shape, np.uint16)
            valid = np.zeros(shape, bool)
            if sample_cells is not None:
                curves = np.full((len(sample_cells), 366, len(cand)), np.nan, np.float32)
            prev = {k: np.zeros(shape, np.float32) for k in clim}
            if harm is not None and not harm:
                harm["count"] = np.zeros(N_PENT)
                harm["sums"] = {k: np.zeros((N_PENT,) + shape, np.float32) for k in clim}
                harm["models"] = dict(clim)
        for i in range(len(times)):
            tmax, tmin, prcp = data["tmax"][i], data["tmin"][i], data["prcp"][i]
            if np.isnan(tmax).all():
                continue
            day += 1
            valid |= ~np.isnan(tmax)
            stepper.step(day, tmax, tmin, prcp)
            for k, e in enumerate(events):
                v = stepper.value(e["model"])
                hit = (doy[k] == 0) & (v >= drivers.event_threshold(e["model"]))
                doy[k][hit] = day
            if curves is not None:
                for c, cm in enumerate(cand):
                    curves[:, day - 1, c] = stepper.value(cm["model"]).ravel()[sample_cells]
            if harm is not None:
                p = min(N_PENT - 1, (day - 1) // 5)
                harm["count"][p] += 1
                for key in clim:
                    v = stepper.accs[key].value
                    harm["sums"][key][p] += np.nan_to_num(v - prev[key])
                    prev[key] = v.copy()
            last_date = str(np.datetime_as_string(times[i], unit="D"))
    if stepper is None:
        print(f"{year}: no nClimGrid data yet, skipped", flush=True)
        return
    complete = (last_date or "").endswith("12-31")
    final = complete and len(urls) == 12 and all(u.endswith("-scaled.nc") for u in urls)
    meta = dict(lat=lat, lon=lon, last_date=last_date, last_doy=day, complete=complete, final=final)
    if mode == "grids":
        doy[:, ~valid] = 0
        out = os.path.join(CACHE, "years", f"{year}.npz")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        gen = np.stack([np.where(valid, stepper.value(g["model"]), np.nan) for g in general]) if general else np.zeros(0)
        np.savez_compressed(out, doy=doy, valid=valid, event_ids=np.array([e["key"] for e in events]),
                            general=gen.astype(np.float32), general_ids=np.array([g["id"] for g in general]), **meta)
        save_state(year, stepper, valid, meta)
    else:
        out = os.path.join(CACHE, "samples", f"{year}.npz")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        np.savez_compressed(out, curves=curves, cells=sample_cells,
                            drivers=np.array([c["id"] for c in cand]), **meta)
    print(f"{year}: {mode} -> {out} (through {last_date})", flush=True)


def save_state(year, stepper, valid, meta):
    """Accumulator states at the last observed day, for nowcast.py."""
    arrays = {"valid": valid}
    keys = []
    for n, (k, a) in enumerate(stepper.accs.items()):
        keys.append(k)
        arrays[f"value_{n}"] = a.value
        arrays[f"started_{n}"] = a.started
        if a.t == "precip":
            arrays[f"window_{n}"] = a.window
    path = os.path.join(CACHE, "state", f"{year}.npz")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez_compressed(path, keys=np.array(keys), **arrays, **meta)


def run_chunk(years, mode, workdir):
    """Worker: a run of years; returns climatology sums for normal-period years.

    Each year is retried a few times (network drops), and its climatology sums
    are merged only once the year has finished.
    """
    harm = {}
    for y in years:
        in_normal = mode == "grids" and NORMAL_YEARS[0] <= y <= NORMAL_YEARS[1]
        for attempt in range(3):
            yh = {} if in_normal else None
            try:
                run_year(y, mode, workdir, yh)
                break
            except Exception as e:  # noqa: BLE001 - retry any transient failure
                if attempt == 2:
                    raise
                print(f"{y}: {e!r}; retrying", flush=True)
        if yh:
            if not harm:
                harm.update(count=yh["count"].copy(), sums=yh["sums"], models=yh["models"])
            else:
                harm["count"] += yh["count"]
                for k in harm["sums"]:
                    harm["sums"][k] += yh["sums"][k]
    return harm


def merge_harm(parts):
    """Fit 4 harmonics to the 1991-2020 pentad-mean daily increments per cell."""
    parts = [p for p in parts if p]
    if not parts:
        return
    count = sum(p["count"] for p in parts)
    centers = np.arange(N_PENT) * 5 + 3.0
    w = np.sqrt(count)
    pinv = np.linalg.pinv(harm_basis(centers) * w[:, None])  # weight pentads by days seen
    coef = {}
    for key in parts[0]["models"]:
        s = sum(p["sums"][key].astype(np.float64) for p in parts)
        mean = s / np.maximum(count, 1)[:, None, None]
        y = (mean * w[:, None, None]).reshape(N_PENT, -1)
        coef[key] = (pinv @ y).reshape((2 * N_HARM + 1,) + s.shape[1:]).astype(np.float32)
    keys = list(coef)
    np.savez_compressed(os.path.join(CACHE, "clim.npz"), keys=np.array(keys),
                        coef=np.stack([coef[k] for k in keys]), n_harm=N_HARM)
    print(f"climatology: {len(keys)} accumulators")


def build_normals(y0=NORMAL_YEARS[0], y1=NORMAL_YEARS[1]):
    """Mean / SD of event DOY and % of years reached, per event."""
    years = list(range(y0, y1 + 1))
    first = np.load(os.path.join(CACHE, "years", f"{y0}.npz"))
    keys = list(first["event_ids"])
    s1 = np.zeros((len(keys),) + first["valid"].shape, np.float64)
    s2 = np.zeros_like(s1)
    n = np.zeros(s1.shape, np.int16)
    for y in years:
        z = np.load(os.path.join(CACHE, "years", f"{y}.npz"))
        if list(z["event_ids"]) != keys:
            raise SystemExit(f"{y}: events differ from {y0}; rebuild grids for all years")
        d = z["doy"].astype(np.float64)
        hit = d > 0
        s1 += np.where(hit, d, 0)
        s2 += np.where(hit, d * d, 0)
        n += hit
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mean = s1 / n
        sd = np.sqrt(np.maximum(0, (s2 - n * mean ** 2) / (n - 1)))
    half = n < len(years) / 2
    mean[half] = np.nan
    sd[half] = np.nan
    np.savez_compressed(os.path.join(CACHE, "normals.npz"), mean=mean.astype(np.float32),
                        sd=sd.astype(np.float32), pct=(100.0 * n / len(years)).astype(np.float32),
                        event_ids=np.array(keys), valid=first["valid"], lat=first["lat"],
                        lon=first["lon"], years=np.array([y0, y1]))
    print(f"normals: {len(keys)} events")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["grids", "samples", "normals"])
    ap.add_argument("years", nargs="*")
    ap.add_argument("--jobs", type=int, default=1)
    a = ap.parse_args()
    if a.mode == "normals":
        return build_normals()
    years = parse_years(a.years)
    chunks = [list(c) for c in np.array_split(years, max(1, min(a.jobs, len(years))))]
    with tempfile.TemporaryDirectory() as tmp:
        if len(chunks) > 1:
            with cf.ProcessPoolExecutor(len(chunks)) as pool:
                parts = list(pool.map(run_chunk, chunks, [a.mode] * len(chunks), [tmp] * len(chunks)))
        else:
            parts = [run_chunk(chunks[0], a.mode, tmp)]
    if a.mode == "grids":
        merge_harm(parts)


if __name__ == "__main__":
    main()
