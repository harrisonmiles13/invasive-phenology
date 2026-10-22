"""Pack model output into the files the web app loads (site/data/).

Binary layers are gzip-compressed little-endian uint16 over the 1/12° grid,
row-major, north row first (see site/js/data.js).

  manifest.json                    species, events, sources, dates, grid
  validation.json                  (fit.py) model checks against iNaturalist
  ev/<sp>/<event>.normal.bin       3 layers: mean DOY x10, SD x10, % of years
  ev/<sp>/<event>.<year>.bin       DOY reached (0 = not reached)
  ev/<sp>/<event>.current.bin      this year: DOY | source << 9, where source is
                                   1 observed, 2 provisional, 3 forecast,
                                   4 projected from normals (0 = not expected)
  general/<driver>.bin             today: value, 1991-2020 normal to date
                                   (°C·days, 65535 = no data) and days ahead of
                                   normal + 1000
  tiles/<ty>_<tx>.bin              every event's normal mean, SD and current
                                   layer for a 32x32 block (point panel)
  obs/<sp>.json, range/<sp>.bin    iNaturalist records and recorded range
  _store/                          pipeline state reused by the daily update

    python pipeline/build_site.py --normals --years 2024 2025   # full
    python pipeline/build_site.py                               # daily
"""
import argparse
import csv
import datetime as dt
import glob
import gzip
import json
import os
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import registry  # noqa: E402
from engine import CACHE, parse_years  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SITE = os.path.join(HERE, "..", "site", "data")
GRID = json.load(open(os.path.join(HERE, "grid.json")))
SHAPE = (GRID["nrows"], GRID["ncols"])
TILE = 32
KEEP_YEARS = 3  # complete years kept for the year selector


def write_u16(path, arr):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    data = np.ascontiguousarray(arr, dtype="<u2").tobytes()
    with open(path, "wb") as f:
        f.write(gzip.compress(data, compresslevel=9, mtime=0))


def read_u16(path, layers=1):
    with open(path, "rb") as f:
        a = np.frombuffer(gzip.decompress(f.read()), dtype="<u2")
    return a.reshape((layers,) + SHAPE) if layers > 1 else a.reshape(SHAPE)


def scaled(a, k):
    return np.where(np.isnan(a), 0, np.clip(np.round(a * k), 0, 65535)).astype(np.uint16)


def ev_path(key, part):
    sp, ev = key.split(".", 1)
    return os.path.join(SITE, "ev", sp, f"{ev}.{part}.bin")


def pack_normals():
    z = np.load(os.path.join(CACHE, "normals.npz"))
    for n, key in enumerate(z["event_ids"]):
        write_u16(ev_path(key, "normal"), np.stack([scaled(z["mean"][n], 10), scaled(z["sd"][n], 10),
                                                    scaled(z["pct"][n], 1)]))
    print(f"normals: {len(z['event_ids'])} events")
    return [int(v) for v in z["years"]]


def pack_year(y):
    z = np.load(os.path.join(CACHE, "years", f"{y}.npz"))
    for n, key in enumerate(z["event_ids"]):
        write_u16(ev_path(key, str(y)), z["doy"][n])
    print(f"{y}: {len(z['event_ids'])} events")
    return {"last_date": str(z["last_date"]), "complete": bool(z["complete"]), "final": bool(z["final"])}


def pack_current():
    z = np.load(os.path.join(CACHE, "current.npz"))
    for n, key in enumerate(z["event_ids"]):
        write_u16(ev_path(key, "current"), z["doy"][n].astype(np.uint16) | (z["src"][n].astype(np.uint16) << 9))
    for n, gid in enumerate(z["general_ids"]):
        v, norm, ahead = z["general"][n]
        enc = lambda a: np.where(np.isnan(a), 65535, np.clip(np.round(a), 0, 65534)).astype(np.uint16)
        write_u16(os.path.join(SITE, "general", f"{gid}.bin"),
                  np.stack([enc(v), enc(norm), enc(np.where(np.isnan(ahead), np.nan, ahead + 1000))]))
    meta = {k: (z[k].item() if z[k].shape == () else z[k].tolist()) for k in
            ("year", "today", "today_doy", "obs_through", "forecast_through_doy", "nbm_cycle")}
    print(f"current: {meta}")
    return meta


def pack_obs():
    """Recorded range (cells within ~25 km of a record) and recent annotated records."""
    reg, _ = registry.load()
    n_obs = {}
    for sp in reg["species"]:
        path = os.path.join(CACHE, "obs", f"{sp['id']}.csv.gz")
        if not os.path.exists(path):
            continue
        rng = np.zeros(SHAPE, bool)
        recs = []
        with gzip.open(path, "rt") as f:
            for r in csv.DictReader(f):
                lat, lon = float(r["lat"]), float(r["lon"])
                i = int((GRID["north"] - lat) // GRID["dlat"])
                j = int((lon - GRID["west"]) // GRID["dlon"])
                if not (0 <= i < SHAPE[0] and 0 <= j < SHAPE[1]):
                    continue
                rng[max(0, i - 3):i + 4, max(0, j - 3):j + 4] = True
                if r["stage"]:
                    d = dt.date.fromisoformat(r["date"][:10])
                    recs.append([round(lon, 3), round(lat, 3), d.year, d.timetuple().tm_yday, r["stage"], int(r["id"])])
        write_u16(os.path.join(SITE, "range", f"{sp['id']}.bin"), rng.astype(np.uint16))
        recs.sort(key=lambda r: -r[5])
        by_stage = {}
        for r in recs:
            if len(by_stage.setdefault(r[4], [])) < 4000:
                by_stage[r[4]].append(r)
        out = [r for v in by_stage.values() for r in v]
        os.makedirs(os.path.join(SITE, "obs"), exist_ok=True)
        with open(os.path.join(SITE, "obs", f"{sp['id']}.json"), "w") as f:
            json.dump({"fields": ["lon", "lat", "year", "doy", "stage", "id"], "records": out}, f, separators=(",", ":"))
        n_obs[sp["id"]] = len(recs)
    return n_obs


def store_state():
    """Keep what the daily update needs alongside the site (not linked from the app)."""
    dst = os.path.join(SITE, "_store")
    os.makedirs(dst, exist_ok=True)
    z = np.load(os.path.join(CACHE, "clim.npz"))
    for n, k in enumerate(z["keys"]):
        np.save(os.path.join(dst, f"clim_{n}.npy"), z["coef"][n].astype(np.float16))
    with open(os.path.join(dst, "clim_keys.json"), "w") as f:
        json.dump([str(k) for k in z["keys"]], f)
    for p in glob.glob(os.path.join(CACHE, "obs", "*.csv.gz")):
        os.makedirs(os.path.join(dst, "obs"), exist_ok=True)
        shutil.copy(p, os.path.join(dst, "obs"))
    print("stored climatology and observations")


def restore_state():
    """Rebuild pipeline/cache/clim.npz from site/data/_store (daily update)."""
    src = os.path.join(SITE, "_store")
    keys = json.load(open(os.path.join(src, "clim_keys.json")))
    coef = np.stack([np.load(os.path.join(src, f"clim_{n}.npy")).astype(np.float32) for n in range(len(keys))])
    os.makedirs(CACHE, exist_ok=True)
    np.savez(os.path.join(CACHE, "clim.npz"), keys=np.array(keys), coef=coef)
    for p in glob.glob(os.path.join(src, "obs", "*.csv.gz")):
        os.makedirs(os.path.join(CACHE, "obs"), exist_ok=True)
        shutil.copy(p, os.path.join(CACHE, "obs"))
    print(f"restored climatology ({len(keys)} drivers)")


def build_tiles(events, species_ids):
    """Per 32x32 block: every event's [normal mean, SD, current], then one
    recorded-range layer per species (1 = recorded nearby, 0 = not, 2 = no data)."""
    n = len(events)
    cube = np.zeros((n, 3) + SHAPE, np.uint16)
    rng = np.full((len(species_ids),) + SHAPE, 2, np.uint16)
    for k, sid in enumerate(species_ids):
        p = os.path.join(SITE, "range", f"{sid}.bin")
        if os.path.exists(p):
            rng[k] = read_u16(p)
    for k, e in enumerate(events):
        p = ev_path(e["key"], "normal")
        if os.path.exists(p):
            cube[k, :2] = read_u16(p, 3)[:2]
        c = ev_path(e["key"], "current")
        if os.path.exists(c):
            cube[k, 2] = read_u16(c)
    land = cube[:, 2].max(0) > 0
    for old in glob.glob(os.path.join(SITE, "tiles", "*.bin")):
        os.remove(old)
    names = []
    for ty in range(0, SHAPE[0], TILE):
        for tx in range(0, SHAPE[1], TILE):
            if not land[ty:ty + TILE, tx:tx + TILE].any():
                continue
            block = np.zeros((n * 3 + len(species_ids), TILE, TILE), np.uint16)
            sub = cube[:, :, ty:ty + TILE, tx:tx + TILE]
            h, w = sub.shape[2], sub.shape[3]
            block[:n * 3, :h, :w] = sub.reshape(n * 3, h, w)
            block[n * 3:, :h, :w] = rng[:, ty:ty + TILE, tx:tx + TILE]
            name = f"{ty // TILE}_{tx // TILE}"
            write_u16(os.path.join(SITE, "tiles", f"{name}.bin"), block)
            names.append(name)
    return names


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--normals", action="store_true")
    ap.add_argument("--years", nargs="*", default=[])
    ap.add_argument("--obs", action="store_true", help="pack iNaturalist records and ranges")
    ap.add_argument("--store", action="store_true", help="save climatology/obs to site/data/_store")
    ap.add_argument("--restore", action="store_true", help="load climatology/obs from site/data/_store")
    a = ap.parse_args()
    if a.restore:
        return restore_state()
    reg, fitted = registry.load()
    events = registry.grid_events()
    mpath = os.path.join(SITE, "manifest.json")
    man = json.load(open(mpath)) if os.path.exists(mpath) else {"years": {}}
    if a.normals:
        man["normal_period"] = pack_normals()
    for y in parse_years(a.years):
        man["years"][str(y)] = pack_year(y)
    cur = pack_current()
    # Keep the newest complete years only, never the current year.
    done = sorted((int(y) for y, v in man["years"].items() if v["complete"] and int(y) != cur["year"]), reverse=True)
    for y in done[KEEP_YEARS:]:
        man["years"].pop(str(y))
        for p in glob.glob(os.path.join(SITE, "ev", "*", f"*.{y}.bin")):
            os.remove(p)
    man["years"].pop(str(cur["year"]), None)
    if a.obs:
        man["obs_counts"] = pack_obs()
    if a.store:
        store_state()
    sp_ids = [sp["id"] for sp in reg["species"]]
    man["tiles"] = build_tiles(events, sp_ids)
    man["tile"] = TILE
    man["tile_events"] = [e["key"] for e in events]
    man["tile_species"] = sp_ids
    man["current"] = cur
    man["grid"] = GRID
    man["built"] = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%MZ")
    man["sources"] = reg["sources"]
    man["drivers"] = reg.get("drivers", [])
    man["species"] = []
    for sp in reg["species"]:
        evs = [e for e in events if e["species"] == sp["id"]]
        if not evs:
            continue
        man["species"].append({
            **{k: sp.get(k) for k in ("id", "common", "scientific", "group", "regions", "inat_name")},
            "has_range": os.path.exists(os.path.join(SITE, "range", f"{sp['id']}.bin")),
            "events": [{k: e.get(k) for k in ("id", "label", "phase", "stage", "quantile", "source", "note", "model")}
                       for e in evs]})
    with open(mpath, "w", encoding="utf-8") as f:
        json.dump(man, f, indent=1, ensure_ascii=False)
    print(f"manifest: {len(man['species'])} species, {len(events)} events, {len(man['tiles'])} tiles")


if __name__ == "__main__":
    main()
