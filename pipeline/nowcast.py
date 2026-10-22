"""Bring the current year up to today and predict ahead.

    python pipeline/nowcast.py [--date YYYY-MM-DD]

Starts from each accumulator's state on the last day of NOAA nClimGrid data
(engine.py grids <this year> writes it), then
  1. fills the 2-4 day gap to today with the National Blend of Models (NBM)
     day-0 forecast issued on each of those days (provisional),
  2. runs today's NBM forecast out ~10 days (daily max/min temperature and
     6-hourly precipitation, regridded to 1/12°),
  3. projects degree-day and chilling drivers to Dec 31 with the 1991-2020
     daily climatology of each driver (clim.npz). Rain-driven events are not
     projected past the forecast.
Every event's date gets a source code: 1 observed, 2 provisional, 3 forecast,
4 projected from normals, 0 not expected this year.

Also writes "today" maps of the general degree-day drivers with their
1991-2020 normal to date and how many days ahead or behind normal they are.
"""
import argparse
import datetime as dt
import os
import sys
import urllib.error
import urllib.request

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import drivers  # noqa: E402
import registry  # noqa: E402
from engine import CACHE, harm_basis  # noqa: E402

NBM_URL = "https://noaa-nbm-grib2-pds.s3.amazonaws.com"
SRC_OBS, SRC_PROV, SRC_FCST, SRC_PROJ = 1, 2, 3, 4
MAX_LEAD = 10
BEYOND = 2000  # days-ahead flag: more heat than a normal full year


# ------------------------------------------------------------------ NBM
def _get(url, rng=None):
    hdr = {"Range": f"bytes={rng[0]}-{rng[1]}"} if rng else {}
    with urllib.request.urlopen(urllib.request.Request(url, headers=hdr), timeout=120) as r:
        return r.read()


class NBM:
    """Daily max/min temperature (°C) and precipitation (mm) from NBM core files."""

    def __init__(self, grid):
        self.g = grid
        self.index = None
        self.idx_cache = {}

    def _file(self, cycle, fh):
        return f"{NBM_URL}/blend.{cycle:%Y%m%d}/00/core/blend.t00z.core.f{fh:03d}.co.grib2"

    def _message(self, cycle, fh, var, period):
        url = self._file(cycle, fh)
        if url not in self.idx_cache:
            self.idx_cache[url] = _get(url + ".idx").decode().splitlines()
        lines = self.idx_cache[url]
        level = "surface" if var == "APCP" else "2 m above ground"
        for k, line in enumerate(lines):
            p = line.split(":")
            if len(p) > 5 and p[3] == var and p[4] == level and p[5] == period and (len(p) < 7 or p[6] == ""):
                end = int(lines[k + 1].split(":")[1]) - 1 if k + 1 < len(lines) else ""
                return _get(url, (int(p[1]), end))
        raise LookupError(f"{var} {period} not in {url}")

    def _regrid(self, msg):
        import eccodes
        h = eccodes.codes_new_from_message(msg)
        try:
            vals = eccodes.codes_get_values(h)
            miss = eccodes.codes_get(h, "missingValue")
            if self.index is None:
                lat = eccodes.codes_get_array(h, "latitudes")
                lon = eccodes.codes_get_array(h, "longitudes")
                lon = np.where(lon > 180, lon - 360, lon)
                i = np.floor((self.g["north"] - lat) / self.g["dlat"]).astype(int)
                j = np.floor((lon - self.g["west"]) / self.g["dlon"]).astype(int)
                ok = (i >= 0) & (i < self.g["nrows"]) & (j >= 0) & (j < self.g["ncols"])
                self.index = np.where(ok, i * self.g["ncols"] + j, -1)
        finally:
            eccodes.codes_release(h)
        n = self.g["nrows"] * self.g["ncols"]
        ok = (self.index >= 0) & (vals != miss)
        s = np.bincount(self.index[ok], vals[ok], n)
        c = np.bincount(self.index[ok], minlength=n)
        with np.errstate(invalid="ignore"):
            out = s / c
        return out.reshape(self.g["nrows"], self.g["ncols"]).astype(np.float32)

    def day(self, cycle, lead):
        """(tmax, tmin, prcp) for the day `lead` days after the 00z `cycle`."""
        h0 = 24 * lead
        tmin = self._regrid(self._message(cycle, h0 + 12, "TMIN", f"{h0}-{h0 + 12} hour min fcst")) - 273.15
        tmax = self._regrid(self._message(cycle, h0 + 24, "TMAX", f"{h0 + 12}-{h0 + 24} hour max fcst")) - 273.15
        prcp = 0
        for h in (6, 12, 18, 24):
            prcp = prcp + self._regrid(self._message(cycle, h0 + h, "APCP", f"{h0 + h - 6}-{h0 + h} hour acc fcst"))
        return tmax, tmin, prcp

    def latest_cycle(self, today):
        for back in range(3):
            c = today - dt.timedelta(days=back)
            try:
                _get(self._file(c, 24 * MAX_LEAD + 24) + ".idx")
                return c
            except urllib.error.HTTPError:
                continue
        return None


# ------------------------------------------------------------------ climatology
def clim_increment(coef, doy):
    """Normal daily increment of a driver on `doy` (never negative)."""
    b = harm_basis(doy)
    return np.maximum(0.0, np.tensordot(b, coef, axes=(0, 0)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--date", help="treat this date as today (default: UTC today)")
    a = ap.parse_args()
    today = dt.date.fromisoformat(a.date) if a.date else dt.datetime.now(dt.timezone.utc).date()
    year = today.year
    today_doy = today.timetuple().tm_yday
    st = np.load(os.path.join(CACHE, "state", f"{year}.npz"))
    grids = np.load(os.path.join(CACHE, "years", f"{year}.npz"))
    clim = np.load(os.path.join(CACHE, "clim.npz"))
    clim_coef = {k: clim["coef"][n] for n, k in enumerate(clim["keys"])}
    import json
    grid = json.load(open(os.path.join(os.path.dirname(__file__), "grid.json")))
    reg, _ = registry.load()
    events = registry.grid_events()
    general = reg.get("drivers", [])
    lat = st["lat"]
    valid = st["valid"]
    shape = valid.shape

    # Restore accumulators.
    by_key = {}
    for m in [e["model"] for e in events] + [g["model"] for g in general]:
        by_key.setdefault(drivers.accumulator_key(m), m)
    accs = {}
    for n, k in enumerate(st["keys"]):
        if k not in by_key:
            continue
        acc = drivers.Accumulator(by_key[k], lat, shape, year)
        acc.value = st[f"value_{n}"].copy()
        acc.started = st[f"started_{n}"].copy()
        if f"window_{n}" in st:
            acc.window = st[f"window_{n}"].copy()
        accs[k] = acc
    missing = [k for k in by_key if k not in accs]
    if missing:
        raise SystemExit(f"state/{year}.npz lacks {len(missing)} accumulators; rerun engine.py grids {year}")

    gid = list(grids["event_ids"])
    doy = np.stack([grids["doy"][gid.index(e["key"])] for e in events]).astype(np.uint16)
    src = np.where(doy > 0, SRC_OBS, 0).astype(np.uint8)
    thr = [drivers.event_threshold(e["model"]) for e in events]
    ekeys = [drivers.accumulator_key(e["model"]) for e in events]

    def record(d, code):
        for k in range(len(events)):
            hit = (doy[k] == 0) & (accs[ekeys[k]].value >= thr[k]) & valid
            doy[k][hit] = d
            src[k][hit] = code

    last_obs = int(st["last_doy"])
    nbm = NBM(grid)
    cycle = None
    try:
        cycle = nbm.latest_cycle(today)
    except Exception as e:
        print(f"NBM not reachable ({e}); projecting from {st['last_date']} with normals", flush=True)
    d = last_obs
    general_today = None
    fc_end = last_obs
    if cycle is not None:
        last_day = min(dt.date(year, 12, 31), cycle + dt.timedelta(days=MAX_LEAD))
        day = dt.date(year, 1, 1) + dt.timedelta(days=last_obs)
        while day <= last_day:
            dd = day.timetuple().tm_yday
            c = day if day < cycle else cycle  # each past day: that day's own 00z run
            try:
                tmax, tmin, prcp = nbm.day(c, (day - c).days)
            except Exception as e:
                print(f"NBM {day} (cycle {c}): {e}; stopping forecast here", flush=True)
                break
            for acc in accs.values():
                acc.step(dd, np.where(valid, tmax, np.nan), np.where(valid, tmin, np.nan), prcp)
            record(dd, SRC_PROV if day < today else SRC_FCST)
            d = fc_end = dd
            if dd == today_doy:
                general_today = {g["id"]: accs[drivers.accumulator_key(g["model"])].value.copy() for g in general}
            print(f"  {day}: {'provisional' if day < today else 'forecast'} (NBM {c})", flush=True)
            day += dt.timedelta(days=1)
    # Normals to Dec 31 (no rain projection).
    n_days = 366 if (year % 4 == 0 and (year % 100 or year % 400 == 0)) else 365
    for dd in range(d + 1, n_days + 1):
        for k, acc in accs.items():
            if acc.t in ("gdd", "cold") and k in clim_coef:
                acc.value = (acc.value + clim_increment(clim_coef[k], dd)).astype(np.float32)
            elif acc.t == "photoperiod":
                acc.step(dd, None, None, None)
        record(dd, SRC_PROJ)
        if dd == today_doy and general_today is None:
            general_today = {g["id"]: accs[drivers.accumulator_key(g["model"])].value.copy() for g in general}
    if general_today is None:  # data already past today (re-runs for past dates)
        general_today = {g["id"]: grids["general"][n] for n, g in enumerate(general)}

    # Normal to date and days ahead (+) / behind (-) for the general drivers.
    gen_out = {}
    for g in general:
        k = drivers.accumulator_key(g["model"])
        inc = np.stack([clim_increment(clim_coef[k], x) for x in range(1, 367)])
        cum = np.cumsum(inc, 0)
        normal_today = cum[today_doy - 1]
        v = general_today[g["id"]]
        reach = (cum < v[None]).sum(0) + 1  # DOY the normal curve reaches today's value
        ahead = np.where(valid & (normal_today > 1), today_doy - reach, np.nan)
        # Past a normal full year's total the comparison has no date: flag it.
        ahead = np.where(valid & (reach > 365), BEYOND, ahead)
        gen_out[g["id"]] = np.stack([np.where(valid, v, np.nan), np.where(valid, normal_today, np.nan),
                                     np.where(valid, ahead, np.nan)]).astype(np.float32)

    doy[:, ~valid] = 0
    out = os.path.join(CACHE, "current.npz")
    np.savez_compressed(out, doy=doy, src=src, event_ids=np.array([e["key"] for e in events]),
                        general=np.stack([gen_out[g["id"]] for g in general]),
                        general_ids=np.array([g["id"] for g in general]), year=year, today=str(today),
                        today_doy=today_doy, obs_through=str(st["last_date"]), forecast_through_doy=fc_end,
                        nbm_cycle=str(cycle) if cycle else "")
    print(f"wrote {out}: observed through {st['last_date']}, forecast through DOY {fc_end}, NBM {cycle}")


if __name__ == "__main__":
    main()
