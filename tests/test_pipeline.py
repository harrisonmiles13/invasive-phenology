"""Pipeline and data checks. Run from the repo root: python tests/test_pipeline.py

Exits non-zero on failure. Unit checks always run; the site-data checks run
when site/data/manifest.json exists (after a build).
"""
import gzip
import json
import math
import os
import shutil
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "pipeline"))
import drivers  # noqa: E402
import fit  # noqa: E402
import registry  # noqa: E402

failures = []


def check(cond, msg):
    if not cond:
        failures.append(msg)


# Degree-day methods (worked by hand; same branches as DDRP's TriDD).
for tmax, tmin, want in [(8, -2, 0.0), (38, 36, 25.0), (20, 12, 6.0), (15, 5, 1.25), (40, 20, 19.375), (40, 5, 12.5)]:
    got = float(drivers.triangle_dd(np.float64(tmax), np.float64(tmin), 10, 35))
    check(math.isclose(got, want, abs_tol=1e-3), f"triangle_dd({tmax},{tmin}) = {got}, want {want}")
check(float(drivers.average_dd(np.float64(20), np.float64(10), 10)) == 5.0, "average_dd")
check(float(drivers.average_dd(np.float64(40), np.float64(20), 10, 30)) == 15.0, "average_dd with cutoff")

# Day length: ~12.2 h at the equinox everywhere, ~15 h at 40°N in June.
check(abs(float(drivers.daylength(40.0, 80)) - 12.2) < 0.25, "equinox day length at 40N")
check(abs(float(drivers.daylength(40.0, 172)) - 15.0) < 0.25, "solstice day length at 40N")
check(float(drivers.daylength(25.0, 172)) < float(drivers.daylength(49.0, 172)), "day length grows with latitude in June")

# Accumulators: precipitation window and photoperiod start.
lat = np.array([40.0])
acc = drivers.Accumulator({"type": "precip", "window": 3, "mm": 10, "start": "01-01", "end": "12-31"}, lat, (1, 1), 2021)
for d, p in enumerate([5, 0, 0, 4, 4, 4], start=1):
    acc.step(d, np.zeros((1, 1)), np.zeros((1, 1)), np.full((1, 1), p))
check(float(acc.value[0, 0]) == 12.0, f"3-day rain window = {acc.value[0, 0]}, want 12")
acc = drivers.Accumulator({"type": "gdd", "base": 10, "upper": None, "start": {"daylength": 13.5}}, lat, (1, 1), 2021)
for d in range(1, 120):
    acc.step(d, np.full((1, 1), 30.0), np.full((1, 1), 20.0), None)
check(acc.started.all() and 0 < float(acc.value[0, 0]) < 15 * 30, "photoperiod biofix starts accumulation in spring")

# Fitting recovers a planted threshold and beats the calendar model.
rng = np.random.default_rng(0)
n, days = 900, np.arange(1, 367)
rate = rng.uniform(4, 12, n)  # each site warms at its own pace
curves = np.zeros((n, 366, len(fit.CANDIDATES)), np.float32)
for k in range(len(fit.CANDIDATES)):
    curves[:, :, k] = (np.maximum(0, days[None] - 60) * rate[:, None]) * (1 + 0.3 * k)
true_doy = np.array([np.argmax(curves[i, :, 1] >= 800) + 1 for i in range(n)])
doy = np.clip(true_doy + rng.normal(0, 4, n).round().astype(int), 1, 366)
year = rng.choice([2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025], n)
r = fit.fit_stage(curves, doy, year, "Flowering")
check(r["publish"], "synthetic fit should publish")
check(r.get("test", {}).get("mae", 99) < r.get("test_null", {}).get("mae", 0), "fitted model beats calendar model")

# Registry integrity.
reg, fitted = registry.load()
ids = [s["id"] for s in reg["species"]]
check(len(ids) == len(set(ids)), "species ids unique")
for e in registry.all_events():
    check(e["model"]["type"] in ("gdd", "cold", "precip", "photoperiod"), f"{e['key']}: driver type")
    check(e["source"] in reg["sources"], f"{e['key']}: unknown source {e['source']}")

# JS helpers.
if shutil.which("node"):
    js = os.path.join(ROOT, "site", "js", "model.js")
    out = subprocess.run(["node", "-e", f"const M=require({json.dumps(js)});console.log(JSON.stringify([M.fmtDoy(60,2024),M.fmtDoy(60),M.decodeCurrent(150|(3<<9)),M.relDays(-2)]))"],
                         capture_output=True, text=True, check=True).stdout
    check(json.loads(out) == ["Feb 29", "Mar 1", {"doy": 150, "src": 3}, "2 days ago"], f"js helpers: {out}")

# Built site data.
data = os.path.join(ROOT, "site", "data")
if os.path.exists(os.path.join(data, "manifest.json")):
    man = json.load(open(os.path.join(data, "manifest.json")))
    g = man["grid"]
    npx = g["nrows"] * g["ncols"]

    def layer(p):
        return np.frombuffer(gzip.decompress(open(os.path.join(data, p), "rb").read()), dtype="<u2")

    for sp in man["species"]:
        for e in sp["events"]:
            check(layer(f"ev/{sp['id']}/{e['id']}.normal.bin").size == 3 * npx, f"{sp['id']}.{e['id']} normal size")
            cur = layer(f"ev/{sp['id']}/{e['id']}.current.bin")
            check(cur.size == npx and ((cur >> 9) <= 4).all(), f"{sp['id']}.{e['id']} current layer")
    for d in man["drivers"]:
        a = layer(f"general/{d['id']}.bin")
        check(a.size == 3 * npx, f"general/{d['id']} size")
    c = man["current"]
    check(c["today_doy"] >= 1 and c["forecast_through_doy"] >= 1, "current dates")
    # Spotted lanternfly hatches earlier in Raleigh than in Reading, PA, in a normal year.
    i = lambda lat, lon: int((g["north"] - lat) // g["dlat"]) * g["ncols"] + int((lon - g["west"]) // g["dlon"])
    hatch = layer("ev/slf/hatch.normal.bin")
    check(0 < hatch[i(35.78, -78.64)] < hatch[i(40.34, -75.93)], "SLF hatch normals ordered south to north")

if failures:
    print("FAILED:\n  " + "\n  ".join(failures))
    sys.exit(1)
print("ok: drivers, day length, accumulators, fitting, registry, JS helpers, site data")
