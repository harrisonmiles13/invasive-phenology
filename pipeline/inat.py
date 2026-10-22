"""Fetch phenology-annotated iNaturalist observations for registry species.

    python pipeline/inat.py                 # all species -> pipeline/cache/obs/
    python pipeline/inat.py --only slf toh  # a subset

For each species and each stage (registry fit_stages plus validate keys) this
pulls research-grade US observations carrying that annotation, plus a capped
sample of all research-grade observations for the species' recorded range.
Only observations with open coordinates, positional accuracy under 5 km and
a CC0 / CC BY / CC BY-NC licence are kept, so they can be republished.

Runs in GitHub Actions (the build sandbox cannot reach the iNaturalist API).
Requests are spaced >= 1 s apart, per the API's recommended practice.
"""
import argparse
import csv
import gzip
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(__file__))
import registry  # noqa: E402

API = "https://api.inaturalist.org/v1"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "cache", "obs")
UA = "invasive-phenology/1.0 (github.com/harrisonmiles13/invasive-phenology)"
LICENSES = "cc0,cc-by,cc-by-nc"
PER_PAGE = 200
MAX_PER_STAGE = 15000  # newest first, so held-out recent years are always included
MAX_RANGE = 2000
_last = [0.0]


def get(path, **params):
    url = f"{API}/{path}?{urllib.parse.urlencode(params)}"
    for attempt in range(6):
        wait = 1.05 - (time.time() - _last[0])
        if wait > 0:
            time.sleep(wait)
        _last[0] = time.time()
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and attempt < 5:
                time.sleep(10 * (attempt + 1))
                continue
            raise
        except (urllib.error.URLError, TimeoutError):
            if attempt < 5:
                time.sleep(10 * (attempt + 1))
                continue
            raise


def controlled_terms():
    """{(term label, value label): (term_id, value_id)} for annotation filters."""
    out = {}
    for t in get("controlled_terms")["results"]:
        for v in t.get("values", []):
            out[(t["label"], v["label"])] = (t["id"], v["id"])
    return out


# Registry stage names -> iNaturalist value labels (iNat renamed the plant
# phenology values, e.g. "Flowering" -> "Flowers"; accept either).
STAGE_LABELS = {
    "Flowering": ("Flowering", "Flowers"),
    "Fruiting": ("Fruiting", "Fruits or Seeds"),
    "Flower Budding": ("Flower Budding", "Flower Buds"),
}


def stage_filter(terms, stage):
    for label in STAGE_LABELS.get(stage, (stage,)):
        for term in ("Life Stage", "Plant Phenology", "Flowers and Fruits", "Leaves"):
            if (term, label) in terms:
                return terms[(term, label)]
    return None


def taxon_id(name):
    res = get("taxa", q=name, is_active="true", per_page=30)["results"]
    for t in res:
        if t["name"].lower() == name.lower():
            return t["id"]
    raise LookupError(f"no iNaturalist taxon named {name!r}")


def observations(params, cap):
    """Yield (id, date, lat, lon), newest first, up to `cap` records."""
    base = dict(params, place_id=1, quality_grade="research", license=LICENSES, geoprivacy="open",
                taxon_geoprivacy="open", acc_below=5000, d1="2008-01-01", per_page=PER_PAGE,
                order_by="id", order="desc")
    last, n = None, 0
    while n < cap:
        page = get("observations", **base, **({"id_below": last} if last else {}))
        rows = page.get("results", [])
        if not rows:
            return
        for o in rows:
            last = o["id"] if last is None else min(last, o["id"])
            loc = o.get("location")
            day = o.get("observed_on")
            if not loc or not day:
                continue
            lat, lon = (float(v) for v in loc.split(","))
            n += 1
            yield o["id"], day, lat, lon
        if len(rows) < PER_PAGE:
            return


def fetch_species(sp, terms, log):
    tid = taxon_id(sp["inat_name"])
    stages = list(dict.fromkeys(sp.get("fit_stages", []) + list(sp.get("validate", {}))))
    path = os.path.join(OUT, f"{sp['id']}.csv.gz")
    counts = {}
    with gzip.open(path, "wt", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "date", "lat", "lon", "stage"])
        for st in stages:
            flt = stage_filter(terms, st)
            if not flt:
                log(f"  {sp['id']}: no annotation for {st!r}")
                continue
            k = 0
            for row in observations({"taxon_id": tid, "term_id": flt[0], "term_value_id": flt[1]}, MAX_PER_STAGE):
                w.writerow([*row, st])
                k += 1
            counts[st] = k
        k = 0
        for row in observations({"taxon_id": tid}, MAX_RANGE):
            w.writerow([*row, ""])
            k += 1
        counts["any"] = k
    log(f"{sp['id']} (taxon {tid}): {counts}")
    return {"taxon_id": tid, "counts": counts}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args()
    reg, _ = registry.load()
    os.makedirs(OUT, exist_ok=True)
    terms = controlled_terms()
    print("annotation terms:", sorted({t for t, _ in terms}), flush=True)
    summary = {}
    for sp in reg["species"]:
        if a.only and sp["id"] not in a.only:
            continue
        try:
            summary[sp["id"]] = fetch_species(sp, terms, lambda m: print(m, flush=True))
        except Exception as e:  # keep going; the summary records the failure
            print(f"{sp['id']}: FAILED {e}", flush=True)
            summary[sp["id"]] = {"error": str(e)}
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)


if __name__ == "__main__":
    main()
