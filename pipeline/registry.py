"""Species registry: literature models plus models fitted from iNaturalist.

registry/species.json   hand-curated species, literature event models,
                        and the iNaturalist stages to fit
registry/fitted.json    written by fit.py: per species and stage, the chosen
                        driver and the onset / peak / end thresholds, with
                        test-year validation; only models with skill are used
"""
import json
import os

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SPECIES = os.path.join(ROOT, "registry", "species.json")
FITTED = os.path.join(ROOT, "registry", "fitted.json")

STAGE_PHASE = {
    "Breaking Leaf Buds": "emergence", "Flower Budding": "reproduction", "Flowering": "reproduction",
    "Fruiting": "reproduction", "Colored Leaves": "senescence", "Larva": "development",
    "Nymph": "development", "Adult": "emergence", "Egg": "reproduction", "Pupa": "development",
}
QUANTS = (("onset", "begins"), ("peak", "peaks"), ("end", "ends"))


def load():
    with open(SPECIES, encoding="utf-8") as f:
        reg = json.load(f)
    fitted = {}
    if os.path.exists(FITTED):
        with open(FITTED, encoding="utf-8") as f:
            fitted = json.load(f)
    return reg, fitted


def stage_slug(stage):
    return stage.lower().replace(" ", "_")


def fitted_events(sp, fitted):
    """Onset / peak / end events for each fitted stage that has skill."""
    out = []
    for stage, fm in fitted.get(sp["id"], {}).items():
        if not fm.get("publish"):
            continue
        for q, verb in QUANTS:
            model = dict(fm["driver_model"], threshold=round(fm["thresholds"][q], 2))
            out.append({"id": f"{stage_slug(stage)}_{q}", "label": f"{stage} {verb}",
                        "phase": STAGE_PHASE.get(stage, "other"), "stage": stage, "quantile": q,
                        "model": model, "source": "inat_fit"})
    return out


def all_events():
    """Every event in the registry, literature first, with species context."""
    reg, fitted = load()
    events = []
    for sp in reg["species"]:
        for e in sp["events"] + fitted_events(sp, fitted):
            events.append({**e, "key": f"{sp['id']}.{e['id']}", "species": sp["id"]})
    return events


def grid_events():
    return all_events()
