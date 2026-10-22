# Invasive Phenology Atlas

**Live app:** https://harrisonmiles13.github.io/invasive-phenology/ ·
**Manager guides:** https://harrisonmiles13.github.io/invasive-phenology/guides.html

Daily degree-days and phenology predictions for invasive insects and plants
across the contiguous United States. Load it any day to see:

- **Today**: heat accumulated since Jan 1 (base 32 °F and 50 °F), against the
  1991–2020 normal, and how many days ahead of or behind normal the season is.
- **Species**: when each modeled event (emergence, flowering, egg laying, fall
  color …) happens this year, marked observed, provisional, forecast (National
  Weather Service, ~10 days) or projected with normal temperatures. Also the
  normal date, its year-to-year variability, the share of years it happens,
  and a season animation.
- **A location**: what is coming up in the next three weeks for every species
  recorded nearby, and a season calendar against normal.

## Drivers

| Driver | Used for |
|---|---|
| Degree-days (single triangle or average, optional upper cutoff), from Jan 1 or a day-length biofix | insect development, plant leaf-out, flowering, fruiting |
| Day length (exact, by latitude) | biofix cues (brown marmorated stink bug), senescence |
| Chilling degree-days after a date or once days shorten | fall leaf color |
| Rolling rainfall | rain-triggered green-up (buffelgrass) |

## Models

- **Published models** (`registry/species.json`): DDRP parameter files
  (Barker & Coop, [bbarker505/ddrp_v3](https://github.com/bbarker505/ddrp_v3)) for spotted
  lanternfly, emerald ash borer, Asian longhorned beetle, Japanese beetle and
  light brown apple moth; the Utah State extension model for brown marmorated
  stink bug; the USA-NPN rain rule for buffelgrass.
- **Fitted models** (`registry/fitted.json`, refit monthly): research-grade iNaturalist
  observations carrying life-stage or plant-phenology annotations
  (CC0 / CC BY / CC BY-NC only). For each species × stage, onset / peak / end
  are the 10th / 50th / 90th percentiles of the driver value on observation
  days. Candidate drivers are degree-days above 0, 5 or 10 °C for every stage,
  plus shortening day length and three chilling rules for stages that peak
  after midsummer; the driver is chosen by year-grouped cross-validation.
  Error is scored per region-year (2° block × year): observed onset (10th
  percentile) and peak (median) dates against predicted ones. On held-out years
  (2024 onward) a model is published only if that error is at most 90% of a
  calendar-only model's, is within four weeks, and its onset–end window holds at
  least 70% of records.
- **Validation** (`site/data/validation.json`): published models are checked
  against iNaturalist records of the matching stage (share observed before the
  predicted onset); fitted models report held-out error, coverage and window
  width against the calendar model. The app shows these beside each prediction.

## Data flow

```
NOAA nClimGrid-Daily (AWS, 1/24°) ──► engine.py ──► per-year event grids (1/12°)
                                           │          normals 1991–2020
                                           │          daily climatology per driver
iNaturalist API ──► inat.py ──► fit.py ────┘
NWS National Blend of Models (AWS) ──► nowcast.py ──► this year: observed → forecast → projected
                                                          │
                                       build_site.py ◄────┘ ──► site/data/ ──► GitHub Pages
```

| Workflow | When | What |
|---|---|---|
| `daily.yml` | every day | current year from nClimGrid + NBM forecast + normals; republish |
| `rebuild.yml` | monthly, on registry change, on demand | fetch iNaturalist, refit, recompute 1991–present, validate, republish |

Generated layers live on the `gh-pages` branch, not `main`.

## Running locally

```bash
pip install -r pipeline/requirements.txt
python pipeline/engine.py grids 1991-2026 --jobs 3   # ~40 min, streams ~25 GB from NOAA
python pipeline/engine.py normals
python pipeline/nowcast.py                            # NBM forecast + projection
python pipeline/build_site.py --normals --years 2023-2025
python pipeline/assemble_page.py
python tests/test_pipeline.py
python -m http.server -d site 8000
```

Fitting needs the iNaturalist API: `python pipeline/inat.py`, then
`fit.py prepare`, `engine.py samples 2008-2026`, `fit.py fit`, then the grid
steps above, then `fit.py validate` and `build_site.py --obs`.

## Adding a species

Add an entry to `registry/species.json`: either `events` with a published
model (driver type, parameters, threshold, source) or `fit_stages` naming
iNaturalist annotations (`Flowering`, `Fruiting`, `Breaking Leaf Buds`,
`Colored Leaves`, `Larva`, `Nymph`, `Adult`, `Egg`). Pushing the change runs the
rebuild.

## Limits

Contiguous US only (nClimGrid domain). 8 km cells average terrain and
microclimate. iNaturalist effort is uneven, and fitted dates describe when people
photograph a stage. Use predictions to plan, and confirm stages in the field.

## Sources

NOAA nClimGrid-Daily (Durre et al. 2022); NOAA National Blend of Models;
Barker, Beyer & Coop 2025 (*Insects* 16:790) and DDRP parameter files;
Utah State University Extension BMSB model and Nielsen et al. 2008; Gerst et al.
2021 (USA-NPN buffelgrass forecast); iNaturalist contributors. Basemap ©
OpenStreetMap contributors, © CARTO. Code MIT.
