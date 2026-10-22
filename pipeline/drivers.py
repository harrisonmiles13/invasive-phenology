"""Environmental drivers of phenology, evaluated day by day on a grid.

A phenological event is the first day a *driver* crosses a threshold. Drivers:

  gdd       growing degree-days above `base` (°C), optional `upper` cutoff,
            single-triangle or simple-average method, accumulated from a
            calendar date or from a photoperiod biofix (day length crossing
            `hours` while days lengthen).
  cold      chilling degree-days below `base` (°C mean), accumulated from a
            calendar date or from the day day length falls below `hours`
            (autumn senescence; after Delpierre et al. 2009).
  precip    rolling precipitation total (mm) over `window` days; the event is
            the first day inside [start, end] the total reaches `mm`.
  photoperiod  the day day length first reaches `hours` (increasing) or falls
            below it (decreasing).

Model specs are plain dicts so the same registry drives the Python grid
engine and the JavaScript point model (site/js/model.js).
"""
import numpy as np


def triangle_dd(tmax, tmin, base, upper=None):
    """Daily degree-days, single triangle, horizontal upper cutoff (DDRP TriDD)."""
    upper = 99.0 if upper is None else upper
    lo = np.minimum(tmax, tmin)
    hi = np.maximum(tmax, tmin)
    hi = np.where(hi == lo, hi + 0.01, hi)
    rng = hi - lo
    t1 = (hi - base) ** 2 / (2 * rng)
    t2 = (hi - upper) ** 2 / (2 * rng)
    mean_above = (hi + lo) / 2 - base
    dd = np.select(
        [hi < base, lo >= upper, (hi < upper) & (lo <= base), (lo <= base) & (hi >= upper),
         (lo > base) & (hi >= upper)],
        [0.0, upper - base, t1, t1 - t2, mean_above - t2],
        default=mean_above)
    return np.where(np.isnan(tmax) | np.isnan(tmin), np.nan, dd)


def average_dd(tmax, tmin, base, upper=None):
    """Simple-average degree-days; tmax capped at `upper` when given."""
    hi = tmax if upper is None else np.minimum(tmax, upper)
    return np.maximum(0.0, (hi + tmin) / 2 - base)


def daylength(lat_deg, doy):
    """Hours of daylight (sunrise to sunset, upper limb + refraction).

    CBM model, Forsythe et al. 1995 (p = 0.8333°). `lat_deg` may be an array.
    """
    theta = 0.2163108 + 2 * np.arctan(0.9671396 * np.tan(0.00860 * (doy - 186)))
    phi = np.arcsin(0.39795 * np.cos(theta))
    lat = np.radians(lat_deg)
    p = np.radians(0.8333)
    x = (np.sin(p) + np.sin(lat) * np.sin(phi)) / (np.cos(lat) * np.cos(phi))
    return 24 - (24 / np.pi) * np.arccos(np.clip(x, -1, 1))


def mmdd_to_doy(mmdd, year):
    import datetime as dt
    m, d = (int(v) for v in mmdd.split("-"))
    return (dt.date(year, m, d) - dt.date(year, 1, 1)).days + 1


def accumulator_key(model):
    """Models that share an accumulator differ only in their threshold."""
    m = dict(model)
    m.pop("threshold", None)
    m.pop("label", None)
    if m["type"] == "precip":
        m.pop("mm", None)
    return repr(sorted((k, repr(v)) for k, v in m.items()))


class Accumulator:
    """Running state of one driver over a grid for one year."""

    def __init__(self, model, lat_rows, shape, year):
        self.m = model
        self.t = model["type"]
        self.year = year
        self.lat = lat_rows[:, None]  # (nrows, 1) broadcast over columns
        self.value = np.zeros(shape, np.float32)
        start = model.get("start", "01-01")
        self.start_doy = mmdd_to_doy(start, year) if isinstance(start, str) else None
        self.photo = start if isinstance(start, dict) else None
        self.started = np.zeros(shape, bool)
        if self.t == "precip":
            self.window = np.zeros((model["window"],) + shape, np.float32)
            self.end_doy = mmdd_to_doy(model.get("end", "12-31"), year)
        self.prev_daylen = None

    def step(self, doy, tmax, tmin, prcp):
        """Advance one day; returns the driver value after this day."""
        t = self.t
        if t == "precip":
            self.window[doy % self.window.shape[0]] = np.nan_to_num(prcp)
            in_season = self.start_doy <= doy <= self.end_doy
            self.value = self.window.sum(0) if in_season else np.zeros_like(self.value)
            return self.value
        # Start: calendar date or photoperiod crossing.
        if self.photo is not None:
            dl = daylength(self.lat, doy)
            h = self.photo["daylength"]
            if self.photo.get("decreasing"):
                hit = (dl < h) & (doy > 172)
            else:
                hit = (dl >= h) & (doy <= 200)
            self.started |= np.broadcast_to(hit, self.started.shape)
        elif doy >= self.start_doy:
            self.started[:] = True
        if t == "gdd":
            f = average_dd if self.m.get("method") == "average" else triangle_dd
            inc = f(tmax, tmin, self.m["base"], self.m.get("upper"))
        elif t == "cold":
            inc = np.maximum(0.0, self.m["base"] - (tmax + tmin) / 2)
        elif t == "photoperiod":
            dl = np.broadcast_to(daylength(self.lat, doy), self.value.shape)
            if self.m.get("direction") == "decreasing":
                self.value = np.where((dl < self.m["hours"]) & (doy > 172), 1.0, 0.0).astype(np.float32)
            else:
                self.value = np.where(dl >= self.m["hours"], 1.0, 0.0).astype(np.float32)
            return self.value
        else:
            raise ValueError(f"unknown driver type {t}")
        inc = np.where(self.started, np.nan_to_num(inc), 0.0)
        self.value = (self.value + inc).astype(np.float32)
        return self.value


def event_threshold(model):
    """The level of the accumulator at which the event happens."""
    t = model["type"]
    if t == "precip":
        return model["mm"]
    if t == "photoperiod":
        return 1.0
    return model["threshold"]
