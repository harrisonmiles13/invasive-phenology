"""NOAA nClimGrid-Daily access (public AWS bucket), aggregated to 1/12°.

Monthly netCDF files hold daily tmax, tmin (°C) and prcp (mm) for the
contiguous US at 1/24°. Final ("scaled") files are preferred; recent months
fall back to "prelim", whose days past the last observation are -999.99.
"""
import os
import time
import urllib.error
import urllib.request
import warnings

import numpy as np
import xarray as xr

BUCKET = "https://noaa-nclimgrid-daily-pds.s3.amazonaws.com"
VARS = ("tmax", "tmin", "prcp")


def month_url(year, month):
    for kind in ("scaled", "prelim"):
        url = f"{BUCKET}/access/grids/{year}/ncdd-{year}{month:02d}-grd-{kind}.nc"
        try:
            with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60):
                return url
        except urllib.error.HTTPError as e:
            if e.code not in (403, 404):
                raise
    return None


def download(url, dest):
    for attempt in range(5):
        try:
            urllib.request.urlretrieve(url, dest)
            return dest
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2 ** (attempt + 1))


def block_mean(a):
    """(T, 596, 1385) -> (T, 298, 692): 2x2 nan-mean, north row first."""
    a = np.where(a < -900, np.nan, a)
    t, h, w = a.shape
    a = a[:, : h - h % 2, : w - w % 2].reshape(t, h // 2, 2, w // 2, 2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        m = np.nanmean(a, axis=(2, 4))
    return m[:, ::-1, :].astype(np.float32)


def grid_coords(ds):
    lat = ds["lat"].values
    lon = ds["lon"].values
    lat = lat[: len(lat) - len(lat) % 2].reshape(-1, 2).mean(1)[::-1]
    lon = lon[: len(lon) - len(lon) % 2].reshape(-1, 2).mean(1)
    return lat.astype(np.float64), lon.astype(np.float64)


def read_month(path):
    with xr.open_dataset(path) as ds:
        lat, lon = grid_coords(ds)
        data = {v: block_mean(ds[v].values) for v in VARS}
        times = ds["time"].values
    return lat, lon, times, data
