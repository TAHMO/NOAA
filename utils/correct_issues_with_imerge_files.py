import xarray as xr
import glob
import pandas as pd
import re

# Path to your NetCDF files
files = sorted(glob.glob("Datasets/IMERG_*.nc"))


def preprocess_with_filename(ds, fname):
    # Extract start and end dates from filename
    match = re.search(r"IMERG_(\d{4}-\d{2}-\d{2})_(\d{4}-\d{2}-\d{2})_", fname)
    if not match:
        raise ValueError(f"Filename format not recognized: {fname}")

    start_date = pd.to_datetime(match.group(1))
    end_date = pd.to_datetime(match.group(2))

    # Actual number of time steps in the file
    n = ds.sizes["time"]

    # 📝 Log the file info here
    expected_n = (end_date - start_date).days + 1
    print(f"{fname} — time steps in file: {n}, expected: {expected_n}")

    # Generate expected dates based on start date and n
    expected_dates = pd.date_range(start=start_date, periods=n, freq="D")

    # Assign the corrected time coordinate
    ds = ds.assign_coords(time=expected_dates)
    return ds



datasets = []
for f in files:
    ds = xr.open_dataset(f)
    ds = preprocess_with_filename(ds, f)
    datasets.append(ds)

# Merge everything along time dimension
ds_merged = xr.concat(datasets, dim="time")

# Optional: sort time just in case
ds_merged = ds_merged.sortby("time")

# Save merged file
ds_merged.to_netcdf("IMERG_Kenya_merged.nc")

print(ds_merged)
