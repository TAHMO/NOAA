import xarray as xr
import pandas as pd
from pathlib import Path
import re
from collections import defaultdict

# 📂 Set this to your IMERG folder
DATA_DIR = Path("Datasets")

# 🧪 Regex to parse start, end date, and country
# Example filename: IMERG_2019-12-31_2020-12-29_Kenya.nc
FILENAME_PATTERN = re.compile(r"IMERG_(\d{4}-\d{2}-\d{2})_(\d{4}-\d{2}-\d{2})_(\w+)\.nc")

def preprocess_with_filename(ds, fname):
    """Preprocess a dataset by assigning correct time coords from filename."""
    fname = Path(fname).name
    match = FILENAME_PATTERN.match(fname)
    if not match:
        print(f"❌ Skipping {fname} — invalid filename")
        return None, None

    start_str, end_str, country = match.groups()
    start_date = pd.to_datetime(start_str)
    end_date = pd.to_datetime(end_str)

    n = ds.sizes["time"]
    expected_n = (end_date - start_date).days + 1

    print(f"📝 {fname} | {country} | {start_date.date()} → {end_date.date()} | expected={expected_n}, actual={n}")

    # 🚨 skip if difference is too large
    if abs(n - expected_n) > 10:
        print(f"⚠️ Skipping {fname} — time mismatch too large")
        return None, country

    new_time = pd.date_range(start=start_date, periods=n, freq="D")
    ds = ds.assign_coords(time=new_time)

    return ds, country

def merge_by_country(data_dir: Path):
    """Merge all IMERG files grouped by country."""
    grouped_files = defaultdict(list)

    # 🧭 First, group files by country
    for f in data_dir.glob("IMERG_*.nc"):
        fname = f.name
        match = FILENAME_PATTERN.match(fname)
        if match:
            _, _, country = match.groups()
            grouped_files[country].append(f)

    print(f"🌍 Found countries: {', '.join(grouped_files.keys())}")

    # 📦 Process each country
    for country, files in grouped_files.items():
        datasets = []
        skipped = []

        for f in sorted(files):
            try:
                with xr.open_dataset(f) as ds:
                    processed, _ = preprocess_with_filename(ds, f)
                    if processed is not None:
                        datasets.append(processed)
                    else:
                        skipped.append(f.name)
            except Exception as e:
                print(f"❌ Error reading {f.name}: {e}")
                skipped.append(f.name)

        if not datasets:
            print(f"🚫 No valid datasets for {country}")
            continue

        # 🪄 Merge along time dimension
        print(f"📦 Merging {len(datasets)} files for {country}...")
        merged = xr.concat(datasets, dim="time").sortby("time")

        # 🧭 Determine overall date range for filename
        start = str(merged.time.values[0])[:10]
        end = str(merged.time.values[-1])[:10]
        output_file = data_dir / f"IMERG_{start}_{end}_{country}.nc"

        merged.to_netcdf(output_file)
        print(f"✅ Saved merged file for {country}: {output_file}")

        # 📝 Log skipped files
        if skipped:
            log_file = data_dir / f"merge_log_{country}.txt"
            with open(log_file, "w") as log:
                log.write("Skipped files:\n")
                for s in skipped:
                    log.write(f"{s}\n")
            print(f"🪵 Log saved: {log_file}")

# 🏁 Run the merging
merge_by_country(DATA_DIR)
