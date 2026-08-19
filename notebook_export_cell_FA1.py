# @title ### 11) Export results to interactive HTML dashboard {"display-mode":"form"}
# @markdown Collects all QC results and generates a self-contained HTML dashboard.
# @markdown Open the downloaded file in any browser — no internet connection needed.

import json
import os
import re
from datetime import date
from google.colab import files

# ── 1. COLLECT STATION DATA FROM NOTEBOOK VARIABLES ───────────────────────────
# Pulls from the variables your existing cells already produce:
#   scores_df   → Step 10a  (station scores)
#   qc_summary  → Step 6c   (QC flags per station)
#   chirps_corr → Step 9    (correlation + bias per station)
#   region_precip_data → Step 6a (for missing % calculation)

station_records = []

try:
    # Build a unified record per station
    for _, row in scores_df.iterrows():
        code = str(row.get('station_id', row.get('code', row.name)))

        # Score
        score = float(row.get('score', row.get('confidence_score', 0)))

        # Correlation + bias from CHIRPS comparison (Step 9 output)
        corr = 0.0
        bias = 0.0
        if 'chirps_corr' in dir() and chirps_corr is not None and not chirps_corr.empty:
            if code in chirps_corr.index:
                corr = float(chirps_corr.loc[code, 'correlation']   if 'correlation'   in chirps_corr.columns else 0)
                bias = float(chirps_corr.loc[code, 'bias']          if 'bias'          in chirps_corr.columns else 0)

        # QC flags from summary (Step 6c output)
        neg_rain  = 0
        spikes    = 0
        flatlines = 0
        if 'qc_summary_df' in dir() and qc_summary_df is not None and not qc_summary_df.empty:
            mask = qc_summary_df['station_id'].astype(str) == code
            if mask.any():
                r = qc_summary_df[mask].iloc[0]
                neg_rain  = int(r.get('neg_rain',  r.get('negative_rain', 0)))
                spikes    = int(r.get('daily_spike', r.get('spikes', 0)))
                flatlines = int(r.get('flatline',   r.get('flatlines', 0)))

        # Missing % from precipitation data
        missing_pct = 0.0
        if 'region_precip_data' in dir() and region_precip_data is not None \
                and not region_precip_data.empty and code in region_precip_data.columns:
            missing_pct = round(region_precip_data[code].isna().mean() * 100, 1)

        station_records.append({
            'code':        code,
            'score':       round(score, 1),
            'corr':        round(corr, 3),
            'bias':        round(bias, 2),
            'missing_pct': missing_pct,
            'neg_rain':    neg_rain,
            'spikes':      spikes,
            'flatlines':   flatlines,
        })

    print(f"✅ Collected data for {len(station_records)} stations.")

except Exception as e:
    print(f"⚠️  Could not collect station data: {e}")
    print("   Dashboard will use demo data. Check that scores_df exists (Step 10a).")
    station_records = []

# ── 2. BUILD DASHBOARD PAYLOAD ─────────────────────────────────────────────────
payload = {
    'country':    COUNTRY,
    'start_date': start_date,
    'end_date':   end_date,
    'generated':  date.today().isoformat(),
    'stations':   station_records,
}

# ── 3. LOAD THE DASHBOARD HTML TEMPLATE ───────────────────────────────────────
# The template is downloaded from the shared Drive alongside the notebook.
template_path = os.path.join(
    '/content/drive/Shareddrives/NOAA-workshop2',
    'dashboard_fa1.html'
)

if not os.path.exists(template_path):
    print(f"⚠️  Template not found at {template_path}.")
    print("   Upload dashboard_fa1.html to the NOAA-workshop2 Shared Drive root.")
    print("   Dashboard not exported.")
else:
    with open(template_path, 'r', encoding='utf-8') as f:
        html = f.read()

    # Inject the data payload by replacing the placeholder line
    data_line  = 'const DASHBOARD_DATA = window.__FA1_DATA__ || null;'
    inject_line = f'const DASHBOARD_DATA = {json.dumps(payload, indent=2)};'
    html = html.replace(data_line, inject_line)

    # Write and download
    out_name = f'FA1_QC_Dashboard_{COUNTRY}_{date.today().isoformat()}.html'
    out_path = os.path.join(BASE_OUT, out_name)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(html)

    print(f"✅ Dashboard exported: {out_name}")
    print(f"   {len(station_records)} stations · {COUNTRY} · {start_date} → {end_date}")
    print("\n📥 Downloading now...")
    files.download(out_path)
    print("✅ Done. Open the downloaded .html file in any browser.")
