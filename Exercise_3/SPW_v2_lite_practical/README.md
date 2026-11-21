# Practical 2: short-term forecasts and planting decisions (SPW v2 lite)

A teaching notebook on **SPW v2**, the rain rules GAP uses to tell farmers when enough rain is likely to plant. It runs
at one example point in Zambia, Malawi or Nigeria for the 2023 and 2024 seasons, compares SPW v2 with the original
40 mm rule, and checks it against what really happened: the rain, the soil and One Acre Fund farm yields in Malawi.
**SPW v2 lite** uses the rain rules only; GAP's full SPW also checks soil moisture and long-term rainfall.

## Contents

| File | What it is |
|---|---|
| `spw_practical_toolkit_v2.ipynb` | The notebook (Sections 1–8), not yet run: run each section in turn, starting with the setup cell |
| `Practical_2_tasks.docx` | The task sheet: questions for each section, five key questions and the group output page |
| `Practical_2_planting_window_sheet.docx` | Printable answer sheet for the group output (Section 8 A and B): space to draft the five key answers and the planting window before typing them into the notebook |
| `04_spw_intro_v2.pptx` | The SPW intro deck, extended with SPW v2 and GROW-Africa slides |
| `data/locations.csv` | The three example points (Chipata, Lilongwe, Kaduna) |
| `data/<country>_observed_rainfall.csv` | TAMSAT daily rainfall at the example point |
| `data/<country>_forecast_daily.csv` | ECMWF daily forecast rainfall at the point (10 members) |
| `data/<country>_forecast_support.csv` | The forecast members SPW v2 uses each day |
| `data/<country>_spw_daily.csv` | Daily original-rule flags and the soil moisture over the next 15 days (Section 5) |
| `data/malawi_farms_2023_24.csv` | One Acre Fund Malawi 2023/24 farms within 100 km of Lilongwe: planting date, yield (t/ha), distance (km). No IDs or coordinates (Section 6) |
| `data/africa_country_boundaries.geojson` | Country outlines for the location map |
| `scripts/extract_mw_farms.py` | How the farm file was made from GROW-Africa. Needs the full GROW-Africa dataset; not needed to run the notebook |

## Running it

1. Create a Python 3.12 environment and install `requirements.txt`.
2. Open `spw_practical_toolkit_v2.ipynb` in JupyterLab, keeping it next to the `data/` folder.
3. Run the setup cell, then choose the country (Zambia, Malawi or Nigeria) and season (2023 or 2024) once in Section 1.
   Every section follows that choice. Section 6 always uses the Malawi 2023/24 farms, so every team does it.

## How it is laid out

Each section has an objective box, at most one settings cell, a chart and table (code folded), a key-result box that
reads the numbers back in plain words, questions to discuss and one key question. Section 7 summarises both seasons.

| Section | Task | Key question |
|---|---|---|
| 1 Setup and your example point | Task 1 | — |
| 2 One rain check | Task 2 | Q1 |
| 3 Through the season: when did SPW v2 first say plant? | Task 3 | Q2 |
| 4 How good was the forecast part? | Task 4 | Q3 |
| 5 What was the soil like after SPW said plant? | Task 5 | Q4 |
| 6 Malawi farms: planting in line with SPW and yields (everyone) | Task 6 | Q5 |
| 7 What we learned | — | — |
| 8 Group output: the five key answers and a planting window | Group output | — |

Section 8 checks the format of the answers (every part filled in, numbers and dates where needed) and saves them as
`planting_window_<country>_<year>_<group>.md` to hand in. Each group needs its own GROUP name. An optional SPW test explorer
with sliders sits after Section 8 for groups who finish early.

## SPW v2 rules used

Rain window: the last 6 observed days + the next 5 forecast days. Inside the planting window from the GAP maize
calendar: chance of more than 30 mm (from a Normal curve fitted to the 10 member totals) at least 80 % → *rain test passed*; otherwise median 30 mm or more → *plant with
caution*; otherwise *wait*. Before the window: *prepare*. After it: a short late-planting check (1 observed + 3 forecast days, 5 mm, 50 %). The 30 mm threshold applies to Zambia, Malawi and Nigeria
(Emily Black, Table 4). GAP's soil-moisture thresholds (SPW v2.5, not used by the rules here) are 45 % of field
capacity for Zambia and Malawi and 70 % for Nigeria; Section 5 shows them for comparison.

## Data sources

- Rainfall: TAMSAT. Forecasts: ECMWF reforecasts (2025 archive).
- Farm yields: One Acre Fund records in GROW-Africa (Geyman et al., *Scientific Data*, 2025), open data.
