# Practical 1: from weather to agricultural risk

A teaching notebook on how weather and seasonal forecasts relate to national maize yield in Zambia, Malawi,
Zimbabwe, Botswana, Mozambique and Nigeria, 1984–2024. Work through all nine sections in order.

## Contents

| File | What it is |
|---|---|
| `yield_weather_practical_restored.ipynb` | The notebook, not yet run: run each section in turn, starting with the setup cell |
| `Practical_1_tasks.docx` | The task sheet: questions for each section |
| `Practical_1_country_brief_sheet.docx` | Printable country brief: the Section 9 questions with space to draft answers |
| `data/<country>.csv` | One row per harvest year: FAOSTAT maize yield, local yield rank, TAMSAT seasonal rainfall, ERA5 temperature, hot-day counts, SPI-3 |
| `data/<country>_monthly_weather.csv` | Monthly rainfall, SPI-3, temperature and hot days (Section 5) |
| `data/<country>_seasonal_forecasts_lead0.csv` | Lead-zero seasonal forecasts of rainfall and temperature, all ensemble members (Sections 3 and 7) |
| `data/oni_three_month.csv` | NOAA Oceanic Niño Index (Section 6) |

## Running it

1. Create a Python environment and install `requirements.txt`.
2. Open the notebook in JupyterLab, keeping it next to the `data/` folder.
3. Run the setup cell, then choose the country once in Section 1. Every section uses that choice.
4. After changing a setting in a section, run every cell below it in that section, so the key result matches the chart.

## How the notebook is laid out

Each section opens with an objective box (what it is for, what you will learn, which task it supports), then
a settings cell, the chart and a coloured table (code folded), a key-result box that reads the numbers back in
plain words, and the questions to discuss. Section 8 summarises what we learned for the chosen country.

## Group output: the country brief

Section 9 is the group output. Each group answers all six questions about its country (about 25 minutes). Answers need evidence
from the notebook (a year and a number), the group's own interpretation (something the notebook does not show),
and what it means for farmers or advice. Run the answer cell, then the check cell: it checks the format (every
part filled in, SMS length) and saves the brief as `country_brief_<country>_<group>.md` to hand in. Each group
needs its own GROUP name. Each group then reports back for 3 minutes. The same questions
are on the last page of the task sheet.
