# NOAA Workshop

The primary objective is to equip national meteorological agencies with the skills to access, explore, and apply high-resolution reanalysis datasets. A specific focus is on understanding severe weather events, including extreme rainfall, high winds, and temperature anomalies.

---

## Focus Area 1 — Ground Observations Monitoring & QC

**Core Objective:** Equip NMHS participants with tools for quality-controlling ground station data and validating it against satellite products, enabling identification of network issues and building confidence in observational networks.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/TAHMO/NOAA/blob/workshop-2/01_Ground_QC_improved_v5.ipynb)

---

## Focus Area 2 — Next-Gen Precipitation for Trends & Extremes

**Core Objective:** Compare next-generation precipitation datasets for historical trends, anomalies, and extreme event performance, enabling NMHS participants to select optimal products for climate monitoring and risk assessment.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/TAHMO/NOAA/blob/workshop-2/NOAA_Rainfall_Skill_Explorer_v2.ipynb)

---

## Focus Area 3 — Temperature Quality & Microclimates

**Core Objective:** Demonstrate the advantages of high-resolution temperature data in capturing microclimates and computing derived metrics such as Potential Evapotranspiration (PET) for improved assessment of heat-related risks.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/TAHMO/NOAA/blob/workshop-2/03_Temperature_Quality_PET.ipynb)

---

## Focus Area 4 — Agro-Ecological Zoning with PyAEZ

**Core Objective:** Synthesise climate monitoring outputs from Focus Areas 1–3 into agricultural interpretation using PyAEZ (FAO/AIT), computing thermal climate classifications, Length of Growing Period (LGP), and crop suitability maps for Kenya, Somalia, Sudan, Ethiopia, Tanzania, Eritrea, and Djibouti.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/TAHMO/NOAA/blob/workshop-2/04_AgroEcological_Zoning_PyAEZ.ipynb)

---

## Data Sources & Resources

### Ground Weather Data
1. **TAHMO (Trans-African Hydro-Meteorological Observatory):** Provides high-resolution weather and climate data collected from a network of ground stations across Africa. [TAHMO Website](https://tahmo.org/)
2. **GHCNd (Global Historical Climatology Network daily):** Integrated database of daily climate summaries from land surface stations across the globe. [NCEI Ground Stations](https://www.ncei.noaa.gov/products/land-based-station/global-historical-climatology-network-daily)

### Satellite & Reanalysis Products
3. **IMERG Data:** High-resolution precipitation estimates provided by NASA.
4. **TAMSAT:** Satellite-based rainfall estimates, soil moisture estimates, and seasonal forecasts tailored for Africa. [TAMSAT Website](https://research.reading.ac.uk/tamsat/)
5. **CHIRPS:** Gridded rainfall time series combining 0.05° resolution satellite imagery with in-situ station data. [CHIRPS Website](https://www.icpac.net/data-center/chirps/)
6. **CBAM:** Analytical tools and data for monitoring climate and environmental conditions.
7. **ERA5 / ERA5-Land:** Comprehensive global atmospheric reanalysis from ECMWF. [ERA5 Website](https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5)

### Weather Forecasting
8. **Google Forecasting:** AI-driven ensemble forecasting models (GraphCast & GenCast) providing 10-day global forecasts.
9. **CBAM Forecasting:** Weather forecasts provided by Tomorrow.io.

---

## Requirements
- API access to TAHMO
- API access to CBAM
- Token access from NCEI
- Google Account
