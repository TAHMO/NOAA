# NOAA Workshop

The primary objective is to equip national meteorological agencies with the skills to access, explore, and apply high-resolution reanalysis datasets. A specific focus is on understanding severe weather events, including extreme rainfall, high winds, and temperature anomalies.

---

## Focus Area 1 — Ground Observations Monitoring & QC

**Core Objective:** Equip NMHS participants with tools for quality-controlling ground station data and validating it against satellite products, enabling identification of network issues and building confidence in observational networks.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Nelson1798/NOAA/blob/main/01_Ground_QC_improved_v5.ipynb)

---

## Focus Area 2 — Next-Gen Precipitation for Trends & Extremes

**Core Objective:** Compare next-generation precipitation datasets for historical trends, anomalies, and extreme event performance, enabling NMHS participants to select optimal products for climate monitoring and risk assessment.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Nelson1798/NOAA/blob/main/NOAA_Rainfall_Skill_Explorer_v2.ipynb)

---

## Focus Area 3 — Temperature Quality & Microclimates

**Core Objective:** Demonstrate the advantages of high-resolution temperature data in capturing microclimates and computing derived metrics such as Potential Evapotranspiration (PET) for improved assessment of heat-related risks.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Nelson1798/NOAA/blob/main/03_Temperature_Quality_PET.ipynb)

---

## Focus Area 4 — Agro-Ecological Zoning with PyAEZ

**Core Objective:** Synthesise climate monitoring outputs from Focus Areas 1–3 into agricultural interpretation using PyAEZ (FAO/AIT), computing thermal climate classifications, Length of Growing Period (LGP), and crop suitability maps for Somalia, Sudan, Ethiopia, Tanzania, Eritrea, and Djibouti.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Nelson1798/NOAA/blob/main/04_AgroEcological_Zoning_PyAEZ.ipynb)

---

<!-- ## Service Account & Shared Drive

All notebooks authenticate with Google Earth Engine using a shared service account and load pre-computed datasets from the NOAA-workshop2 Shared Drive, so no manual data upload or API configuration is required. -->

<!-- | Resource | Details |
|---|---|
| Service account | `gee-notebook-user@natural-notch-435413-j3.iam.gserviceaccount.com` |
| Key file (Drive) | [Download](https://drive.google.com/file/d/181IKB3sJ3iUn6ZOZbg50htgH2JKcxFkT/view?usp=sharing) |
| Shared Drive | `Shareddrives/NOAA-workshop2` | -->
<br>
<br>
<br>

<!--  
Use the Google Colab as the notebook to be able to run the activities of the workshop<br>

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1hkuUbjgdmrFCER621sb3dSwkrW26EZGN?usp=sharing)

###  Objectives
1. Understanding and Accessing the data in Ghana, Kenya and ZImbabwe (CBAM, ERA5, CHIRPS, TAMSAT, TAHMO and GHCNd ground stations)
2. Applying the data to Severe weather forecasting
<br>

## Understanding the data
### Ground Weather Data
1. TAHMO (Trans-African Hydro-Meteorological Observatory)
  - Provides high-resolution weather and climate data collected from a network of ground stations across Africa. <br>
  [TAHMO Website](https://tahmo.org/)

2. Global Historical Climatology Network daily (GHCNd)
  - The Global Historical Climatology Network daily (GHCNd) is an integrated database of daily climate summaries from land surface stations across the globe
  - GHCNd contains records from more than 100,000 stations in 180 countries and territories.<br>
  [NCEI Ground Stations](https://www.ncei.noaa.gov/products/land-based-station/global-historical-climatology-network-daily)

### Satellite Products
3. IMERG Data
     - (Precipiation data provided by NASA)
5. TAMSAT (Tropical Applications of Meteorology using SATellite data)
    - Enhances the capacity of African meteorological agencies and other organisations by providing and supporting the use of satellite-based rainfall estimates, soil moisture estimates and forecasts, and related data products.<br>
  [TAMSAT Website](https://research.reading.ac.uk/tamsat/)

5. CHIRPS (Climate Hazards Group InfraRed Precipitation with Station Data)
    - CHIRPS incorporates 0.05° resolution satellite imagery with insitu station data to create gridded rainfall time series for trend analysis and seasonal drought monitoring. <BR>
  [CHIRPS Website](https://www.icpac.net/data-center/chirps/)

### Reanalysis Datasets

6. CBAM
    - Offers data and analytical tools for monitoring and assessing climate and environmental conditions.

7. ECMWF Reanalysis v5 (ERA5) and ERA5-Land
    - Provides a comprehensive, global reanalysis of historical weather and climate data from the European Centre for Medium-Range Weather Forecasts.<br>
  [ERA5 Website](https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5)

### Weather Forecasting

8. Google Forecasting 
   - An ensemble of Google's weather forecasting (GraphCast and GenCast)
   - Provides a 10 day forecast
9. CBAM Forecasting
    - This is provided by Tomorrow.io

### Requirements
- API access to TAHMO
- API Access to CBAM
- Token Access from NCEI
- Google Account
=======

