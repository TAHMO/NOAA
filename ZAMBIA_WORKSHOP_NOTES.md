# Southern Africa Preparation Notes (Zambia workshop, 24–29 Nov)

Branch: `zambia-workshop`, created from `main`. The East Africa configuration on
`main` and the Madagascar work on `madagascar-work` are untouched.

Countries added: **Zimbabwe, Malawi, Mozambique, Botswana, Democratic Republic
of the Congo (DRC)**. All four notebooks keep their existing workflow. On this
branch the dropdowns list the new countries first and default to **Zimbabwe**.
The East African countries are still selectable and behave as on `main`.

## 1. Selected date ranges

| Notebook | Period | Why |
|---|---|---|
| 01 Ground QC, 02 Rainfall Skill Explorer | **2023-07-01 → 2025-06-30** | Two full austral hydrological years (Jul–Jun) |
| 02 Case studies | Event windows (below) | All inside the main period |
| 03 Temperature & PET | **2023-10-01 → 2023-10-31** (comparison window **2023-10-25 → 2023-10-29**) | Late-Oct 2023 heatwave |
| 03 Step 7 / 04 AEZ long-term | Unchanged | Climatology, not event driven |

### Reasoning

* Southern Africa's main rainy season runs from **November to March**. A July–June
  window keeps each season whole and starts and ends in the dry season.
* The two seasons give a **drought-versus-flood contrast**:
  * **2023/24** was a severe El Niño season. Zimbabwe, Zambia, Malawi, Mozambique
    and Botswana got less than 20% of normal February rainfall, and Zimbabwe,
    Zambia and Malawi declared national disasters.
  * **2024/25** was wet, with Cyclone Chido and the Gaborone and Kinshasa floods.

  This is directly relevant to a Zambia-hosted workshop.
* Every country gets a well-documented event inside the window. A 2021–2023
  window (as used for Madagascar) had no comparable Zimbabwe or Botswana events.
* The window is 731 days long, the same as the East Africa window.

### IMERG V06 → V07 (code change)

GEE's `NASA/GPM_L3/IMERG_V06` collection ends on **2024-06-02** and is deprecated,
so it can't cover this window. On this branch, Notebook 02 now pulls
`NASA/GPM_L3/IMERG_V07`. The units are the same (half-hourly mm/hr, ×0.5 to get
mm), the band is now `precipitation`, and the result is still saved as
`precipitationCal`, so nothing downstream changes. V07 permanent products run
to 2025-09-30.

## 2. Events investigated

### Rainfall (Notebook 02 `EVENTS`, 30 km radius)

| Country | Event | Peak | Window | Point | Evidence |
|---|---|---|---|---|---|
| Zimbabwe | El Niño mid-season **drought** | Feb 2024 | 02-01 → 02-29 | Harare (-17.83, 31.05) | <20% of normal Feb rain; national disaster 3 Apr 2024; 2.7 M people needing food aid ([WWA](https://www.worldweatherattribution.org/?p=2673), [AP](https://www.nbcnews.com/news/world/zimbabwe-el-nino-africa-drought-disaster-humanitarian-aid-rcna146347)) |
| Malawi | Nkhotakota floods | 26–27 Feb 2024 | 02-23 → 03-02 | Nkhotakota (-12.93, 34.30) | ~80,000 affected, 8 deaths ([Save the Children](https://resourcecentre.savethechildren.net/document/2024-save-the-children-malawi-nkhotakota-floods-situation-report), [IFRC](https://go-api.ifrc.org/api/downloadfile/87312/MDRMW020dfr)) |
| Mozambique | Tropical Cyclone Chido | 15 Dec 2024 | 12-12 → 12-19 | Mecúfi, Cabo Delgado (-13.28, 40.57) | ~260 km/h winds, ~250 mm in a day; 94–120 deaths; ~110,000 homes destroyed ([Al Jazeera](https://www.aljazeera.com/amp/news/2024/12/19/cyclone-chido-death-toll-in-mozambique-rises-to-73)) |
| Botswana | Gaborone floods | 16–20 Feb 2025 | 02-15 → 02-22 | Gaborone (-24.65, 25.91) | ≥9 deaths, 5,000+ displaced; ~40-yr return period ([WWA](https://www.worldweatherattribution.org/increasing-extreme-rainfall-and-rapid-urbanisation-major-drivers-behind-gaborones-deadly-floods/)) |
| DRC | Kinshasa floods | 4–5 Apr 2025 | 04-02 → 04-09 | Ndjili River area (-4.38, 15.33) | Ndjili River overflowed; ≥33 deaths ([IOL](https://iol.co.za/news/world/2025-04-08-congos-capital-kinshasa-faces-crisis-as-flooding-death-toll-rises-to-33/)) |

The Zimbabwe case is deliberately a **drought**: it tests whether the products
agree on near-zero rain, not on a flood peak.

Considered but not used: Cyclone Freddy in Malawi/Mozambique (March 2023) and the
Kinshasa floods of 13 December 2022. Both are strong events but fall before the
chosen window. Freddy is already configured on the `madagascar-work` branch.

### Heat (Notebook 03 `HEAT_EVENTS`)

**Late-October 2023 El Niño heatwave.** Botswana's Department of Meteorological
Services issued a nationwide heatwave warning with forecast maxima of 39–42 °C
until at least 27–30 Oct 2023 ([Crisis24](https://crisis24.garda.com/alerts/2023/10/botswana-namibia-heatwave-forecast-across-much-of-both-countries-through-at-least-oct-30)).
Regional outlooks gave 35–40 °C, 2–8 °C above average, across Botswana,
Zimbabwe and Mozambique. All five countries use this event; the East African
countries keep March 2024.

## 3. Datasets

Same products as the existing workflow:

* **Ground:** TAHMO (default), GHCN-Daily, or custom national met-service CSVs.
* **Rainfall:** CHIRPS pentad, TAMSAT (PERSIANN-CDR fallback), **IMERG V07**,
  ERA5-Land.
* **Temperature:** ERA5 / ERA5-Land, CBAM.
* **Region:** OSM country polygon (Notebook 01; DRC is queried as "Democratic
  Republic of the Congo") and the bboxes below (Notebooks 02–04).

| Country | bbox [lon_min, lat_min, lon_max, lat_max] | Capital |
|---|---|---|
| Zimbabwe | [25.2, -22.5, 33.1, -15.6] | Harare |
| Malawi | [32.6, -17.2, 36.0, -9.3] | Lilongwe |
| Mozambique | [30.2, -26.9, 41.0, -10.4] | Maputo |
| Botswana | [19.9, -27.0, 29.4, -17.7] | Gaborone |
| DRC | [12.0, -13.5, 31.4, 5.4] | Kinshasa |

### GHCN coverage (checked against NOAA's public archive)

PRCP days out of 731 in the window:

| Country | GHCN stations | With any PRCP | Best station |
|---|---|---|---|
| Mozambique | 19 | 9 | 265 |
| Botswana | 21 | 5 | 192 |
| Malawi | 20 | 2 | 182 |
| Zimbabwe | 20 | 3 | 129 |
| DRC | 13 | 5 | 41 |

**GHCN cannot carry the station QC on its own for any of these countries.**
The 2021–2023 window is no better (Zimbabwe has no GHCN rainfall at all).
Ground data needs to come from **TAHMO** or from **national met services**
via the `custom` upload.

## 4. Assumptions and limitations

* **TAHMO coverage for the five countries was not verified** from here. If
  Notebook 01 reports 0 stations, use `"custom"` with national data. GHCN is
  too sparse (see above).
* **Zambia is not in the requested list** although the workshop is in Zambia. It
  can be added the same way (bbox, capital, event) if needed.
* **DRC is very large**, and its equatorial north and centre have a bimodal rainy
  season with no single Nov–Mar season. The bbox pull covers the whole country
  (the export is still small at product resolution). The October 2023 heatwave
  mainly affected the southern DRC.
* **PERSIANN-CDR** (the TAMSAT fallback) is released with a delay of a few
  months. Check that it reaches June 2025 when you run it.
* **Notebook 04** was empty on `main`. The populated version from `tahmo-work`
  was brought in, and its crop list is unchanged (maize, sorghum, cassava and
  groundnut matter most in these countries).
* **Saved cell outputs** are from earlier East Africa runs. Re-run in Colab to
  regenerate them. Cached files go to `Datasets/<Country>/` on the Shared Drive.
