
'''
Datasets Required: (Stored at Shared Drives)
- IMERG ('2024-04-20','2024-05-05'), ('2022-07-28','2022-08-02'), ('2023-05-01','2023-05-08')
- TAMSAT ('2024-04-20','2024-05-05'), ('2022-07-28','2022-08-02'), ('2023-05-01','2023-05-08')
- ERA5 ('2024-04-20','2024-05-05'), ('2022-07-28','2022-08-02'), ('2023-05-01','2023-05-08')
- CHIRPS ('2024-04-20','2024-05-05'), ('2022-07-28','2022-08-02'), ('2023-05-01','2023-05-08')
- TAHMO Ground data ('2024-04-20','2024-05-05'), ('2022-07-28','2022-08-02'), ('2023-05-01','2023-05-08') #filtered to high confidence scores in kenya, uganda and Rwanda

'''
import ee
from utils.CHIRPS_helpers import get_chirps_pentad_gee
import json
import geopandas as gpd
import os
from utils.IMERG_helpers import get_imerg_raw
from utils.ERA5_helpers import era5_data_extracts, era5_var_handling
from utils.helpers import get_region_geojson
# from uti

# authenticate and initialize if not done during the runtime
try:
    ee.Initialize(project='noaa-tahmo')
except Exception as e:
    ee.Authenticate()
    ee.Initialize(project='noaa-tahmo')


with open('config.json') as f:
    config = json.load(f)

api_key = config['apiKey']
api_secret = config['apiSecret']
location_keys = config['location_keys']

def xmin_ymin_xmax_ymax(polygon):
    lons = [pt[0] for pt in polygon]
    lats = [pt[1] for pt in polygon]
    return min(lons), min(lats), max(lons), max(lats)

def fetch_region_google(query):
    """Primary: Fetch polygon geometry via Google Maps API"""
    if not location_key:
        raise RuntimeError("Missing Google Maps API key.")
    region_geom = get_region_geojson(query, location_key)['geometry']['coordinates'][0]
    return region_geom

def fetch_region_osm(query):
    """Fallback: Fetch geometry from OSM (Nominatim) via GeoPandas"""
    url = f"https://nominatim.openstreetmap.org/search?country={query}&format=geojson&polygon_geojson=1"
    gdf = gpd.read_file(url)
    if gdf.empty:
        raise ValueError("No OSM data found for that query.")
    geom = gdf.iloc[0].geometry
    if geom.geom_type == "Polygon":
        return list(geom.exterior.coords)
    elif geom.geom_type == "MultiPolygon":
        return list(list(geom.geoms)[0].exterior.coords)
    else:
        raise ValueError("Unsupported geometry type from OSM.")


# Extract CHIRPS in Kenya
def extract_product_data(start_date, end_date, product, region_query, output_location, variable=None):
  '''
  Accepted products: CHIRPS, TAMSAT, ERA5, IMERG, TAHMO
  '''
  # get the polygon for the region
  region_geom = fetch_region_google(region_query)
  xmin, ymin, xmax, ymax = xmin_ymin_xmax_ymax(region_geom)

  roi = ee.Geometry.Polygon(region_geom)

  # get the product

  # check to see if it is CHIRPS
  if product == 'CHIRPS':
    chirps_ds = get_chirps_pentad_gee(
        start_date=start_date,
        end_date=end_date,
        region=roi,
        export_path=f'{output_location}/CHIRPS_{start_date}_{end_date}_{region_query}.nc'
    )

    # Removing the imputed values resetting to NaN
    chirps_ds = chirps_ds.where(chirps_ds != -9999)
    chirps_ds.to_netcdf(f'{output_location}/CHIRPS_{region_query}.nc')
  elif product == 'ERA5':
    # ERA5 helper expects ISO-like datetime strings with time component (%Y-%m-%dT%H:%M:%S)
    iso_start_date = f"{start_date}T00:00:00"
    iso_end_date = f"{end_date}T23:59:59"

    era5_eac_data = era5_data_extracts(iso_start_date, iso_end_date, era5_l=False,
                                      polygon=region_geom)

    # te_ds_era5_ke = era5_var_handling(era5_ke_data, 'temperature_2m', xarray_ds=True)
    eac_pr_ds_era5 = era5_var_handling(era5_eac_data, 'total_precipitation', xarray_ds=True)

    # save to xarray
    eac_pr_ds_era5.to_netcdf(f'{output_location}/ERA5_precip_{start_date}_{end_date}_{region_query}.nc')

  elif product == 'IMERG':
    ds_raw = get_imerg_raw(start_date,
                        end_date,
                       region=roi,
                       export_path=f'{output_location}/IMERG_{start_date}_{end_date}_{region_query}.nc')

  elif product == 'TAMSAT':
    # !git clone https://github.com/TAMSAT/tamsat_download_extraction_api > /dev/null 2>&1

    file_path = "tamsat_download_extraction_api/tamsat_download_extract_api.py"

    # Read the file
    with open(file_path, "r") as file:
        content = file.read()

    # Replace np.float with np.float64
    content = content.replace("np.float", "np.float64")

    # Write back the modified content
    with open(file_path, "w") as file:
        file.write(content)

    # print("Successfully replaced np.float with np.float64.")

    import sys
    sys.path.append('tamsat_download_extraction_api')

    from tamsat_download_extract_api import download, extract
    import shutil

    # sys.path.append('utils/tamsat_download_extraction_api/tamsat_download_extract_api.py')

    # from utils.tamsat_download_extraction_api import download, extract
    curr_dir = os.getcwd()

    print(f'Beginning TAMSAT data extraction for {region_query}')
    # Download data
    download({
        "timestep": 'daily' ,
        "resolution": 0.0375,
        "start_date": start_date,
        "end_date": end_date,
        "version": 3.1,
        "localdata_dir": '/home/user/scripts/tamsat_api/data'
        })

    # check if the file exists
    if not os.path.exists('/home/user/scripts/tamsat_api/data'):
      print('❌ Errors occurred during extraction. Please restart the runtime and try again.')
    else:
      print('✅ TAMSAT data extraction completed successfully.')


    extract({
        "extract_type": 'domain',
        "S": ymin,
        "N": ymax,
        "W": xmin,
        "E": xmax,
        "timestep": 'daily',
        "start_date": start_date,
        "end_date": end_date,
        "version": 3.1,
        "localdata_dir": '/home/user/scripts/tamsat_api/data',
        "resolution": 0.0375
    })

    # save to drive
    file_path = f'/home/user/scripts/tamsat_api/data/extracted_data/domain/TAMSATv3.1_daily_0.0375_{ymax}_{ymin}_{xmin}_{xmax}_{start_date}_{end_date}.nc'

    # rename file path to TAMSAT_precip_{region_query}.nc
    shutil.move(file_path, f'{output_location}/TAMSAT_precip_{start_date}_{end_date}_{region_query}.nc')

  elif product == 'TAHMO':
    print("Extracting TAHMO data...")
    from utils.filter_stations import RetrieveData

    api_key = config['apiKey']
    api_secret = config['apiSecret']

    # Initialize the class
    rd = RetrieveData(apiKey=api_key,
                      apiSecret=api_secret)
    info = rd.get_stations_info()
    info = info[(info['location.longitude'] >= xmin) &
                            (info['location.longitude'] <= xmax) &
                            (info['location.latitude'] >= ymin) &
                            (info['location.latitude'] <= ymax)]
    region_precip = rd.multiple_measurements(stations_list=info['code'].tolist(),
                                     startDate=start_date,
                                     endDate=end_date,
                                     variables=['pr'],
                                         csv_file=f'{output_location}/tahmo_precip_{start_date}_{end_date}_{region_query}.csv')
    region_precip.to_csv(f'{output_location}/tahmo_precip_{start_date}_{end_date}_{region_query}.csv')
  else:
    raise ValueError(f"Product {product} not supported, either IMERG, TAMSAT, ERA5, CHIRPS or TAHMO")


# List out the events
EVENTS={'Kenya': ('2024-04-20','2024-05-05'),
        'Uganda': ('2022-07-28','2022-08-02'),
        'Rwanda': ('2023-05-01','2023-05-08')
        }

# extract for each country for all the products
products_list = ['ERA5', 'IMERG', 'TAHMO']
output_location = 'Datasets'

# check if file path exists
if not os.path.exists(output_location):
      os.makedirs(output_location, exist_ok=True)

# ------------------- Data Extraction Loop
for product in products_list:
  print(f"Extracting {product} data...")
  for country, dates in EVENTS.items():
    print(f"Extracting {product} data for {country} and dates {dates[0]} to {dates[1]}")
    extract_product_data(dates[0], dates[1], product, country, output_location)


