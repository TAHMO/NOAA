"""Smoke tests for the Conda/Jupyter environment defined in environment.yml."""

import importlib
import unittest


REQUIRED_MODULES = (
    "cartopy",
    "contextily",
    "ee",
    "folium",
    "geemap",
    "geopandas",
    "ipykernel",
    "ipywidgets",
    "matplotlib",
    "netCDF4",
    "numpy",
    "pandas",
    "rioxarray",
    "scipy",
    "seaborn",
    "shapely",
    "xarray",
)


class EnvironmentSmokeTests(unittest.TestCase):
    def test_required_modules_import(self):
        for module_name in REQUIRED_MODULES:
            with self.subTest(module=module_name):
                importlib.import_module(module_name)

    def test_notebook_helpers_import(self):
        from utils import CHIRPS_helpers, ERA5_helpers, IMERG_helpers, notebook_compat

        self.assertFalse(notebook_compat.is_colab())
        self.assertTrue(CHIRPS_helpers)
        self.assertTrue(ERA5_helpers)
        self.assertTrue(IMERG_helpers)


if __name__ == "__main__":
    unittest.main()
