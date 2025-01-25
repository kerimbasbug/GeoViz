import json
import pandas as pd
import geopandas as gpd
import os

import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

import folium
from shapely.geometry import Polygon, LineString
from folium.plugins import HeatMap, MarkerCluster

class CountyGeoData:
    """
    A class to analyze the relationship between parcel data and the closest wetlands 
    and Category 1 (C1) wetlands within a specific county.

    This class retrieves and processes geospatial data, including parcel data, wetland 
    data, and Category 1 wetland data, to determine spatial relationships. Due to the 
    computational complexity and time required for these analyses, the resulting relationships 
    are saved as a CSV file. The saved data can then be reused for various tasks by the 
    County class, avoiding the need for repeated calculations.
    
    Attributes:
        county_name (str): The name of the county.
        county_poly (Polygon): The polygon geometry representing the county boundary.
        parcel_df (GeoDataFrame): A GeoDataFrame containing parcel data for the county.
        wetland_df (GeoDataFrame): A GeoDataFrame containing wetland data for the county.
        c1_df (GeoDataFrame): A GeoDataFrame containing Category 1 wetland data for the county.
        geo_df (GeoDataFrame): A GeoDataFrame containing additional geographic data for the county.

    Methods:
        get_county_polygon(): Retrieves the polygon geometry of the county.
        get_parcel_data(): Retrieves parcel data for the county.
        get_wetland_data(): Retrieves wetland data for the county.
        get_c1_data(): Retrieves Category 1 wetland data for the county.
        get_geo_data(): Retrieves additional geographic data for the county.
        save_relationship_csv(): Saves the relationships between parcels, wetlands, 
                                 and C1 wetlands as a CSV file for later use.
    """
    def __init__(self, county_name):
        self.county_name = county_name
        self.county_poly = self.get_county_polygon()
        self.parcel_df = self.get_parcel_data()
        self.wetland_df = self.get_wetland_data()
        self.c1_df = self.get_c1_data()
        self.geo_df = self.get_geo_data()
    
    def get_county_polygon(self):
        """
        Load and return the polygon of the county from a GeoJSON file.

        Returns:
            shapely.geometry.Polygon: Polygon representing the county boundary.
        """
        with open('Data/CountyShapes/'+self.county_name.lower()+'.json', 'r') as file:
            data = json.load(file)
        county_coord = data[0]['geo_shape']['geometry']['coordinates'][0]
        county_poly = Polygon(county_coord)
        return county_poly

    def get_parcel_data(self):
        """
        Retrieves and processes parcel data for the specified county.

        Returns:
            GeoDataFrame: A GeoDataFrame containing parcel data for the county, with 
            columns standardized and geometry converted to EPSG:4326.
        """
        parcel_file = 'Data/CountyParcels/'+self.county_name.capitalize()+'County.gdb'
        parcel_df = gpd.read_file(parcel_file, layer=None)
        parcel_df['geometry'] = parcel_df['geometry'].to_crs("EPSG:4326")
        parcel_df = parcel_df.rename({'geometry':'PARCEL_GEOMETRY'}, axis=1)
        parcel_df.columns = parcel_df.columns.str.upper()
        parcel_df = parcel_df.reset_index(drop=True)
        return parcel_df
    
    def get_wetland_data(self):
        """
        Load and process wetland data for the county.

        Returns:
            geopandas.GeoDataFrame: GeoDataFrame containing wetland data within the county.
        """
        county_poly = self.get_county_polygon()
        wet_land_shapefile_path = 'Data/Wetlands_of_New_Jersey/Wetlands_of_New_Jersey_(from_Land_Use_Land_Cover_2012_Update).shp'
        wetland_df = gpd.read_file(wet_land_shapefile_path)
        if wetland_df.crs != "EPSG:4326":
            wetland_df = wetland_df.to_crs("EPSG:4326")
        wetland_df['is_in_county'] = wetland_df['geometry'].apply(lambda x: x.within(county_poly))
        wetland_df = wetland_df[wetland_df['is_in_county']]
        wetland_df = wetland_df.drop('is_in_county', axis=1)
        wetland_df = wetland_df.rename({'geometry':'WETLAND_GEOMETRY', 'OBJECTID': 'WETLAND_OBJECTID'}, axis=1)
        wetland_df.columns = wetland_df.columns.str.upper()
        wetland_df = wetland_df.reset_index(drop=True)
        return wetland_df

    def get_c1_data(self):
        """
        Load and process Category 1 waters data for the county.

        Returns:
            geopandas.GeoDataFrame: GeoDataFrame containing C1 waters data within the county.
        """
        county_poly = self.get_county_polygon()
        c1_waters_shapefile_path = 'Data/C1_Waters_of_New_Jersey/Category_One_(C1)_Waters_of_New_Jersey.shp'
        c1_df = gpd.read_file(c1_waters_shapefile_path)
        c1_df = c1_df.drop('DATE_EFFEC', axis=1)
        c1_df = gpd.GeoDataFrame(c1_df, geometry="geometry", crs="EPSG:3857")  # Assuming EPSG:3857
        c1_df = c1_df.to_crs("EPSG:4326") 
        c1_df['is_in_county'] = c1_df['geometry'].apply(lambda x: x.within(county_poly))
        c1_df = c1_df[c1_df['is_in_county']]
        c1_df = c1_df.drop('is_in_county', axis=1)
        c1_df = c1_df.rename({'geometry':'C1_geometry', 'OBJECTID': 'C1_OBJECTID'}, axis=1)
        c1_df.columns = c1_df.columns.str.upper()
        c1_df = c1_df.reset_index(drop=True)
        return c1_df

    def get_geo_data(self):
        """
        Load the GeoData file for the county.

        Returns:
            pandas.DataFrame: DataFrame containing geographic data for the county.
        """
        directory = "Data/CountyGeoData/" + self.county_name.capitalize() + "/"
        os.makedirs(directory, exist_ok=True)
        filename = self.county_name.capitalize()+"GeoData.csv"
        file_path = os.path.join(directory, filename)
        
        if os.path.isfile(file_path):
            print("GeoData for " + self.county_name.capitalize() + " County Exists!")
            df = pd.read_csv(directory+filename)
            return df
        else:
            print("Generating GeoData for " + self.county_name.capitalize() + " County...")
            columns = ['PAMS_PIN']
            parcels_gdf = gpd.GeoDataFrame(self.parcel_df, geometry="PARCEL_GEOMETRY", crs="EPSG:4326").to_crs(epsg=3857)
            parcels_gdf['PARCEL_CENTROID'] = parcels_gdf.apply(lambda x: x['PARCEL_GEOMETRY'].centroid, axis=1)
            
            if not self.wetland_df.empty:
                wetlands_gdf = gpd.GeoDataFrame(self.wetland_df, geometry="WETLAND_GEOMETRY", crs="EPSG:4326").to_crs(epsg=3857)
                wetlands_gdf['WETLAND_CENTROID'] = wetlands_gdf.apply(lambda x: x['WETLAND_GEOMETRY'].centroid, axis=1)
                distances_p_wl = parcels_gdf["PARCEL_CENTROID"].apply(lambda point: wetlands_gdf.distance(point))
                parcels_gdf["CLOSEST_WETLAND_ID"] = distances_p_wl.apply(lambda x: wetlands_gdf.loc[x.idxmin(), "WETLAND_OBJECTID"], axis=1)
                columns.append("CLOSEST_WETLAND_ID")

            if not self.c1_df.empty:
                c1_gdf = gpd.GeoDataFrame(self.c1_df, geometry="C1_GEOMETRY", crs="EPSG:4326").to_crs(epsg=3857)
                c1_gdf['C1_CENTROID'] = c1_gdf.apply(lambda x: x['C1_GEOMETRY'].centroid, axis=1)
                distances_p_c1 = parcels_gdf["PARCEL_CENTROID"].apply(lambda point: c1_gdf.distance(point))
                parcels_gdf["CLOSEST_C1_ID"] = distances_p_c1.apply(lambda x: c1_gdf.loc[x.idxmin(), "C1_OBJECTID"], axis=1)
                columns.append("CLOSEST_C1_ID")

            parcels_gdf = parcels_gdf[columns]
            parcels_gdf.to_csv(directory+filename, index=False)
            return parcels_gdf

if __name__ == '__main__':
    sussex = CountyGeoData(county_name='Sussex')