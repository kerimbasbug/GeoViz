import os
import json
import pandas as pd
import numpy as np
import geopandas as gpd

import folium
from folium.plugins import HeatMap, MarkerCluster
from shapely.geometry import Polygon, LineString, Point
from shapely.ops import nearest_points
from geopy.distance import geodesic

import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


class County:
    def __init__(self, county_name):
        """
        Initialize a County object with data related to the specified county.

        Parameters:
            county_name (str): Name of the county to analyze.
        """
        self.county_name = county_name
        self.county_poly = self.get_county_polygon()
        self.parcel_df = self.get_parcel_data()
        self.wetland_df = self.get_wetland_data()
        self.c1_df = self.get_c1_data()
        self.geo_df = self.get_geo_data() # Contains parcel-wetland and parcel-C1 relations.
        
        self.listing_df = self.get_listing_data()
        self.listing_df = self.fill_pams_pin() # Fills PAMS_PINs to identify the parcel

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

    def get_geo_data(self):
        """
        Load the GeoData file for the county.

        Returns:
            pandas.DataFrame: DataFrame containing geographic data for the county.
        """
        directory = "Data/CountyGeoData/" + self.county_name.capitalize() + "/"
        filename = self.county_name.capitalize()+"GeoData.csv"
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path):
            print('GeoData File Exists!')
            geo_df = pd.read_csv('Data/CountyGeoData/'+self.county_name.capitalize()+'/'+self.county_name.capitalize()+'GeoData.csv')
            return geo_df
        else:
            print('Please run CountyGeoData for ' + self.county_name + ' County!')
            return None

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

    def get_listing_data(self):
        """
        Load and process property listing data for the county.

        Returns:
            pandas.DataFrame: DataFrame containing processed property listing data.
        """
        listing_df = pd.read_csv('Data/Listings/listings.csv')
        listing_df = listing_df.dropna(subset=['PRICE'])
        listing_df = listing_df[listing_df['COUNTY'] == self.county_name]
        listing_df = listing_df.drop(['COUNTY','INTERESTED','FAVORITE','NEXT OPEN HOUSE END TIME','NEXT OPEN HOUSE START TIME'], axis=1)
        listing_df = listing_df.rename({'URL (SEE https://www.redfin.com/buy-a-home/comparative-market-analysis FOR INFO ON PRICING)':'URL'},axis=1)
        listing_df['LISTING_GEOMETRY'] = listing_df.apply(lambda x: Point(x['LONGITUDE'], x['LATITUDE']), axis=1)
        listing_df = listing_df.reset_index(drop=True)
        return listing_df
    
    def get_parcel_data(self):
        """
        Load and process parcel data for the county.

        Returns:
            geopandas.GeoDataFrame: GeoDataFrame containing parcel data.
        """
        parcel_file = 'Data/CountyParcels/'+self.county_name.capitalize()+'County.gdb'
        parcel_df = gpd.read_file(parcel_file, layer=None)
        parcel_df['geometry'] = parcel_df['geometry'].to_crs("EPSG:4326")
        parcel_df = parcel_df.rename({'geometry':'PARCEL_GEOMETRY'}, axis=1)
        parcel_df.columns = parcel_df.columns.str.upper()
        parcel_df = parcel_df.reset_index(drop=True)
        return parcel_df
        
    def fill_pams_pin(self):
        """
        Match property listings with parcels by spatial join and assign PAMS_PIN values.

        Returns:
            geopandas.GeoDataFrame: GeoDataFrame containing property listings with PAMS_PIN values.
        """
        parcels_gdf = gpd.GeoDataFrame(self.parcel_df, geometry="PARCEL_GEOMETRY", crs="EPSG:4326")
        listings_gdf = gpd.GeoDataFrame(self.listing_df, geometry="LISTING_GEOMETRY", crs="EPSG:4326")
        listings_with_pams_pin = gpd.sjoin(
            listings_gdf, parcels_gdf[["PAMS_PIN", "PARCEL_GEOMETRY"]], 
            how="left", 
            predicate="within"
        )
        listings_with_pams_pin.drop(columns=["index_right"], inplace=True)
        listings_with_pams_pin
        return listings_with_pams_pin
    
    def get_closest_wetland_c1(self, pams_pin):
        """
        Retrieve the closest wetland and C1 water geometry for a given parcel PAMS_PIN.

        Parameters:
            pams_pin (str): Unique identifier for the parcel.

        Returns:
            tuple: Closest wetland geometry and closest C1 water geometry.
        """
        closest_wl_id = self.geo_df[self.geo_df['PAMS_PIN'] == pams_pin]['CLOSEST_WETLAND_ID'].values[0]
        closest_wl_geo = self.wetland_df[self.wetland_df['WETLAND_OBJECTID'] == closest_wl_id]['WETLAND_GEOMETRY']
        
        closest_c1_id = self.geo_df[self.geo_df['PAMS_PIN'] == pams_pin]['CLOSEST_C1_ID'].values[0]
        closest_c1_geo = self.c1_df[self.c1_df['C1_OBJECTID'] == closest_c1_id]['C1_GEOMETRY']

        return closest_wl_geo, closest_c1_geo

    def show_county(self, save_to_html=False):
        """
        Generate a folium map visualizing the county with specific geographic layers and a heatmap of property listings.

        Parameters:
        -----------
        save_to_html : bool, optional
            If True, saves the map as an HTML file in the 'htmls/' directory with the county name in the filename.

        Returns:
        --------
        folium.Map
            A folium map object displaying the county, including:
            - Category 1 waters (in blue)
            - Wetlands (in red)
            - Property listings (as markers with links and tooltips)
            - A heatmap of property prices.
        """
        m = folium.Map(location=[self.county_poly.centroid.y, self.county_poly.centroid.x], zoom_start=10)
        folium.GeoJson(
                gpd.GeoDataFrame(self.c1_df, geometry='C1_GEOMETRY'),
                name="Category 1 Waters",
                style_function=lambda x: {
                    "color": "blue",
                    "weight": 2,
                    "opacity": 0.7,
                },
            ).add_to(m)

        folium.GeoJson(
                gpd.GeoDataFrame(self.wetland_df, geometry='WETLAND_GEOMETRY'),
                name="Wetlands",
                style_function=lambda x: {
                    "color": "red",
                    "weight": 2,
                    "opacity": 0.7,
                },
            ).add_to(m)

        heat_data = [[row['LATITUDE'], row['LONGITUDE'], row['PRICE']] for _, row in self.listing_df.iterrows()]
        heatmap_layer = folium.FeatureGroup(name='Heatmap Layer')

        HeatMap(heat_data, radius=10, blur=10, max_zoom=1).add_to(heatmap_layer)
        marker_cluster = MarkerCluster(name="Property Listings").add_to(m)

        heatmap_layer.add_to(m)

        for _, row in self.listing_df.iterrows():
            url_html = f'<a href="{row["URL"]}" target="_blank">Click for Details</a>'
            folium.Marker(
                location=[row['LATITUDE'], row['LONGITUDE']],
                popup=folium.Popup(url_html, max_width=300),  # Add the clickable HTML to the popup
                tooltip=f"Lat: {row['LATITUDE']}, Lon: {row['LONGITUDE']}",
            ).add_to(marker_cluster)

        folium.LayerControl().add_to(m)
        
        if save_to_html:
            m.save('htmls/'+self.county_name+'_map.html')
        
        return m
    
    def show_listing(self, ind=200, show_wl=True, show_c1=True, show_near_listings=True, circle_radius=1000, filters=None, save_to_html=False):
        """
        Generate a detailed folium map for a specific property listing, including surrounding features and optional filters.

        Parameters:
        -----------
        ind : int, optional
            Index of the property listing in the dataset to display. Defaults to 200.
        show_wl : bool, optional
            If True, display the nearest wetland and its distance to the property. Defaults to True.
        show_c1 : bool, optional
            If True, display the nearest Category 1 water and its distance to the property. Defaults to True.
        show_near_listings : bool, optional
            If True, display nearby listings within the defined circle radius. Defaults to True.
        circle_radius : int, optional
            Radius (in feet) to draw around the property for displaying nearby listings. Defaults to 1000 feet.
        filters : dict, optional
            A dictionary of filters to apply to nearby listings. Keys are column names, and values can be specific conditions (e.g., values or functions).
        save_to_html : bool, optional
            If True, saves the map as 'listing_view.html'. Defaults to False.

        Returns:
        --------
        tuple
            A tuple containing:
            - folium.Map: The generated map with layers and markers, including:
                - The selected property listing (in red)
                - The nearest wetland (if show_wl=True)
                - The nearest Category 1 water (if show_c1=True)
                - Nearby property listings (if show_near_listings=True)
                - A circle highlighting the area defined by `circle_radius`.
            - pandas.DataFrame: Filtered dataframe of nearby listings.
        """
        m = folium.Map(location=[self.listing_df.iloc[ind]['LISTING_GEOMETRY'].centroid.y, self.listing_df.iloc[ind]['LISTING_GEOMETRY'].centroid.x], zoom_start=15)
        parcel = self.parcel_df[self.parcel_df['PAMS_PIN'] == self.listing_df.iloc[ind]['PAMS_PIN']]['PARCEL_GEOMETRY']
        prcl = gpd.GeoDataFrame(geometry=[parcel.values[0]], crs='EPSG:4326')

        closest_wl, closest_c1 = self.get_closest_wetland_c1(self.listing_df.iloc[ind]['PAMS_PIN'])
        folium.GeoJson(parcel,
                        name="Parcel",
                        style_function=lambda x: {
                            "color": "blue",
                            "weight": 2,
                            "opacity": 0.7,
                        },
                    ).add_to(m)

        if show_wl:
            wl = gpd.GeoDataFrame(geometry=[closest_wl.values[0]], crs='EPSG:4326')
            point1, point2 = nearest_points(prcl.geometry[0], wl.geometry[0])
            closest_line_prcl_wl = LineString([point1, point2])
            gdf_closest_line = gpd.GeoDataFrame(geometry=[closest_line_prcl_wl], crs='EPSG:4326')#.to_crs(epsg=4326)
            closest_line_epsg4326 = gdf_closest_line.geometry[0]
            distance = geodesic(closest_line_epsg4326.coords[0], closest_line_epsg4326.coords[1]).feet
            folium.GeoJson(closest_wl,
                            name="Wetland",
                            style_function=lambda x: {
                                "color": "blue",
                                "weight": 2,
                                "opacity": 0.7,
                            },
                        ).add_to(m)
            folium.GeoJson(closest_line_prcl_wl,
                            name="Wetland Line",
                            style_function=lambda x: {
                                "color": "blue",
                                "weight": 2,
                                "opacity": 1,
                                "dashArray": "5, 5",
                            },
                            tooltip=f"Distance : {np.round(distance,2)} ft",
                        ).add_to(m)
        if show_c1:
            c1 = gpd.GeoDataFrame(geometry=[closest_c1.values[0]], crs='EPSG:4326')
            point3, point4 = nearest_points(prcl.geometry[0], c1.geometry[0])
            closest_line_prcl_c1 = LineString([point3, point4])
            gdf_closest_line = gpd.GeoDataFrame(geometry=[closest_line_prcl_c1], crs='EPSG:4326')#.to_crs(epsg=4326)
            closest_line_epsg4326 = gdf_closest_line.geometry[0]
            distance = geodesic(closest_line_epsg4326.coords[0], closest_line_epsg4326.coords[1]).feet
            folium.GeoJson(closest_c1,
                            name="C1",
                            style_function=lambda x: {
                                "color": "red",
                                "weight": 2,
                                "opacity": 0.7,
                            },
                        ).add_to(m)
            folium.GeoJson(closest_line_prcl_c1,
                            name="C1 Line",
                            style_function=lambda x: {
                                "color": "red",
                                "weight": 2,
                                "opacity": 1,
                                "dashArray": "5, 5",
                            },
                            tooltip=f"Distance : {np.round(distance,2)} ft",
                        ).add_to(m)


        data = {'LISTING_GEOMETRY': [self.listing_df.iloc[ind]['LISTING_GEOMETRY'].centroid]}  # Example point in Denver, CO
        gdf = gpd.GeoDataFrame(data, geometry='LISTING_GEOMETRY', crs="EPSG:4326")

        # Convert to a projected CRS (e.g., UTM Zone 13N for Denver)
        gdf_projected = gdf.to_crs(epsg=32613)  # UTM Zone 13N

        # Define radius in meters
        radius_meters = circle_radius*2 * 0.3048  # Convert feet to meters

        # Create the buffer in the projected CRS
        gdf_projected['BUFFER'] = gdf_projected['LISTING_GEOMETRY'].buffer(radius_meters)

        # Reproject the buffer back to EPSG:4326
        gdf['BUFFER'] = gdf_projected['BUFFER'].to_crs(epsg=4326)

        # Extract the circle geometry
        circle_geometry = gdf.iloc[0]['BUFFER']

        listings_in_circle = self.listing_df[self.listing_df['LISTING_GEOMETRY'].apply(lambda x: x.within(circle_geometry))]
        folium.GeoJson(circle_geometry,
                        name="Listing Circle",
                        style_function=lambda x: {'fillColor': 'green', 'color': 'green', 'fillOpacity': 0.2}
                        ).add_to(m)

        filtered_df = listings_in_circle.copy()
        if filters:
            for column, condition in filters.items():
                if callable(condition):
                    filtered_df = filtered_df[filtered_df[column].apply(condition) | pd.isna(filtered_df[column])]
                else:
                    #filtered_df = filtered_df[filtered_df[column] == condition]
                    filtered_df = filtered_df[(filtered_df[column] == condition) | pd.isna(filtered_df[column])]

            filtered_df = pd.concat([filtered_df, self.listing_df.iloc[[ind]]], ignore_index=True)
            filtered_df = filtered_df.drop_duplicates(subset=['URL'])

        if show_near_listings:
            for _, row in filtered_df.iterrows():
                if row.name != ind:
                    folium.Marker(
                        location=[row['LISTING_GEOMETRY'].centroid.y, row['LISTING_GEOMETRY'].centroid.x],
                        popup=folium.Popup(
                            f"""
                            <b>Price:</b> {row.get('PRICE', 'N/A')}<br>
                            <b>Beds:</b> {row.get('BEDS', 'N/A')}<br>
                            <b>Baths:</b> {row.get('BATHS', 'N/A')}<br>
                            <b>SQ Feet:</b> {row.get('SQUARE FEET', 'N/A')}<br>
                            <b>$/SQUARE FEET:</b> {row.get('$/SQUARE FEET', 'N/A')}<br>
                            <a href="{row.get('URL', '#')}" target="_blank">View Listing</a>
                            """,
                            max_width=300
                        ),
                        icon=folium.Icon(color="blue")
                    ).add_to(m)

        folium.Marker(
            location=[self.listing_df.iloc[ind]['LISTING_GEOMETRY'].centroid.y, self.listing_df.iloc[ind]['LISTING_GEOMETRY'].centroid.x],
            popup=folium.Popup(
                f"""
                <b>Price:</b> {self.listing_df.iloc[ind].get('PRICE', 'N/A')}<br>
                <b>Beds:</b> {self.listing_df.iloc[ind].get('BEDS', 'N/A')}<br>
                <b>Baths:</b> {self.listing_df.iloc[ind].get('BATHS', 'N/A')}<br>
                <b>SQ Feet:</b> {self.listing_df.iloc[ind].get('SQUARE FEET', 'N/A')}<br>
                <b>$/SQUARE FEET:</b> {self.listing_df.iloc[ind].get('$/SQUARE FEET', 'N/A')}<br>
                <a href="{self.listing_df.iloc[ind].get('URL', '#')}" target="_blank">View Listing</a>
                """,
                max_width=300
            ),
            icon=folium.Icon(color="red", icon="info-sign")
        ).add_to(m)

        folium.LayerControl().add_to(m)
        
        if save_to_html:
            m.save('listing_view.html')
        return m, filtered_df
    