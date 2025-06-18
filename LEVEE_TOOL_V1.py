import pandas as pd
import rasterio 
import numpy as np
import geopandas as gpd 
import os 
from shapely.geometry import Point, MultiPoint, LineString
from scipy.spatial import cKDTree
import math 

Levee_filepath = r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\Levee\4310_MooneePonds_Ck_LeveeAlign.shp"
Lidar_filepath = r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\LiDAR\LiDAR_Merged_Tiles\LiDAR_Merge.tif"
flood_extent = r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_EXG_100y_657_max_h.asc"
waterway = r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\Waterway Alignment\4310_Channel_Centreline.shp"

def sample_lidar_values(points_gdf, Lidar_filepath):
    
    # get coordinates from the GeoDataFrame
    coords = [(geom.x, geom.y) for geom in points_gdf.geometry]
    # Open raster 
    with rasterio.open(Lidar_filepath) as src:
        # set an empty attribute 

        
        # Loop through each point and extract the raster value
        for index, point in points_gdf.iterrows():
            x, y = point.geometry.x, point.geometry.y
            row, col = src.index(x, y)  
            points_gdf.at[index, 'Elev(mAHD)'] = src.read(1)[row, col]



def sample_flood_extent(points_gdf, Lidar_filepath):
    coords = [(geom.x, geom.y) for geom in points_gdf.geometry]
    with rasterio.open(Lidar_filepath) as src:
        # this will need some processing to automate the column heading rather than hard code it 

        for index, point in points_gdf.iterrows():
            x, y = point.geometry.x, point.geometry.y
            row, col = src.index(x, y)  # Get the row and column of the raster
            points_gdf.at[index, '1% AEP'] = src.read(1)[row, col]


### STEP 1 #####
# READ in levee shapefile 
levee_shp = gpd.read_file(Levee_filepath)

# point on the levee every 10m? create a dataframe 
points_list = []
chainage_dataframe = gpd.GeoDataFrame({'Levee Identifier': [], 'Chainage': []})
# this is m 
distance = 10 

unique_levee_names = levee_shp['LOCATION_I'].unique()
looped_chain_count = 0
levee_count = 0

# this first step takes over 15 mins to complete lines to points
for line in levee_shp.geometry:
    # total length 
    length = line.length
    # create the points, length / total distance between points 

    points = int(length // distance)

    for point_loc in range(points + 1):  # +1 to include both endpoints
        # extract vertices of the line 
        # append to list 
        # pick start and end point of a line 
        point = line.interpolate(point_loc * distance) 
        # append to list 
        # Point is an attribute of shapely 
        points_list.append(Point(point.x, point.y)) 

        # updating chainage based on individual levee length
        unique_chainage = point_loc * 10
        chainage_dataframe.at[(looped_chain_count + point_loc), 'Chainage'] = unique_chainage
        chainage_dataframe.at[(looped_chain_count + point_loc), 'Levee Identifier'] = unique_levee_names[(levee_count)]
        
    looped_chain_count = looped_chain_count + point_loc + 1
    levee_count = levee_count + 1

# Create a new GeoDataFrame from the points
# crs same as input shapefile 
points_gdf = gpd.GeoDataFrame(geometry=points_list, crs=levee_shp.crs)

points_gdf = pd.concat([points_gdf, chainage_dataframe], axis = 1)
print("point sampling complete") 


# call the function to sample lidar 
lidar_values = sample_lidar_values(points_gdf, Lidar_filepath)

# sample flood results in the same location of the levee 
flood_extent_values = sample_flood_extent(points_gdf, flood_extent)

#### QUERY THIS DATA eg if this value is higher than the LiDAR -> attribute it 
for row in points_gdf:
    ## change this to -999 values 
    points_gdf['1% AEP'] = np.where(points_gdf['Elev(mAHD)'] < points_gdf['1% AEP'], points_gdf['1% AEP'], 0)


# filter where there is no data (0) in the dataframe 
zero_values =  points_gdf[points_gdf['1% AEP'] == 0]


# Load raster
with rasterio.open(flood_extent) as raster:
    raster_data = raster.read(1)
    transform = raster.transform
    nodata = raster.nodata

    # You can also build your valid raster coords/values if needed:
    rows, cols = np.indices(raster_data.shape)
    xs, ys = rasterio.transform.xy(transform, rows, cols, offset='center')
    xs_flat = np.array(xs).flatten()
    ys_flat = np.array(ys).flatten()
    values_flat = raster_data.flatten()
    valid_mask = values_flat != nodata
    valid_coords = np.column_stack((xs_flat[valid_mask], ys_flat[valid_mask]))
    valid_values = values_flat[valid_mask]



################## FLOOD EXTENT SEARCH ##########################

###  Not sure if this will make the code slow

# create a line (vector) from th point to the nearest waterway 
# sample raster along the line and store the value 
# identify the closest raster cell between the point and the waterway (exlcuding null values)


    closest_coords = []
    closest_values = []
    trace_lines = []

    # Extract the geometry column (GeoSeries)
    waterway_gdf = gpd.read_file(waterway)



    # check co or system
    if points_gdf.crs != waterway_gdf.crs:
        waterway_gdf = waterway_gdf.to_crs(points_gdf.crs)


    # Filter points where '1% AEP' == 0
    subset = points_gdf[points_gdf['1% AEP'] == 0]
    print("subset complete")

    # Extract just the geometries of the waterway
    waterway_lines = waterway_gdf.geometry

  

    for idx, point in subset.geometry.items():
        # Ensure valid Point - can probs delete this once checked? 
        if not isinstance(point, Point):
            point = Point(point)

        # Find the nearest waterway line
        distances = waterway_lines.distance(point)
        nearest_line_geom = waterway_lines.loc[distances.idxmin()]

        # Nearest point on the waterway line
        # project functin allows us to look along the waterway without converting waterway to points 
        nearest_point_on_water = nearest_line_geom.interpolate(nearest_line_geom.project(point))


        # this is forcing co ordinates to become float. 
        # something wasn't working with shapely, had valid geometries and points 
        p1 = (float(point.x), float(point.y))
        p2 = (float(nearest_point_on_water.x), float(nearest_point_on_water.y))
        trace_line = LineString([p1, p2])
        trace_lines.append(trace_line)
        


        # Sample along the trace line - edit from 100 to check 
        num_samples = 100
        sampled_points = [trace_line.interpolate(d, normalized=True) for d in np.linspace(0, 1, num_samples)]
        coords = [(p.x, p.y) for p in sampled_points]

        # Sample raster at those coordinates
        sampled_values = list(raster.sample(coords))
        sampled_values = [val[0] for val in sampled_values]

        # Get first non-nodata raster cell
        found = False
        for val, coord in zip(sampled_values, coords):
            # Handle band values
            v = val[0] if isinstance(val, (list, tuple, np.ndarray)) else val
            if v != raster.nodata:
                closest_coords.append(Point(coord))
                closest_values.append(float(v))
                found = True
                break

        if not found:
            closest_coords.append(None)
            closest_values.append(None)


    # Add to original GeoDataFrame 
    points_gdf.loc[subset.index, 'closest_cell_value'] = closest_values


for row in points_gdf:
    ## change this to -999 values 
    points_gdf['1% AEP'] = np.where(points_gdf['1% AEP'] == 0, points_gdf['closest_cell_value'], points_gdf['1% AEP'])

points_gdf = points_gdf.drop(columns=['closest_cell_value'])
points_gdf['1% AEP Diff'] = points_gdf['Elev(mAHD)'] - points_gdf['1% AEP']
print(points_gdf)

trace_lines_gdf = gpd.GeoDataFrame(geometry=trace_lines, crs=subset.crs)
trace_lines_gdf.to_file(r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\04 Scripting_Code\Python\Processing\trace_lines_4310_v3.shp")

code = str(levee_shp.iloc[0]['LOCATION_I'])
description = str(levee_shp.iloc[0]['DESCRIPTIO'])

code = code.replace('/','_')
# Save points to a new shapefile
#points_gdf.to_file(f"J:\\IE\\Projects\\03_Southern\\IA5000TK\\07 Technical\\04 Scripting_Code\\Python\\Processing\{code}.shp", driver='ESRI Shapefile')
points_gdf.to_file(f"J:\\IE\\Projects\\03_Southern\\IA5000TK\\07 Technical\\04 Scripting_Code\\Python\\Processing\{code}_check_TEST.shp", driver='ESRI Shapefile')