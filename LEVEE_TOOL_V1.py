import pandas as pd
import rasterio 
import numpy as np
import geopandas as gpd 
import os 
from shapely.geometry import Point, MultiPoint, LineString


# Define file paths, these will eventually be user inputs 
Levee_filepath = r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\0503\03 Exports\0503_Eumemmering_Ck_LeveeAlign_v2.shp"
Lidar_filepath = r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\0503\03 Exports\LiDAR_Merge.tif"
waterway = r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\0503\02 Build Data\03 Shapefiles\0503_Channel_Centreline.shp"
# please update these paths and labels to flood extent rasters using the same format as below. Comment out any extra filepaths.
flood_extent = {
    '1%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\0503\03 Exports\0503_EumemmeringCk_3m_aep1_h_Max.tif",
    # '1%CC_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_CC_C_100y_657_max_h.asc",
    '2%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\0503\03 Exports\0503_EumemmeringCk_3m_aep2_h_Max.tif",
    '5%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\0503\03 Exports\0503_EumemmeringCk_3m_aep5_h_Max.tif",
    '10%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\0503\03 Exports\0503_EumemmeringCk_3m_aep10_h_Max.tif",
    # '10%CC_AEP' : r'J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_CC_C_010y_657_max_h.asc',
    '20%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\0503\03 Exports\0503_EumemmeringCk_3m_aep20_h_Max.tif",
    #'20%CC_AEP' : r'J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_CC_C_005y_657_max_h.asc'
 
}

############# FUNCTIONS #####################
# SAMPLE LIDAR VALUES 
def sample_lidar_values(points_gdf, Lidar_filepath):
    # get coordinates from the GeoDataFrame
    coords = [(geom.x, geom.y) for geom in points_gdf.geometry]
    # Open raster 
    with rasterio.open(Lidar_filepath) as src:
        # Loop through each point and extract the raster value
        for index, point in points_gdf.iterrows():
            x, y = point.geometry.x, point.geometry.y
            row, col = src.index(x, y)  
            points_gdf.at[index, 'Elev(mAHD)'] = src.read(1)[row, col]



def sample_flood_extent(points_gdf, flood_extent, label):
    coords = [(geom.x
               , geom.y) for geom in points_gdf.geometry]
    with rasterio.open(raster_path) as src:
        for index, point in points_gdf.iterrows():
            x, y = point.geometry.x, point.geometry.y
            row, col = src.index(x, y)  
            # edit to pass label (1%AEP etc as set above). This function isn't doing anything until we call it further down. After we set label
            points_gdf.at[index, label] = src.read(1)[row, col]

##################################################

########################### STEP 1 - READ IN LEVEE SHAPEFILE ###########################
# Read in levee shapefile 
levee_shp = gpd.read_file(Levee_filepath)

# point on the levee every 10m? create a dataframe 
points_list = []
chainage_dataframe = gpd.GeoDataFrame({'Levee Identifier': [], 'Chainage': []})
# this is m 
distance = 10 

unique_levee_names = levee_shp['LOCATION_I'].unique()
looped_chain_count = 0
levee_count = 0

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



######################### STEP 2 - Call Functions for LiDAR and Flood Extent Samples ################################
# Sample lidar 
lidar_values = sample_lidar_values(points_gdf, Lidar_filepath)

# Sample flood extents 
# Loop through raster list of flood extents 
for label, raster_path in flood_extent.items():
    # sample flood results in the same location of the levee 
    flood_extent_values = sample_flood_extent(points_gdf, flood_extent, label)
    print(f"Flood extent sampling complete for {label}")

    #### QUERY THIS DATA eg if this value is higher than the LiDAR -> attribute it 
    for row in points_gdf:
        ## change this to -999 values 
        points_gdf[label] = np.where(points_gdf['Elev(mAHD)'] < points_gdf[label], points_gdf[label], 0)

    # filter where there is no data (0) in the dataframe 
    zero_values =  points_gdf[points_gdf[label] == 0]


################## STEP 3 FLOOD EXTENT SEARCH ##########################

# create a line (vector) from th point to the nearest waterway 
# sample raster along the line and store the value 
# identify the closest raster cell between the point and the waterway (exlcuding null values)

# Load raster
# Loop through raster list of flood extents 
for label, raster_path in flood_extent.items():
    with rasterio.open(raster_path) as raster:
        raster_data = raster.read(1)
        transform = raster.transform
        nodata = raster.nodata

        rows, cols = np.indices(raster_data.shape)
        xs, ys = rasterio.transform.xy(transform, rows, cols, offset='center')
        xs_flat = np.array(xs).flatten()
        ys_flat = np.array(ys).flatten()
        values_flat = raster_data.flatten()
        valid_mask = values_flat != nodata
        valid_coords = np.column_stack((xs_flat[valid_mask], ys_flat[valid_mask]))
        valid_values = values_flat[valid_mask]

      

        # Filter points where sampled value is 0 or nodata
        subset = points_gdf[(points_gdf[label] == 0)]
        print("subset complete")

        # create empty lists within the loop, so they are generated for each raster
        closest_coords = []
        trace_lines = []
        closest_values = []

        # Extract the geometry column (GeoSeries)
        waterway_gdf = gpd.read_file(waterway) 
        # check co or system
        if points_gdf.crs != waterway_gdf.crs:
            waterway_gdf = waterway_gdf.to_crs(points_gdf.crs)
        
        waterway_lines = waterway_gdf.geometry


        for idx, point in subset.geometry.items():
            if not isinstance(point, Point):
                point = Point(point)

            # Find nearest waterway line
            distances = waterway_lines.distance(point)
            nearest_line_geom = waterway_lines.loc[distances.idxmin()]

            # Nearest point on the waterway line
            # project functin allows us to look along the waterway without converting waterway to points 
            nearest_point_on_water = nearest_line_geom.interpolate(nearest_line_geom.project(point))


            # this is forcing co ordinates to become float. 
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
        points_gdf.loc[subset.index, f'{label}_closest_cell_value'] = closest_values

        # replace 0 values with the closest cell value 
        points_gdf[label] = np.where(points_gdf[label] == 0,
        points_gdf[f'{label}_closest_cell_value'], points_gdf[label])

        # we can now drop this column
        points_gdf = points_gdf.drop(columns=[f'{label}_closest_cell_value'])
        points_gdf[f'{label}_Diff'] = points_gdf['Elev(mAHD)'] - points_gdf[label]
        print(points_gdf)

        # for anaylsis of correct raster cells, output generated trace lines 
        #trace_lines_gdf = gpd.GeoDataFrame(geometry=trace_lines, crs=subset.crs)
        #trace_lines_gdf.to_file(f'J:\\IE\\Projects\\03_Southern\\IA5000TK\\07 Technical\\02 Hydraulics\\0503\\04 Code Processing\\{label}_trace_lines_0503_v2.shp')

############# STEP 4 - LoS Processing ##########################    

# Specifying the diff columns to search for the critical LoS

check_cols = []
numeric_cols = ["Elev(mAHD)"]

for label in flood_extent.keys():
    check_cols.append(f'{label}_Diff')
    numeric_cols.append(f'{label}')
    numeric_cols.append(f'{label}_Diff')

freeboard_events = points_gdf[check_cols].where(points_gdf[check_cols] > 0)
points_gdf['LoS'] = freeboard_events.idxmin(axis=1)
points_gdf.loc[freeboard_events.min(axis=1).isna(), 'LoS'] = ">20%"
points_gdf['LoS'] = points_gdf['LoS'].str.split('_').str[0]

# Converting data columns to numeric to allow for rounding

points_gdf[numeric_cols] = points_gdf[numeric_cols].apply(pd.to_numeric, errors='coerce').round(3)
    
############# STEP 5 - SAVE THE OUTPUT ###############        

code = str(levee_shp.iloc[0]['LOCATION_ID'])
description = str(levee_shp.iloc[0]['DESCRIPTIO'])

code = code.replace('/','_')
# Save points to a new shapefile
#points_gdf.to_file(f"J:\\IE\\Projects\\03_Southern\\IA5000TK\\07 Technical\\04 Scripting_Code\\Python\\Processing\{code}.shp", driver='ESRI Shapefile')
points_gdf = points_gdf.round(3)
points_gdf.to_file(f"J:\\IE\\Projects\\03_Southern\\IA5000TK\\07 Technical\\02 Hydraulics\\0503\\04 Code Processing\\{code}_Levee_LoS_Processing_v3.shp", driver='ESRI Shapefile')