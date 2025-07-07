import pandas as pd
import rasterio 
import numpy as np
import geopandas as gpd 
import os 
from shapely.geometry import Point, MultiPoint, LineString


Levee_filepath = r"J:\IE\Projects\03_Southern\IA339800\06 Technical\SurfaceWater\03 Spatial\Levee_Flood_Levels\Leeve_shp\SimsStreet.shp"
#Lidar_filepath = r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\LiDAR\DEM_Z_Engeny\TIFs\AMC_CC2100_100y_180m_tp29_011_MPC_DEM_Z.tif"
#waterway = r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\Waterway Alignment\4310_Channel_Centreline_GDA94.shp"
# please update these paths and labels to flood extent rasters using the same format as below. Comment out any extra filepaths.
flood_extent = {
    '1%_AEP' : r"J:\IE\Projects\03_Southern\IA339800\06 Technical\SurfaceWater\03 Spatial\Levee_Flood_Levels\LMAR_FloodExtents\LMAR_1AEP_T1AEP_Keilor_QT_base_076_h_Max.flt",
    '1%CC_AEP' : r"J:\IE\Projects\03_Southern\IA339800\06 Technical\SurfaceWater\03 Spatial\Levee_Flood_Levels\LMAR_FloodExtents\LMAR_1AEPCC_T1AEPCC_Keilor_QT_base_076_h_Max.flt",
    #'2%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_EXG_050y_657_max_h.asc",
    #'5%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_EXG_020y_657_max_h.asc",
    #'10%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_EXG_010y_657_max_h.asc",
    # '10%CC_AEP' : r'J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_CC_C_010y_657_max_h.asc',
    #'20%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_EXG_005y_657_max_h.asc",
    #'20%CC_AEP' : r'J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_CC_C_005y_657_max_h.asc'
 }


# Sample Flood Extent 
def sample_flood_extent(points_gdf, flood_extent, label):
    coords = [(geom.x
               , geom.y) for geom in points_gdf.geometry]
    with rasterio.open(raster_path) as src:
        for index, point in points_gdf.iterrows():
            x, y = point.geometry.x, point.geometry.y
            row, col = src.index(x, y)  
            # edit to pass label (1%AEP etc as set above). This function isn't doing anything until we call it further down. After we set label
            points_gdf.at[index, label] = src.read(1)[row, col]


########################### STEP 1 - READ IN LEVEE SHAPEFILE ###########################
# Read in levee shapefile 
levee_shp = gpd.read_file(Levee_filepath)

# point on the levee every 10m? create a dataframe 
points_list = []
chainage_dataframe = gpd.GeoDataFrame({'Levee Identifier': [], 'Chainage': []})
# this is m 
distance = 10 

unique_levee_names = levee_shp['ID'].unique()
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
#lidar_values = sample_lidar_values(points_gdf, Lidar_filepath)

# Sample flood extents 
# Loop through raster list of flood extents 
for label, raster_path in flood_extent.items():
    # sample flood results in the same location of the levee 
    flood_extent_values = sample_flood_extent(points_gdf, flood_extent, label)
    

print(points_gdf)


# Converting data columns to numeric to allow for rounding

    
############# STEP 5 - SAVE THE OUTPUT ###############        

code = str(levee_shp.iloc[0]['ID'])
name = str(levee_shp.iloc[0]['Name'])

code = code.replace('/','_')
# Save points to a new shapefile
#points_gdf.to_file(f"J:\\IE\\Projects\\03_Southern\\IA5000TK\\07 Technical\\04 Scripting_Code\\Python\\Processing\{code}.shp", driver='ESRI Shapefile')


points_gdf.to_file(f"J:\\IE\Projects\\03_Southern\\IA339800\\06 Technical\\SurfaceWater\\03 Spatial\\Levee_Flood_Levels\\LMAR_Levee_FloodLevels\\{code}_{name}_FloodLevel_1pct_1pctCC.shp", driver='ESRI Shapefile')