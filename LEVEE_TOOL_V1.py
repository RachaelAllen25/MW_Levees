# Importing necessary python packages
# These will need to be installed on the users PC to run the script

import pandas as pd
import rasterio
import geopandas as gpd
from shapely.geometry import Point, LineString

# Define file paths
Levee_filepath = r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\Levee\4310_MooneePonds_Ck_LeveeAlign_GDA94_v2.shp"
Lidar_filepath = r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\LiDAR\DEM_Z_Engeny\TIFs\AMC_CC2100_100y_180m_tp29_011_MPC_DEM_Z.tif"
waterway = r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\Waterway Alignment\4310_Channel_Centreline_GDA94.shp"
 
# Define flood result file paths
# Flooding results loaded in with label function
# Comment out any extra filepaths/used AEP events.
flood_extent = {
    '1%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_EXG_100y_657_max_h.asc",
    '1%CC_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_CC_C_100y_657_max_h.asc",
    '2%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_EXG_050y_657_max_h.asc",
    '5%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_EXG_020y_657_max_h.asc",
    '10%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_EXG_010y_657_max_h.asc",
    # '10%CC_AEP' : r'J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_CC_C_010y_657_max_h.asc',
    '20%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_EXG_005y_657_max_h.asc",
    #'20%CC_AEP' : r'J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_CC_C_005y_657_max_h.asc'
 }

########## - FUNCTIONS - ##########

# SAMPLE LIDAR VALUES 
def sample_lidar_values(points_gdf, Lidar_filepath):
    
    # Get coordinates from the GeoDataFrame
    coords = [(geom.x, geom.y) for geom in points_gdf.geometry]
    
    # Open raster 
    with rasterio.open(Lidar_filepath) as src:
        
        # Loop through each point and extract the raster value
        # Saving sample value under 'Elev' field
        for index, point in points_gdf.iterrows():
            x, y = point.geometry.x, point.geometry.y
            row, col = src.index(x, y)  
            points_gdf.at[index, 'Elev(mAHD)'] = src.read(1)[row, col]

# SAMPLE FLOOD DATA
# The results from this function are not included in the code output
# Function was left in code for reference and potential repurposing at later development stages
def sample_flood_extent(points_gdf, flood_extent, label):
    
    # Get coordinates from the GeoDataFrame
    coords = [(geom.x, geom.y) for geom in points_gdf.geometry]
    
    # Open raster
    with rasterio.open(raster_path) as src:
        
        # Loop through each point and extract the raster value
        for index, point in points_gdf.iterrows():
            x, y = point.geometry.x, point.geometry.y
            row, col = src.index(x, y)

            # Label index used to assign raster sample to the right AEP event
            points_gdf.at[index, f'{label}_L'] = src.read(1)[row, col]

# LEVEE REORDERING
def extract_parts(code):

    # Returns the last entry of a sting
    suffix = code[-1]

    # int function is used to store data as an interger - this allows for numeric listing
    # .split function splits the string at the "LEV" input
    # [1] returns the second portion, i.e. after the LEV
    # [:-1] removes the last character (the suffix)
    number = int(code.split("LEV")[1][:-1])

    # Function print
    return suffix, number

########## - STEP 1 - READ IN LEVEE SHAPEFILE - ##########

# Read in levee shapefile 
levee_shp = gpd.read_file(Levee_filepath)

# Point on the levee - create an empty list 
points_list = []
chainage_dataframe = gpd.GeoDataFrame({'Levee Identifier': [], 'Chainage': []})

# This is in metres
distance = 10 

# Create unique list of levee names based on MW data field
code_levee_names = levee_shp['LOCATION_I'].unique()

# Empty dictionary to store reordered data
levee_order = {}

# For loop to cycle through each unique name
for code in code_levee_names:
        
    # Earlier function is called to generate the levee name parts
    suffix, number = extract_parts(code)

    # This line is grouping the levee names into their respective parts
    # By using .sefdefault, if suffix is not already a key in the dictionary, a new key (suffix) and empty list pairing is added
    # If suffix is already a key, the existing value associated with that key is returned (i.e. no change) and the new number/code pairing is appended
    # Assigning the number as the first dictionary input is input to the next line of code
    levee_order.setdefault(suffix, []).append((number, code))

# This is the equivalent to a nested for loop
# First the suffix (dictionary key) are iterated over - using sorted(levee_order) groups the E, W, etc. alphabetically
# The second for loop iterates over the levee code based on the number entry, e.g. 1, 2 etc.
# The _, represents a variable that is being iterated over in the for loop but is not used
# print statement to show reordered codes
reordered_codes = [code for suffix in sorted(levee_order) for _, code in sorted(levee_order[suffix])]
print(reordered_codes)

# Reindexing the levee shp based on new index order
# print statement to show reordered dataframe
levee_shp = levee_shp.set_index('LOCATION_I').loc[reordered_codes]
levee_shp.reset_index(inplace = True)
print(levee_shp)

# Initialising datasets for chainage reordering
unique_levee_names = levee_shp['LOCATION_I'].unique()
looped_chain_count = 0
levee_count = 0

# For loop to work through each line feature in the input shapefile - each specific levee asset
for line in levee_shp.geometry:
    
    # Total length 
    length = line.length
    
    # Create the points, length / total distance between points 
    points = int(length // distance)

    # Indented for loop to work through each point generated from the one unique line feature
    for point_loc in range(points + 1): # +1 to include both endpoints
       
        # Extract vertices of the line 
        # Append to list 
        # Pick start and end point of a line 
        point = line.interpolate(point_loc * distance) 
       
        # Append to list 
        # Point is an attribute of shapely 
        points_list.append(Point(point.x, point.y)) 

        # Updating chainage based on individual levee length
        unique_chainage = point_loc * 10
        chainage_dataframe.at[(looped_chain_count + point_loc), 'Chainage'] = unique_chainage
        chainage_dataframe.at[(looped_chain_count + point_loc), 'Levee Identifier'] = unique_levee_names[(levee_count)]

    # Resetting loop counter for next individual levee asset
    looped_chain_count = looped_chain_count + point_loc + 1
    levee_count = levee_count + 1

# Create a new GeoDataFrame from the points
# Crs same as input shapefile 
points_gdf = gpd.GeoDataFrame(geometry=points_list, crs=levee_shp.crs)

# Connecting points and chainage dataframe
points_gdf = pd.concat([points_gdf, chainage_dataframe], axis = 1)
print("point sampling complete") 

########## - STEP 2 - CALL FUNCTIONS - ##########

# Sample lidar 
lidar_values = sample_lidar_values(points_gdf, Lidar_filepath)

# Sample flood extents 
# Loop through raster list of flood extents 
for label, raster_path in flood_extent.items():
    
    # Sample flood results in the same location of the levee 
    flood_extent_values = sample_flood_extent(points_gdf, flood_extent, label)
    print(f"Flood extent sampling complete for {label}")

########## - STEP 3 - FLOOD EXTENT SEARCH - ##########

# create a line (vector) from th point to the nearest waterway 
# sample raster along the line and store the value 
# identify the closest raster cell between the point and the waterway (exlcuding null values)

# Load raster
# Loop through raster list of flood extents 
for label, raster_path in flood_extent.items():
    with rasterio.open(raster_path) as raster:
       
        # Create empty lists within the loop, so they are generated for each raster
        trace_lines = []
        central_waterway_values = []

        # Extract the geometry column (GeoSeries)
        # Reading the provided waterway shapefile into a geodataframe - this includes all of the associated shapefile attribute data
        waterway_gdf = gpd.read_file(waterway) 
        
        # check co or system
        if points_gdf.crs != waterway_gdf.crs:
            waterway_gdf = waterway_gdf.to_crs(points_gdf.crs)
        
        # This command strips the attribute data from the geodataframe and leaves only the geometries
        waterway_lines = waterway_gdf.geometry

        # This for loop interates over both the index and the geometry (points) with the geodata frame
        # Using the .geometry isolates the geometry
        # Using .items allow iteration over both idx and points
        for idx, point in points_gdf.geometry.items():
            
            # This check assess whether each geometry (point) is a Point class
            if not isinstance(point, Point):
                # Converts to Point if not
                point = Point(point)

            # Find nearest waterway line
            # These lines of code only find the nearest individual line feature (as provided in the shapefile)
            distances = waterway_lines.distance(point)
            nearest_line_geom = waterway_lines.loc[distances.idxmin()]

            # Nearest point on the waterway line
            # project functin allows us to look along the waterway without converting waterway to points
            # The .project function returns the distance along the nearest_line_geom that is closest to the levee reference point
            # The .interpolate function converts the distance output into a coordinate reference 
            nearest_point_on_water = nearest_line_geom.interpolate(nearest_line_geom.project(point))

            # This is forcing co ordinates to become float
            # Completed for both the levee reference point and the nearest coordinate point on the waterway line
            p1 = (float(point.x), float(point.y))
            p2 = (float(nearest_point_on_water.x), float(nearest_point_on_water.y))

            # Creation of line string between the two points
            trace_line = LineString([p1, p2])
            trace_lines.append(trace_line)

            centre_sample_value = float(next(raster.sample([(nearest_point_on_water.x, nearest_point_on_water.y)]))[0])
            central_waterway_values.append(centre_sample_value)

        # Add to original GeoDataFrame 
        points_gdf[f'{label}_W'] = central_waterway_values

        # replace 0 values with the closest cell value 
        # points_gdf[label] = np.where(points_gdf[label] == 0,
        # points_gdf[f'{label}_Strm'], points_gdf[label])

        # we can now drop this column
        # points_gdf = points_gdf.drop(columns=[f'{label}_central_cell_value'])
        points_gdf[f'{label}_D'] = points_gdf['Elev(mAHD)'] - points_gdf[f'{label}_W']
        print(points_gdf)

        # for anaylsis of correct raster cells, output generated trace lines 
        trace_lines_gdf = gpd.GeoDataFrame(geometry=trace_lines, crs=points_gdf.crs)
        trace_lines_gdf.to_file(f'J:\\IE\\Projects\\03_Southern\\IA5000TK\\07 Technical\\02 Hydraulics\\0503\\04 Code Processing\\{label}_trace_lines_central_sample.shp')

########## - STEP 4 - LOS PROCESSING - ##########

# Specifying the diff columns to search for the critical LoS
# The [] command specifies an empty python list
# A python list is one-dimensional, i.e. just a row of values
check_cols = []
numeric_cols = ["Elev(mAHD)"]

# For loop to group numeric datasets
for label in flood_extent.keys():

    # Check cols includes _Diff values for subsequent LoS testing
    check_cols.append(f'{label}_D')
    numeric_cols.append(f'{label}_L')
    numeric_cols.append(f'{label}_W')
    numeric_cols.append(f'{label}_D')

# Creating a subset of the points_gdf from the check_cols
# This is using the list as an index reference for points_gdf
# All other values not > 0 become NaN
freeboard_events = points_gdf[check_cols].where(points_gdf[check_cols] > 0)

# Taking the LoS classification corresponding to the minimum freeboard event
# idxmin is used rather than min to return the the label/position of the column rather than the actual data entry.
# axis = 1 locks to search to be per row
# idxmin does not search and return on NaN values
points_gdf['LoS'] = freeboard_events.idxmin(axis=1)

# Code to assign the levee points with no observed freeboard a LoS rank of >20%
# .min(axis=1) is used to search for the minimum value along each row of the freeboard_events gdf
# .isna() checks if the minimum row value is NaN. Based on the last row of code, if NaN is found, the whole row will be NaN
# Result is a series of true/false - one per row
# If .isna() condition met, >20% assiged to the LoS column
points_gdf.loc[freeboard_events.min(axis=1).isna(), 'LoS'] = ">20%"

# Dropping the _Diff suffix from LoS entry
points_gdf['LoS'] = points_gdf['LoS'].str.split('_').str[0]

# Converting data columns to numeric to allow for rounding
points_gdf[numeric_cols] = points_gdf[numeric_cols].apply(pd.to_numeric, errors='coerce').round(3)
    
########## - STEP 5 - SAVE THE OUTPUT - ##########

# Extracting asset name to use in file export
code = str(levee_shp.iloc[0]['LOCATION_I'])
code = code.split("/")[0]

# Save points to a new shapefile
points_gdf = points_gdf.round(3)
points_gdf.to_file(f"J:\\IE\\Projects\\03_Southern\\IA5000TK\\07 Technical\\04 Scripting_Code\\Python\\Processing\\{code}_Levee_LoS_Processing_Central_Sample_v2.shp", driver='ESRI Shapefile')