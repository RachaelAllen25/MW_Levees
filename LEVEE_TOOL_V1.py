# Importing necessary python packages
# These will need to be installed on the users PC to run the script
import os
import pandas as pd
import rasterio
import numpy as np
import geopandas as gpd
from shapely.geometry import Point, LineString

########## - STEP 1 - DEFINE FUNCTIONS - ##########

def check_file(path, expected_ext=None):
    """
    Validates the existence and file extension of a given file path.

    Parameters:
    ----------
    path : str
        The full file path to check.

    expected_ext : str, optional
        The expected file extension (e.g., '.tif', '.shp'). If provided, the function will verify
        that the file has this extension.

    Raises:
    -------
    FileNotFoundError
        If the file does not exist at the specified path.

    ValueError
        If the file exists but does not match the expected extension.

    Notes:
    ------
    - This function is useful for pre-validating inputs before attempting file I/O operations.
    - The extension check is case-insensitive.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")
    if expected_ext and not path.lower().endswith(expected_ext):
        raise ValueError(f"Unexpected file format for {path}. Expected {expected_ext}")

def get_crs(filepath):
    """
    Determines the coordinate reference system (CRS) of a spatial file.

    Supports vector files (.shp) and raster files (.tif, .asc, .flt).

    Parameters:
    ----------
    filepath : str
        The full path to the spatial file whose CRS is to be retrieved.

    Returns:
    --------
    CRS object
        The CRS of the file, as returned by GeoPandas (for shapefiles) or Rasterio (for raster files).

    Raises:
    -------
    ValueError
        If the file extension is not supported for CRS checking.

    RuntimeError
        If the file cannot be read or the CRS cannot be determined.

    Notes:
    ------
    - For raster formats like .asc and .flt, ensure accompanying metadata files (.prj, .hdr) are present
      to allow proper CRS detection.
    - Result files are assumed to be in a .tif, .asc or .flt format.
    """
    try:
        if filepath.lower().endswith('.shp'):
            gdf = gpd.read_file(filepath)
            return gdf.crs
        elif filepath.lower().endswith(('.tif', '.asc', '.flt')):
            with rasterio.open(filepath) as src:
                return src.crs
        else:
            raise ValueError(f"Unsupported file type for CRS check: {filepath}")
    except Exception as e:
        raise RuntimeError(f"Failed to read CRS from {filepath}: {e}")

def sample_lidar_values(points_gdf, Lidar_filepath):
    """
    Samples elevation values from a LiDAR raster and assigns them to point geometries in a GeoDataFrame.

    Parameters:
    ----------
    points_gdf : geopandas.GeoDataFrame
        A GeoDataFrame containing point geometries. This object will be modified in place by adding a new column
        'Elev(mAHD)' with elevation values sampled from the raster.

    Lidar_filepath : str
        File path to the LiDAR raster (e.g., a .tif file) from which elevation values will be extracted.

    Notes:
    ------
    - This function modifies `points_gdf` in place. To preserve the original data, consider passing a copy.
    - Elevation values are sampled using the raster index corresponding to each point's coordinates.
    - The function assumes the raster and point geometries are in the same coordinate reference system (CRS).

    Returns:
    --------
    None
    """
    
    # Open raster
    with rasterio.open(Lidar_filepath) as src:
        
        # Loop through each point and extract the raster value
        for index, point in points_gdf.iterrows():
            x, y = point.geometry.x, point.geometry.y
            row, col = src.index(x, y)
            points_gdf.at[index, 'Elev(mAHD)'] = src.read(1)[row, col]

def sample_flood_extent(points_gdf, raster_path, label):
    """
    Samples flood extent values from a raster and assigns them to point geometries in a GeoDataFrame.

    Parameters:
    ----------
    points_gdf : geopandas.GeoDataFrame
        A GeoDataFrame containing point geometries. This object will be modified in place by adding a new column
        named '{label}_L' with values sampled from the raster.

    raster_path : str
        File path to the flood extent raster (e.g., a .tif file) from which values will be extracted.

    label : str
        A label used to name the new column in the GeoDataFrame. For example, if label is 'AEP_1%', the new column
        will be 'AEP_1%_L'.

    Notes:
    ------
    - This function modifies `points_gdf` in place. To preserve the original data, consider passing a copy.
    - The `flood_extent` parameter is unused and should be removed unless needed for future logic.
    - The function assumes the raster and point geometries are in the same coordinate reference system (CRS).
    - Consider refactoring with `sample_lidar_values` into a generic raster sampling function to reduce duplication.

    Returns:
    --------
    None
    """
    # Open raster
    with rasterio.open(raster_path) as src:
        
        # Loop through each point and extract the raster value
        for index, point in points_gdf.iterrows():
            x, y = point.geometry.x, point.geometry.y
            row, col = src.index(x, y)

            # Assign raster value to a column named with the label
            points_gdf.at[index, f'{label}_L'] = src.read(1)[row, col]

def levee_identifier(levee_code):
    """
    Extracts the suffix and numeric identifier from a levee code string.

    Parameters:
    ----------
    levee_code : str
        A string representing a levee identifier, expected to follow the format 'LEV<number><suffix>'.
        Example: 'LEV12A' → number: 12, suffix: 'A'

    Returns:
    --------
    tuple
        A tuple containing:
        - suffix (str): The last character of the levee code.
        - number (int): The numeric portion of the levee code, extracted from between 'LEV' and the suffix.

    Notes:
    ------
    - The function assumes the levee code always starts with 'LEV' and ends with a single-character suffix.
    - If the format varies, consider adding error handling to validate the input.
    """
    suffix = levee_code[-1]
    number = int(levee_code.split("LEV")[1][:-1])
    return suffix, number

def trace_and_sample(point, waterway_lines, raster, nodata_value=-9999):
    """
    Traces a point to the nearest waterway line and samples a raster at the nearest location.

    Parameters:
        point (shapely.geometry.Point): The reference point.
        waterway_lines (GeoSeries): Geometry of waterway lines.
        raster (rasterio.DatasetReader): Open raster object.
        nodata_value (float): Value to return if raster sample is missing or invalid.

    Returns:
        tuple: (trace_line: LineString, sampled_value: float)
    """
    if not isinstance(point, Point):
        point = Point(point)

    # Find nearest waterway line
    distances = waterway_lines.distance(point)
    nearest_line_geom = waterway_lines.loc[distances.idxmin()]

    # Nearest point on the waterway line
    nearest_point_on_water = nearest_line_geom.interpolate(
        nearest_line_geom.project(point)
    )

    # Create trace line
    p1 = (float(point.x), float(point.y))
    p2 = (float(nearest_point_on_water.x), float(nearest_point_on_water.y))
    trace_line = LineString([p1, p2])

    # Sample raster with nodata handling
    sample = next(raster.sample([(nearest_point_on_water.x, nearest_point_on_water.y)]))[0]
    sampled_value = float(sample) if sample is not None and not np.isnan(sample) else nodata_value

    return trace_line, sampled_value

########## - STEP 2 - DEFINE FILEPATHS & CALL DATA CHECK FUNCTIONS - ##########

# Define file paths
Levee_filepath = r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\04 Scripting_Code\Python\Review\Levee_Align\4420_Merri_Ck_LeveeAlign_v2.shp"
Lidar_filepath = r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\04 Scripting_Code\Python\Review\LiDAR\LiDAR_Merge.tif"
waterway = r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\04 Scripting_Code\Python\Review\Waterway_Align\4420_Channel_Centreline.shp"

# Extract lidar extension for input into check_file function
# Levee and waterway files assumed to be in .shp format
_, lidar_ext = os.path.splitext(Lidar_filepath)

# Define flood result file paths
# Flooding results loaded in with label function - this creates a dictionary of titles and filepaths
# Comment out any extra filepaths/used AEP events
flood_extent = {
    '1%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\04 Scripting_Code\Python\Review\Flood_Results\4420_MerriCreek_3m_aep1_h_Max.tif",
    # '1%CC_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_CC_C_100y_657_max_h.asc",
    '2%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\04 Scripting_Code\Python\Review\Flood_Results\4420_MerriCreek_3m_aep2_h_Max.tif",
    '5%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\04 Scripting_Code\Python\Review\Flood_Results\4420_MerriCreek_3m_aep5_h_Max.tif",
    '10%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\04 Scripting_Code\Python\Review\Flood_Results\4420_MerriCreek_3m_aep10_h_Max.tif",
    # '10%CC_AEP' : r'J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_CC_C_010y_657_max_h.asc',
    '20%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\04 Scripting_Code\Python\Review\Flood_Results\4420_MerriCreek_3m_aep20_h_Max.tif",
    #'20%CC_AEP' : r'J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_CC_C_005y_657_max_h.asc'
 }

# Validate input files
check_file(Levee_filepath, '.shp')
check_file(Lidar_filepath, lidar_ext)
check_file(waterway, '.shp')

# Printing data projections
levee_crs = get_crs(Levee_filepath)
print(f"Levee is projected in {levee_crs}")
lidar_crs = get_crs(Lidar_filepath)
print(f"LiDAR is projected in {lidar_crs}")
waterway_crs = get_crs(waterway)
print(f"Waterway is projected in {waterway_crs}")

# Check CRS compatibility
if lidar_crs != levee_crs:
    print("⚠️ Warning: Levee shapefile and LiDAR raster have different coordinate systems")
if waterway_crs != levee_crs:
    print("⚠️ Warning: Levee shapefile and waterway shapefile have different coordinate systems")

for label, path in flood_extent.items():
    check_file(path, '.tif')
    flood_crs = get_crs(path)
    print(f"The {label} result is projected in {flood_crs}")
    if flood_crs != levee_crs:
        print(f"⚠️ Warning: CRS mismatch between levee shapefile and flood extent raster for {label}")

########## - STEP 3 - READ IN LEVEE SHAPEFILE & REFORMAT - ##########

# Read in levee shapefile 
levee_shp = gpd.read_file(Levee_filepath) 

# Point on the levee - create an empty list 
points_list = []
chainage_dataframe = gpd.GeoDataFrame({'Levee Identifier': [], 'Chainage': []})

# This is in metres
distance = 10 

# Create unique list of levee names based on MW data field
code_levee_names = levee_shp['LOCATION_I'].unique() #NF: This line assumes that LOCATION_I is the field containing the levee names. It that always the case? If not, consider making it a parameter of the function or checking the field name before using it.

# Empty dictionary to store reordered data
levee_order = {}

# For loop to cycle through each unique name
for asset_code in code_levee_names: #NF: code is a bit generic. Consider other name. 
    #NF: Consider not using the same variable name as the function parameter in the loop, as it can lead to confusion and potentially overwrite the function parameter.
        
    # Earlier function is called to generate the levee name parts
    suffix, number = levee_identifier(asset_code)

    # This line is grouping the levee names into their respective parts
    # By using .sefdefault, if suffix is not already a key in the dictionary, a new key (suffix) and empty list pairing is added
    # If suffix is already a key, the existing value associated with that key is returned (i.e. no change) and the new number/code pairing is appended
    # Assigning the number as the first dictionary input is input to the next line of code
    levee_order.setdefault(suffix, []).append((number, asset_code)) #NF: not sure I understand the reason why the levees are split by suffix?

# This is the equivalent to a nested for loop
# First the suffix (dictionary key) are iterated over - using sorted(levee_order) groups the E, W, etc. alphabetically
# The second for loop iterates over the levee code based on the number entry, e.g. 1, 2 etc.
# The _, represents a variable that is being iterated over in the for loop but is not used
# print statement to show reordered codes
reordered_levee_codes = [code for suffix in sorted(levee_order) for _, code in sorted(levee_order[suffix])] #NF: what is the purpose of creating this sorted list? 
print(reordered_levee_codes) #NF: This print statement is useful for debugging, but consider removing it or replacing it with a logging statement in production code.

# Reindexing the levee shp based on new index order
# print statement to show reordered dataframe
levee_shp_reorder = levee_shp.copy()
levee_shp_reorder = levee_shp_reorder.set_index('LOCATION_I').loc[reordered_levee_codes] #NF: I would usualy to this on a copy of the dataframe to avoid modifying the original data.
levee_shp_reorder.reset_index(inplace = True)
# print(levee_shp_reorder) #NF: This print statement is useful for debugging, but consider removing it or replacing it with a logging statement in production code.

# Initialising datasets for chainage reordering
# unique_levee_names = levee_shp['LOCATION_I'].unique() #NF: This is the same as code_levee_names, consider removing one of them to avoid redundancy.
looped_chain_count = 0
levee_count = 0

# For loop to work through each line feature in the input shapefile - each specific levee asset
for line in levee_shp_reorder.geometry:
    
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
        chainage_dataframe.at[(looped_chain_count + point_loc), 'Levee Identifier'] = code_levee_names[(levee_count)]

        # If the point_loc is the last point, we can add the end point of the levee line using
        if point_loc == points:
             # Append the last point of the levee line
             end_point = line.interpolate(length)
             points_list.append(Point(end_point.x, end_point.y))
            
             # Update chainage for the end point
             unique_chainage = length
             chainage_dataframe.at[(looped_chain_count + point_loc + 1), 'Chainage'] = unique_chainage
             chainage_dataframe.at[(looped_chain_count + point_loc + 1), 'Levee Identifier'] = code_levee_names[(levee_count)]
        
    # Resetting loop counter for next individual levee asset
    looped_chain_count = looped_chain_count + point_loc + 2
    levee_count = levee_count + 1

# Create a new GeoDataFrame from the points
# Crs same as input shapefile 
points_gdf = gpd.GeoDataFrame(geometry=points_list, crs=levee_shp_reorder.crs)

# Connecting points and chainage dataframe
points_gdf = pd.concat([points_gdf, chainage_dataframe], axis = 1)
print("Point delineation complete") 

########## - STEP 4 - CALL SAMPLE FUNCTIONS - ##########

# Sample lidar 
lidar_values = sample_lidar_values(points_gdf, Lidar_filepath)
print("LiDAR point sampling complete")

# Sample flood extents 
# Loop through raster list of flood extents 
for label, raster_path in flood_extent.items():
    
    # Sample flood results in the same location of the levee 
    flood_extent_values = sample_flood_extent(points_gdf, raster_path, label) # NF: See comments for function
    print(f"Flood extent sampling complete for {label}")

########## - STEP 5 - FLOOD EXTENT SEARCH - ##########

# create a line (vector) from the levee point to the nearest waterway 
# sample raster along the line and store the value 
# identify the closest raster cell between the point and the waterway (excluding null values)

# Load raster
# Loop through raster list of flood extents 

waterway_gdf = gpd.read_file(waterway)
if points_gdf.crs != waterway_gdf.crs:
    waterway_gdf = waterway_gdf.to_crs(points_gdf.crs) 

for label, raster_path in flood_extent.items():
    with rasterio.open(raster_path) as raster:
       
        # Create empty lists within the loop, so they are generated for each raster
        trace_lines = []
        central_waterway_values = []

        # check coordinate system
        if points_gdf.crs != waterway_gdf.crs:
            waterway_gdf = waterway_gdf.to_crs(points_gdf.crs) 
        
        # This command strips the attribute data from the geodataframe and leaves only the geometries
        waterway_lines = waterway_gdf.geometry

        # This for loop interates over both the index and the geometry (points) with the geodata frame
        # Using the .geometry isolates the geometry
        # Using .items allow iteration over both idx and points
        for idx, point in points_gdf.geometry.items():
            trace_line, value = trace_and_sample(point, waterway_gdf.geometry, raster)
            trace_lines.append(trace_line)
            central_waterway_values.append(value)
            
        # Add to original GeoDataFrame 
        points_gdf[f'{label}_W'] = central_waterway_values

        # we can now drop this column NF: no column is dropped in the current code as it is commented out
        points_gdf[f'{label}_D'] = points_gdf['Elev(mAHD)'] - points_gdf[f'{label}_W'] #NF: this only compared the water level in the river and the elevation. What about the water level at the levee?
        print(points_gdf) #NF: This print statement is useful for debugging, but consider removing it or replacing it with a logging statement in production code.

# for anaylsis of correct raster cells, output generated trace lines 
trace_lines_gdf = gpd.GeoDataFrame(geometry=trace_lines, crs=points_gdf.crs) #NF: you generate tracelines for each event, but they don't change, as the points in the levee are the same and the waterway is the same. Just make one traceline output.
trace_lines_gdf.to_file(f'J:\\IE\\Projects\\03_Southern\\IA5000TK\\07 Technical\\04 Scripting_Code\\Python\\Review\\Outputs\\{label}_trace_lines_central_sample_ReviewTest.shp')

########## - STEP 6 - LOS PROCESSING - ##########

# Specifying the diff columns to search for the critical LoS
# The [] command specifies an empty python list
# A python list is one-dimensional, i.e. just a row of values
check_cols = []
numeric_cols = ["Elev(mAHD)"]

# For loop to group numeric datasets
for label in flood_extent.keys():

    # Check cols includes _Diff values for subsequent LoS testing
    check_cols.append(f'{label}_D')
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
    
########## - STEP 7 - SAVE THE OUTPUT - ##########

# Extracting asset name to use in file export
code = str(levee_shp_reorder.iloc[0]['LOCATION_I'])
code = code.split("/")[0]

# Save points to a new shapefile
points_gdf = points_gdf.round(3)
points_gdf.to_file(f"J:\\IE\\Projects\\03_Southern\\IA5000TK\\07 Technical\\04 Scripting_Code\\Python\\Review\\Outputs\\{code}_Levee_LoS_Processing_Central_Sample_ReviewTest.shp", driver='ESRI Shapefile')