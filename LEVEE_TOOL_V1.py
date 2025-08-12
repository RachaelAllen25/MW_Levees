# Importing necessary python packages
# These will need to be installed on the users PC to run the script
# Refer to requirements file for package details
import os
import pandas as pd
import rasterio
import numpy as np
import geopandas as gpd
from shapely.geometry import Point, LineString

########## - STEP 1 - DEFINE FUNCTIONS - ##########

def check_file(path, expected_ext=None):
    """
    Validates the existence and file extension of a given file path

    Parameters:
    ----------
    path : str
        The full file path to check

    expected_ext : str, optional
        The expected file extension (e.g., '.tif', '.shp'). If provided, the function will verify
        that the file has this extension

    Raises:
    -------
    FileNotFoundError
        If the file does not exist at the specified path

    ValueError
        If the file exists but does not match the expected extension

    Notes:
    ------
    - This function is useful for pre-validating inputs before attempting file operations
    - The extension check is case-insensitive
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")
    if expected_ext and not path.lower().endswith(expected_ext):
        raise ValueError(f"Unexpected file format for {path}. Expected {expected_ext}")

def get_crs(filepath):
    """
    Determines the coordinate reference system (CRS) of a spatial file
    Supports vector files (.shp) and raster files (.tif, .asc, .flt)

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
      to allow proper CRS detection
    - Result files are assumed to be in a .tif, .asc or .flt format
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
    Samples elevation values from a LiDAR raster and assigns them to point geometries in a GeoDataFrame

    Parameters:
    ----------
    points_gdf : geopandas.GeoDataFrame
        A GeoDataFrame containing point geometries. This object will be modified and replaced by a new 'sample' gdf
        'Elev(mAHD)' with elevation values sampled from the raster

    Lidar_filepath : str
        File path to the LiDAR raster (e.g., a .tif file) from which elevation values will be extracted

    Notes:
    ------
    - Elevation values are sampled using the raster index corresponding to each point's coordinates
    - The function assumes the raster and point geometries are in the same coordinate reference system (CRS)

    Returns:
    --------
    sampled_points_gdf : geopandas.GeoDataFrame
        A GeoDataFrame with sampled LiDAR values appended
    """
    # Creates copy of input data
    sampled_points_gdf = points_gdf.copy()
    sampled_points_gdf["Elev(mAHD)"] = None # Initate column

    # Open raster
    with rasterio.open(Lidar_filepath) as src:
        
        # Loop through each point and extract the raster value
        for index, point in sampled_points_gdf.iterrows():
            x, y = point.geometry.x, point.geometry.y
            row, col = src.index(x, y)
            sampled_points_gdf.at[index, 'Elev(mAHD)'] = src.read(1)[row, col]

    return sampled_points_gdf            

def sample_flood_extent(points_gdf, raster_path, label):
    """
    Samples flood extent values from a raster and assigns them to point geometries in a GeoDataFrame

    Parameters:
    ----------
    points_gdf : geopandas.GeoDataFrame
        A GeoDataFrame containing point geometries. This object will be modified in place by adding a new column
        named '{label}_L' with values sampled from the raster

    raster_path : str
        File path to the flood extent raster (e.g., a .tif file) from which values will be extracted.

    label : str
        A label used to name the new column in the GeoDataFrame. For example, if label is 'AEP_1%', the new column
        will be 'AEP_1%_L'

    Notes:
    ------
    - The function assumes the raster and point geometries are in the same coordinate reference system (CRS)
    - This function is not called in the current version of the LEVEE_TOOL, however has been left in for future iterations/changes

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
    Extracts the suffix and numeric identifier from a levee code string

    Parameters:
    ----------
    levee_code : str
        A string representing a levee identifier, expected to follow the format 'LEV<number><suffix>'
        Example: 'LEV12A' → number: 12, suffix: 'A'

    Returns:
    --------
    tuple
        A tuple containing:
        - suffix (str): The last character of the levee code
        - number (int): The numeric portion of the levee code, extracted from between 'LEV' and the suffix

    Notes:
    ------
    - The function assumes the levee code always starts with 'LEV' and ends with a single-character suffix
    - If the format varies, consider adding error handling to validate the input
    - This function was developed to help ease post processing ordering of levee assets
    """
    suffix = levee_code[-1]
    number = int(levee_code.split("LEV")[1][:-1])
    return suffix, number

def trace_and_sample(point, waterway_lines, raster, nodata_value=-9999):
    """
    Traces a point to the nearest waterway line and samples a raster at the nearest location

    Parameters:
        point (shapely.geometry.Point): The reference point
        waterway_lines (GeoSeries): Geometry of waterway lines
        raster (rasterio.DatasetReader): Open raster object
        nodata_value (float): Value to return if raster sample is missing or invalid

    Returns:
        tuple: (trace_line: LineString, sampled_value: float)

    Notes:
    ------
    - -9999 is an assigned nodata value
    - Users can update this based on their own preference
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

    # Function outputs
    return trace_line, sampled_value

########## - STEP 2 - DEFINE FILEPATHS & CALL DATA CHECK FUNCTIONS - ##########

# Define file paths
Levee_filepath = r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\0503\03 Exports\0503_Eumemmering_Ck_LeveeAlign_v3.shp"
Lidar_filepath = r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\0503\03 Exports\LiDAR_Merge.tif"
waterway = r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\0503\02 Build Data\03 Shapefiles\0503_Channel_Centreline.shp"

# Extract lidar extension for input into check_file function
# Levee and waterway files assumed to be in .shp format
_, levee_ext = os.path.splitext(Levee_filepath)
_, lidar_ext = os.path.splitext(Lidar_filepath)
_, waterway_ext = os.path.splitext(waterway)

# Define flood result file paths
# Flooding results loaded in with label function - this creates a dictionary of titles and filepaths
# Comment out any extra filepaths/used AEP events
flood_extent = {
    '1%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\0503\03 Exports\Rev1\0503_EumemmeringCk_3m_aep1_h_Max.tif",
    # '1%CC_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\0621\01 Build Data\TUFLOW_results_H_max\GDA2020_Reprojected\Rasters\0621_100yr_CC_GDA2020.tif",
    '2%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\0503\03 Exports\Rev1\0503_EumemmeringCk_3m_aep2_h_Max.tif",
    '5%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\0503\03 Exports\Rev1\0503_EumemmeringCk_3m_aep5_h_Max.tif",
    '10%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\0503\03 Exports\Rev1\0503_EumemmeringCk_3m_aep10_h_Max.tif",
    # '10%CC_AEP' : r'J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_CC_C_010y_657_max_h.asc',
    '20%_AEP' : r"J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\0503\03 Exports\Rev1\0503_EumemmeringCk_3m_aep20_h_Max.tif",
    #'20%CC_AEP' : r'J:\IE\Projects\03_Southern\IA5000TK\07 Technical\02 Hydraulics\4310\01 Build Data\TUFLOW_results_H_max\Arden_fail_CC_C_005y_657_max_h.asc'
 }

# Validate input files
check_file(Levee_filepath, levee_ext)
check_file(Lidar_filepath, lidar_ext)
check_file(waterway, waterway_ext)

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

# Loop to manage flood dictionary
for label, path in flood_extent.items():
    _, flood_ext = os.path.splitext(path)
    check_file(path, flood_ext)
    flood_crs = get_crs(path)
    print(f"The {label} result is projected in {flood_crs}")
    if flood_crs != levee_crs:
        print(f"⚠️ Warning: CRS mismatch between levee shapefile and flood extent raster for {label}")

########## - STEP 3 - READ IN LEVEE SHAPEFILE & REFORMAT - ##########

# Read in levee shapefile 
levee_shp = gpd.read_file(Levee_filepath) 

# Point on the levee - create an empty list 
points_list = []
chainage_dataframe = gpd.GeoDataFrame({
    'Levee Identifier': pd.Series(dtype='object'),
    'Chainage': pd.Series(dtype = 'float64')})

# This is in metres
distance = 10 

# Create unique list of levee names based on MW data field
code_levee_names = levee_shp['LOCATION_I'].unique() # Code assume LOCATION_I is a field within the input levee shapefile as per MW geodatabase for levees

# Empty dictionary to store reordered data
levee_order = {}

# For loop to cycle through each unique name
for asset_code in code_levee_names:
        
    # Earlier function is called to generate the levee name parts
    suffix, number = levee_identifier(asset_code)

    # This line is grouping the levee names into their respective parts
    levee_order.setdefault(suffix, []).append((number, asset_code))

# This is the equivalent to a nested for loop and is to sort the levees in an order convenient for post processing
reordered_levee_codes = [code for suffix in sorted(levee_order) for _, code in sorted(levee_order[suffix])]
print(f"Levee assets within input shapefile include {reordered_levee_codes}")

# Reindexing the levee shp based on new index order
levee_shp_reorder = levee_shp.copy()
levee_shp_reorder = levee_shp_reorder.set_index('LOCATION_I').loc[reordered_levee_codes]
levee_shp_reorder.reset_index(inplace = True)

# Initialising datasets for chainage reordering
looped_chain_count = 0
levee_count = 0

# For loop to work through each line feature in the input shapefile, i.e. each specific levee asset
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
        chainage_dataframe.at[(looped_chain_count + point_loc), 'Levee Identifier'] = reordered_levee_codes[(levee_count)]

        # If the point_loc is the last point, we can add the end point of the levee line using
        if point_loc == points:
             
             # Append the last point of the levee line
             end_point = line.interpolate(length)
             points_list.append(Point(end_point.x, end_point.y))
            
             # Update chainage for the end point
             unique_chainage = length
             chainage_dataframe.at[(looped_chain_count + point_loc + 1), 'Chainage'] = unique_chainage
             chainage_dataframe.at[(looped_chain_count + point_loc + 1), 'Levee Identifier'] = reordered_levee_codes[(levee_count)]
        
    # Resetting loop counter for next individual levee asset
    looped_chain_count = looped_chain_count + point_loc + 2 # +2 to move past end point
    levee_count = levee_count + 1

# Create a new GeoDataFrame from the points
# Crs same as input shapefile 
points_gdf = gpd.GeoDataFrame(geometry=points_list, crs=levee_shp_reorder.crs)

# Connecting points and chainage dataframe
points_gdf = pd.concat([points_gdf, chainage_dataframe], axis = 1)
print("Point delineation complete") 

########## - STEP 4 - CALL SAMPLE FUNCTIONS - ##########

# Sample lidar 
sampled_points_gdf = sample_lidar_values(points_gdf, Lidar_filepath)
print("LiDAR point sampling complete")

# Central waterway sampling
# create a line (vector) from the levee point to the nearest waterway 
# Sample flood result at location of closest point on the waterway
# Loop through raster list of flood extents 
waterway_gdf = gpd.read_file(waterway)
if sampled_points_gdf.crs != waterway_gdf.crs:
    waterway_gdf = waterway_gdf.to_crs(sampled_points_gdf.crs) 

for label, raster_path in flood_extent.items():
    with rasterio.open(raster_path) as raster:
       
        # Create empty lists within the loop, so they are generated for each raster
        trace_lines = []
        central_waterway_values = []

        # check coordinate system
        if sampled_points_gdf.crs != waterway_gdf.crs:
            waterway_gdf = waterway_gdf.to_crs(sampled_points_gdf.crs) 
        
        # This command strips the attribute data from the geodataframe and leaves only the geometries
        waterway_lines = waterway_gdf.geometry

        # This for loop interates over both the index and the geometry (points) with the geodata frame
        for idx, point in sampled_points_gdf.geometry.items():

            # Waterway/result sample function called
            trace_line, value = trace_and_sample(point, waterway_gdf.geometry, raster)
            trace_lines.append(trace_line)
            central_waterway_values.append(value)
            
        # Add to original GeoDataFrame 
        sampled_points_gdf[f'{label}_W'] = central_waterway_values

        # Computing difference/freeboard analysis
        sampled_points_gdf[f'{label}_D'] = sampled_points_gdf['Elev(mAHD)'] - sampled_points_gdf[f'{label}_W']
        
# Progress print statement
print('Central waterway sampling complete')

########## - STEP 6 - LOS PROCESSING - ##########

# Specifying the diff columns to search for the critical LoS
check_cols = []
numeric_cols = ["Elev(mAHD)"]

# For loop to group numeric datasets
for label in flood_extent.keys():

    # Check cols includes _Diff values for subsequent LoS testing
    check_cols.append(f'{label}_D')
    numeric_cols.append(f'{label}_W')
    numeric_cols.append(f'{label}_D')

# Creating a subset of the points_gdf from the check_cols
# All other values not > 0 become NaN
freeboard_events = sampled_points_gdf[check_cols].where(sampled_points_gdf[check_cols] > 0)

# Taking the LoS classification corresponding to the minimum freeboard event
sampled_points_gdf['LoS'] = freeboard_events.idxmin(axis=1)

# Code to assign the levee points with no observed freeboard a LoS rank of >20%
sampled_points_gdf.loc[freeboard_events.min(axis=1).isna(), 'LoS'] = ">20%"

# Dropping the _Diff suffix from LoS entry
sampled_points_gdf['LoS'] = sampled_points_gdf['LoS'].str.split('_').str[0]

# Converting data columns to numeric to allow for rounding
sampled_points_gdf[numeric_cols] = sampled_points_gdf[numeric_cols].apply(lambda col: col.astype('float64') if col.dtype != 'float64' else col)
sampled_points_gdf[numeric_cols] = sampled_points_gdf[numeric_cols].apply(pd.to_numeric, errors='coerce').round(3)
    
########## - STEP 7 - SAVE THE OUTPUT - ##########

# Extracting asset name to use in file export
code = str(levee_shp_reorder.iloc[0]['LOCATION_I'])
code = code.split("/")[0]

# Output generated trace lines exported for anaylsis/checking of correct raster cells
trace_lines_gdf = gpd.GeoDataFrame(geometry=trace_lines, crs=sampled_points_gdf.crs)
trace_lines_gdf.to_file(f'J:\\IE\\Projects\\03_Southern\\IA5000TK\\07 Technical\\02 Hydraulics\\{code}\\04 Code Processing\\Post_Review\\{code}_Sample_Trace_Lines.shp')

# Save points to a new shapefile
sampled_points_gdf = sampled_points_gdf.round(3)
sampled_points_gdf.to_file(f"J:\\IE\\Projects\\03_Southern\\IA5000TK\\07 Technical\\02 Hydraulics\\{code}\\04 Code Processing\\Post_Review\\{code}_Levee_LoS_Matrix.gpkg", driver='GPKG')
print('LoS processing complete and exports generated')