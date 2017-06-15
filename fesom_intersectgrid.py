from numpy import *
from netCDF4 import Dataset
from fesom_grid import *
from unrotate_vector import *

# Classes and routines to zonally average the FESOM grid between given longitude
# bounds, creating a regular latitude x depth grid

r = 6.371e6
deg2rad = pi/180.0


# IntersectElement object containing 2 longitudes, latitude, side length dx in
# metres, and an array of  values for the variable (chosen earlier by the user)
# corresponding to each depth in depth_vals (defined in the main function).
# At depths which  represent non-ocean points (i.e. seafloor or ice shelf), var
# will be NaN.
# These IntersectElements are interpolated between the original grid Nodes and
# represent the intersections of the original grid Elements with the specified
# latitude.
class IntersectElement:

    # Initialise with location
    def __init__ (self, lon, lat):

        self.lon = lon  # List of size 2
        self.dx = r*cos(lat*deg2rad)*abs(lon[0]-lon[1])*deg2rad
        self.lat = lat
        self.var = []


# Function to read FESOM files and FESOM output, and zonally average between
# the given longitude bounds.
# Input:
# mesh_path = path to directory containing grid files
# file_path = path to FESOM output file
# var_name = string containing name of variable in FESOM output file
# tstep = time index in FESOM output file
# lon_min, lon_max = bounds on latitudes to process
# lat_min, lat_max = bounds on latitudes to process
# depth_min, depth_max = bounds on depths to process (depth_min is deepest),
#                        negative, in metres
# num_lat, num_depth = number of regular intervals to split
# Output:
# lat_vals = latitude values on the regular lat x depth interpolated grid
# depth_vals = depth values on the regular grid (negative, in metres)
# data_reg = zonally averaged data on the regular grid
def fesom_intersectgrid (mesh_path, file_path, var_name, tstep, lon_min, lon_max, lat_min, lat_max, depth_min, depth_max, num_lat, num_depth):

    if lon_min == -180 and lon_max == 180:
        lon_bounds = False
    else:
        lon_bounds = True

    # Build the regular FESOM grid
    elements = fesom_grid(mesh_path, cross_180=False)

    # Read data
    id = Dataset(file_path, 'r')
    data = id.variables[var_name][tstep-1,:]
    # Check for vector variables that need to be unrotated
    if var_name in ['u', 'v']:
        # Read the rotated lat and lon
        fid = open(mesh_path + 'nod3d.out', 'r')
        fid.readline()
        lon = []
        lat = []
        for line in fid:
            tmp = line.split()
            lon_tmp = float(tmp[1])
            lat_tmp = float(tmp[2])
            if lon_tmp < -180:
                lon_tmp += 360
            elif lon_tmp > 180:
                lon_tmp -= 360
            lon.append(lon_tmp)
            lat.append(lat_tmp)
        fid.close()
        lon = array(lon)
        lat = array(lat)
        if var_name == 'u':
            u_data = data[:]
            v_data = id.variables['v'][tstep-1,:]
            u_data_lonlat, v_data_lonlat = unrotate_vector(lon, lat, u_data, v_data)
            data = u_data_lonlat[:]
        elif var_name == 'v':
            v_data = data[:]
            u_data = id.variables['u'][tstep-1,:]
            u_data_lonlat, v_data_lonlat = unrotate_vector(lon, lat, u_data, v_data)
            data = v_data_lonlat[:]
    id.close()

    # Build the regular grid
    lat_vals = linspace(lat_min, lat_max, num_lat)
    # Make depth positive to match the "depth" attribute in grid Nodes
    depth_vals = -1*linspace(depth_min, depth_max, num_depth)

    # Set up array of NaNs to overwrite with zonally averaged data
    data_reg = zeros((num_depth, num_lat))
    data_reg[:,:] = NaN

    # Process one latitude value at a time
    for j in range(num_lat):
        ielm_list = []
        # Loop over 2D grid Elements
        for elm in elements:
            # Select elements which intersect the current latitude, and which
            # fall entirely between the longitude bounds
            if lon_bounds:
                keep = any(elm.y <= lat_vals[j]) and any(elm.y >= lat_vals[j]) and all(elm.x >= lon_min) and all(elm.x <= lon_max)
            else:
                # No bounds on longitude
                keep = any(elm.y <= lat_vals[j]) and any(elm.y >= lat_vals[j])
            if keep:
                # Create an IntersectElement
                ielm = create_ielm(elm, lat_vals[j], depth_vals, data)
                # Check for cases where the Element intersected the given
                # latitude at exactly one corner; these aren't useful
                if ielm is not None:                    
                    ielm_list.append(ielm)
        # Zonally average at each depth
        for k in range(num_depth):
            # Set up integrals of var*dx and dx
            int_vardx = 0
            int_dx = 0
            for ielm in ielm_list:
                # Check if data exists at the current depth level
                if ielm.var[k] is not NaN:
                    int_vardx += ielm.var[k]*ielm.dx
                    int_dx += ielm.dx
            if int_dx > 0:
                data_reg[k,j] = int_vardx/int_dx

    # Convert depth back to negative for plotting
    depth_vals = -1*depth_vals

    return lat_vals, depth_vals, data_reg



# Given an Element which intersects the line latitude=lat0, create an
# IntersectElement and linearly interpolate the model output at each given
# depth level.
# Input:
# elm = Element crossing lat0
# lat0 = latitude to interpolate to 
# depth_vals = array of depth values (positive, in metres) to interpolate to
# data = FESOM output on original grid
# Output: ielm = new IntersectElement
def create_ielm (elm, lat0, depth_vals, data):

    # Special case where Nodes (corners) of the Element are exactly at lat0
    if any(elm.y == lat0):
        # If exactly one of the corners is at lat0, don't proceed; this element
        # only touches lat0 at one point
        if count_nonzero(elm.y == lat0) == 1:
            # Return an empty IntersectElement
            ielm = None
        else:
            # Impossible for all three corners to be at lat0
            # So two of the corners are. Find them.
            index = nonzero(elm.y == lat0)
            nodes = elm.nodes[index]
            node1 = nodes[0]
            node2 = nodes[1]
            # Now an entire side of the Element lies along lat0, which
            # makes things easy.
            # Initialise the IntersectElement
            ielm = IntersectElement([node1.lon, node2.lon], lat0)
            # Interpolate data to each depth value
            for depth in depth_vals:
                id1A, id1B, coeff1A, coeff1B = node1.find_depth(depth)
                id2A, id2B, coeff2A, coeff2B = node2.find_depth(depth)
                if any(isnan(array([id1A, id1B, coeff1A, coeff1B, id2A, id2B, coeff2A, coeff2B]))):
                    # At least one of the two nodes doesn't have ocean data
                    # here (i.e. seafloor or ice shelf); save a NaN
                    ielm.var.append(NaN)
                else:
                    # Linearly interpolate for each node
                    var1 = coeffA*data[idA] + coeffB*data[idB]
                    var2 = coeffC*data[idC] + coeffD*data[idD]
                    # Save the average
                    ielm.var.append(0.5*(var1+var2))
    # Regular case where no Nodes are exactly at lat0
    else:
        # Figure out which two sides of the Element intersect lat0
        # For each side, interpolate data to the intersection at each depth
        var1 = None
        var2 = None
        if any(array([elm.y[0], elm.y[1]]) < lat0) and any(array([elm.y[0], elm.y[1]]) > lat0):
            lon1, var1 = interp_intersection(elm.nodes[0], elm.nodes[1], lat0, depth_vals, data)
        if any(array([elm.y[1], elm.y[2]]) < lat0) and any(array([elm.y[1], elm.y[2]]) > lat0):
            if var1 is None:
                lon1, var1 = interp_intersection(elm.nodes[1], elm.nodes[2], lat0, depth_vals, data)
            else:
                lon2, var2 = interp_intersection(elm.nodes[1], elm.nodes[2], lat0, depth_vals, data)
        if any(array([elm.y[0], elm.y[2]]) < lat0) and any(array([elm.y[0], elm.y[2]]) > lat0):
            lon2, var2 = interp_intersection(elm.nodes[0], elm.nodes[2], lat0, depth_vals, data)
        # Initialise IntersectElement for these intersections
        ielm = IntersectElement([lon1, lon2], lat0)
        # Save data for each depth value
        for k in range(len(depth_vals)):
            if var1[k] is NaN or var2[k] is NaN:
                # At least one of the two nodes doesn't have ocean data here
                # (i.e. seafloor or ice shelf); save a NaN
                ielm.var.append(NaN)
            else:
                # Save the average
                ielm.var.append(0.5*(var1[k] + var2[k]))

    return ielm


# Given two Nodes where the straight line (in lon-lat space) between them
# intersects the line latitude=lat0, calculate the longitude of this
# intersection and linearly interpolate the model output at this intersection,
# for each given depth.
# Input:
# node1, node2 = Nodes at the endpoints of this line
# lat0 = latitude to interpolate to
# depth_vals = array of depth values to interpolate to (positive, in metres)
# data = FESOM output on original grid
# Output:
# lon0 = longitude of intersection
# var = data interpolated to the intersection for each depth
def interp_intersection (node1, node2, lat0, depth_vals, data):

    # Calculate longitude at the intersection using basic equation of a line
    lon0 = (lat0 - node1.lat)*(node2.lon - node1.lon)/(node2.lat - node1.lat) + node1.lon
    # Calculate distances from intersection to node1 (d1) and node2 (d2)
    d1 = sqrt((lon0 - node1.lon)**2 + (lat0 - node1.lat)**2)
    d2 = sqrt((lon0 - node2.lon)**2 + (lat0 - node2.lat)**2)

    # Interpolate model output to each given depth
    var = []
    for depth in depth_vals:
        id1A, id1B, coeff1A, coeff1B = node1.find_depth(depth)
        id2A, id2B, coeff2A, coeff2B = node2.find_depth(depth)
        if any(isnan(array([id1A, id1B, coeff1A, coeff1B, id2A, id2B, coeff2A, coeff2B]))):
            # At least one of the Nodes has no ocean data here (i.e. it is
            # seafloor or ice shelf); save NaN
            var.append(NaN)
        else:
            # Linearly interpolate to each given Node at this depth
            var1 = coeff1A*data[id1A] + coeff1B*data[id1B]
            var2 = coeff2A*data[id2A] + coeff2B*data[id2B]
            # Linearly interpolate to the intersection
            var.append(var1 + (var2 - var1)/(d2 + d1)*d1)

    return lon0, var
    
