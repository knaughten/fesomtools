from numpy import *
from netCDF4 import Dataset
from fesom_grid import *
from unrotate_vector import *

# Classes and routines to zonally average the FESOM grid between given longitude
# bounds, creating a regular latitude x depth grid


# IntersectNode object containing longitude, latitude, and an array of values
# for the variable (chosen earlier by the user) corresponding to each depth in
# depth_vals (defined in the main function). At depths which represent non-ocean
# points (i.e. seafloor or ice shelf), var will be NaN.
# These IntersectNodes are interpolated between the original grid Nodes and
# represent the intersections of the original grid Elements with the specified
# latitude.
class IntersectNode:

    # Initialise with location
    def __init__ (self, lon, lat):

        self.lon = lon
        self.lat = lat
        self.var = []


# Function to read FESOM files and FESOM output, and zonally average between
# the given longitude bounds.
# Input:
# mesh_path = path to directory containing grid files
# file_path = path to FESOM output file
# var_name = string containing name of variable in FESOM output file
# tstep = time index in FESOM output file
# lon_min, lon_max = bounds on longitude for the zonal average
# lat_min, lat_max = bounds on latitudes to process (generally -90 to -50)
# depth_min, depth_max = bounds on depths to process (depth_min is deepest),
#                        negative, in metres
# num_lat, num_depth = number of regular intervals to split 
# Output:
# lat_vals = latitude values on the regular lat x depth interpolated grid
# depth_vals = depth values on the regular grid (negative, in metres)
# data_reg = zonally averaged data on the regular grid
def fesom_intersectgrid (mesh_path, file_path, var_name, tstep, lon_min, lon_max, lat_min, lat_max, depth_min, depth_max, num_lat, num_depth):

    # Build the regular FESOM grid
    elements = fesom_grid(mesh_path)

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
        inodes_lat = []
        # Loop over 2D grid Elements
        for elm in elements:
            # Select elements which intersect the current latitude, and which
            # fall entirely between the longitude bounds
            if any(elm.y <= lat_vals[j]) and any(elm.y >= lat_vals[j]) and all(elm.x >= lon_min) and all(elm.x <= lon_max):
                # Special case where nodes (corners) of the element are exactly
                # at lat_vals[j]
                if any(elm.y == lat_vals[j]):
                    # If exactly one of the corners is at lat_vals[j], ignore
                    # it; this element only touches lat_vals[j] at one point
                    # If two of the corners are at lat_vals[j], an entire side
                    # of the element lies along the line lat_vals[j]
                    if count_nonzero(elm.y == lat_vals[j]) == 2:
                        # Select these two Nodes
                        index = nonzero(elm.y == lat_vals[j])
                        nodes = elm.nodes[index]
                        node1 = nodes[0]
                        node2 = nodes[1]
                        # Convert to IntersectNodes and add them to inodes_lat
                        inodes_lat.append(coincide_inode(node1, depth_vals, data))
                        inodes_lat.append(coincide_inode(node2, depth_vals, data))
                    # Impossible for all three corners to be at lat_vals[j]
                else:
                    # Regular case
                    # Find the two sides of the triangular element which
                    # intersect lat_vals[j]
                    # For each such side, interpolate an IntersectNode between
                    # the two endpoint nodes, and add them to inodes_lat
                    if any(array([elm.y[0], elm.y[1]]) < lat_vals[j]) and any(array([elm.y[0], elm.y[1]]) > lat_vals[j]):
                        inodes_lat.append(interp_inode(elm.nodes[0], elm.nodes[1], lat_vals[j], depth_vals, data))
                    if any(array([elm.y[1], elm.y[2]]) < lat_vals[j]) and any(array([elm.y[1], elm.y[2]]) > lat_vals[j]):
                        inodes_lat.append(interp_inode(elm.nodes[1], elm.nodes[2], lat_vals[j], depth_vals, data))
                    if any(array([elm.y[0], elm.y[2]]) < lat_vals[j]) and any(array([elm.y[0], elm.y[2]]) > lat_vals[j]):
                        inodes_lat.append(interp_inode(elm.nodes[0], elm.nodes[2], lat_vals[j], depth_vals, data))

        # Sort inodes_lat by longitude (ascending)
        inodes_lat.sort(key=lambda inode: inode.lon)

        # Interpolate the variable values at each depth
        for k in range(num_depth):
            valid_lon = []
            valid_var = []
            for inode in inodes_lat:
                # Select all IntersectNodes where data exists at the current
                # depth level
                if inode.var[k] is not nan:
                    # Only continue if an identical inode (same longitude)
                    # hasn't already been added to valid_lon and valid_var
                    # (this will happen on adjacent elements which share a side)
                    if inode.lon not in valid_lon:
                        # Save longitude and variable values
                        valid_lon.append(inode.lon)
                        valid_var.append(inode.var[k])
            # Convert to numpy arrays so we can do math with them
            valid_lon = array(valid_lon)
            valid_var = array(valid_var)
            if len(valid_lon) == 0:
                # No valid data; leave data_reg[k,j] as NaN
                pass
            elif len(valid_lon) == 1:
                # Only one valid data point; save to data_reg
                data_reg[k,j] = valid_var[0]
            else:
                # Average over longitude
                # Trapezoidal rule for integration
                dlon = valid_lon[1:] - valid_lon[0:-1]
                var_centres = 0.5*(valid_var[0:-1] + valid_var[1:])
                # Divide integral of var_centres*dlon by integral of dlon
                # to get average; save to data_reg
                var_avg = sum(var_centres*dlon)/sum(dlon)
                data_reg[k,j] = var_avg

    # Convert depth back to negative for plotting
    depth_vals = -1*depth_vals

    return lat_vals, depth_vals, data_reg


# Process the special case where an entire side of a triangular Element lies
# on the line of constant latitude. Convert the given endpoint to an
# IntersectNode.
# Input:
# node = endpoint Node from this Element
# depth_vals = array of depth values to interpolate to (positive, in metres)
# data = FESOM output on original grid
# Output: inode = new IntersectNode
def coincide_inode (node, depth_vals, data):

    # Initialise an IntersectNode
    inode = IntersectNode(node.lon, node.lat)

    # Interpolate data at this Node to each depth value
    for depth in depth_vals:
        idA, idB, coeffA, coeffB = node.find_depth(depth)
        if any(isnan(array([idA, idB, coeffA, coeffB]))):
            # No ocean data here (i.e. seafloor or ice shelf);
            # save a NaN
            inode.var.append(NaN)
        else:
            # Linearly interpolate and save to inode.var
            var0 = coeffA*data[idA] + coeffB*data[idB]
            inode.var.append(var0)
    return inode


# Given two Nodes where the straight line (in lon-lat space) between them
# intersects the line latitude=lat0, calculate the longitude of this
# intersection and linearly interpolate the model output at this intersection.
# Input:
# node1, node2 = Nodes at the endpoints of this line
# lat0 = latitude to interpolate to
# depth_vals = array of depth values to interpolate to (positive, in metres)
# data = FESOM output on original grid
# Output: inode = new IntersectNode
def interp_inode (node1, node2, lat0, depth_vals, data):

    # Calculate longitude at the intersection using basic equation of al ine
    lon0 = (lat0 - node1.lat)*(node2.lon - node1.lon)/(node2.lat - node1.lat) + node1.lon
    # Initialise an IntersectNode
    inode = IntersectNode(lon0, lat0)
    # Calculate distances from intersection to node1 (d1) and node2 (d2)
    d1 = sqrt((lon0 - node1.lon)**2 + (lat0 - node1.lat)**2)
    d2 = sqrt((lon0 - node2.lon)**2 + (lat0 - node2.lat)**2)

    # Interpolate model output at inode to each given depth
    for depth in depth_vals:
        # Interpolate model output at each original Node
        id1A, id1B, coeff1A, coeff1B = node1.find_depth(depth)
        id2A, id2B, coeff2A, coeff2B = node2.find_depth(depth)
        if any(isnan(array([id1A, id1B, coeff1A, coeff1B, id2A, id2B, coeff2A, coeff2B]))):
            # At least one of the Nodes has no ocean data here (i.e. it is
            # seafloor or ice shelf); save NaN to inode's var array
            inode.var.append(NaN)
        else:
            # Linearly interpolate to each given Node
            var1 = coeff1A*data[id1A] + coeff1B*data[id1B]
            var2 = coeff2A*data[id2A] + coeff2B*data[id2B]
            # Linearly interpolate to inode and save to its var array
            var0 = var1 + (var2 - var1)/(d2 + d1)*d1
            inode.var.append(var0)

    return inode

    


                                             
                 
                   
                   
        
