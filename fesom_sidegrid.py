from numpy import *
from netCDF4 import Dataset
from fesom_grid import *

# Classes and routines to extract a zonal slice (depth vs latitude) of the
# FESOM grid


# SideNode object containing longitude, latitude, depth, value of variable
# (chosen earlier by the user), and the node directly below it (if there is one)
# These SideNodes are interpolated between the original grid Nodes and represent
# the intersections of the original grid Elements with the specified longitude.
class SideNode:

    # Initialise with location and variable data
    def __init__ (self, lon, lat, depth, var):

        self.lon = lon
        self.lat = lat
        self.depth = depth
        self.var = var
        self.below = None

    # Save the node directly below the current node
    def set_below (self, snode_below):

        self.below = snode_below


# SideNodePair object containing two SideNodes which will later be used for
# the boundaries of SideElements.
class SideNodePair:

    # Initialise with two SideNodes
    def __init__ (self, snode1, snode2):

        # Figure out which is further south
        if snode1.lat < snode2.lat:
            self.south = snode1
            self.north = snode2
        else:
            self.south = snode2
            self.north = snode1


# SideElement object containing the four SideNodes making up the trapezoidal
# element (intersection of an original grid Element, and its 3D extension down
# through the water column, with the user-defined zonal slice).
class SideElement:

    # Initialise with four nodes (assumed to trace continuously around the
    # border of the SideElement, i.e. not jump between diagonal corners)
    def __init__ (self, snode1, snode2, snode3, snode4):

        self.snodes = array([snode1, snode2, snode3, snode4])
        lat = array([snode1.lat, snode2.lat, snode3.lat, snode4.lat])
        depth = array([snode1.depth, snode2.depth, snode3.depth, snode4.depth])
        self.y = lat
        # Make depth negative
        self.z = -depth
        # Set the value of the user-defined variable to be the mean of the
        # values at each corner (not quite mathematically correct but this
        # is just a visualisation, and much easier than integrating around a
        # trapezoid!!)
        self.var = mean(array([snode1.var, snode2.var, snode3.var, snode4.var]))


# Function to read FESOM grid files and FESOM output, and build SideElement mesh
# Input:
# mesh_path = path to directory containing grid files
# file_path = path to FESOM output file
# var_name = string containing name of variable in FESOM output file
# tstep = time index in FESOM output file
# lon0 = longitude to use for zonal slice
# lat_max = northernmost latitude to plot (generally -50, depends on your
#           definition of the Southern Ocean)
# Output:
# selement = array of SideElements making up the zonal slice
def fesom_sidegrid (mesh_path, file_path, var_name, tstep, lon0, lat_max):

    # Build the regular FESOM grid
    elm2D = fesom_grid(mesh_path)

    # Read data
    id = Dataset(file_path, 'r')
    data = id.variables[var_name][tstep-1,:]
    id.close()    

    snode_pairs = []
    for elm in elm2D:
        # Don't consider elements outside the Southern Ocean
        if any(elm.y <= lat_max):
            # Select elements which intersect lon0
            if any(elm.x <= lon0) and any(elm.x >= lon0):
                # Special case where nodes (corners) of the element are exactly
                # at longitude lon0
                if any(elm.x == lon0):
                    # If exactly one of the corners is at lon0, ignore it; this
                    # element only touches lon0 at one point
                    # If two of the corners are at lon0, an entire side of the
                    # element lies along the line lon0
                    if count_nonzero(elm.x == lon0) == 2:
                        # Select these two Nodes
                        index = nonzero(elm.x == lon0)
                        nodes = elm.nodes[index]
                        node1 = nodes[0]
                        node2 = nodes[1]
                        # Convert to SideNodes and add them to snode_pairs
                        coincide_snode(node1, node2, data, snode_pairs)
                    # Impossible for all three corners to be at lon0
                else:
                    # Regular case
                    snodes_curr = []
                    # Find the two sides of the triangular element which
                    # intersect longitude lon0
                    # For each such side, interpolate a SideNode between the
                    # two endpoint Nodes.
                    if any(array([elm.x[0], elm.x[1]]) < lon0) and any(array([elm.x[0], elm.x[1]]) > lon0):
                        snodes_curr.append(interp_snode(elm.nodes[0], elm.nodes[1], lon0, data))
                    if any(array([elm.x[1], elm.x[2]]) < lon0) and any(array([elm.x[1], elm.x[2]]) > lon0):
                        snodes_curr.append(interp_snode(elm.nodes[1], elm.nodes[2], lon0, data))
                    if any(array([elm.x[0], elm.x[2]]) < lon0) and any(array([elm.x[0], elm.x[2]]) > lon0):
                        snodes_curr.append(interp_snode(elm.nodes[0], elm.nodes[2], lon0, data))
                    # Add the two resulting SideNodes to snode_pairs
                    snode_pairs.append(SideNodePair(snodes_curr[0], snodes_curr[1]))

    selements = []
    # Build the trapezoidal SideElements
    for pair in snode_pairs:
        # Start at the surface
        snode1_top = pair.south
        snode2_top = pair.north
        while True:
            # Select the SideNodes directly below
            snode1_bottom = snode1_top.below
            snode2_bottom = snode2_top.below
            if snode1_bottom is None or snode2_bottom is None:
                # Reached the bottom, so stop
                break
            # Make a SideElement from these four SideNodes
            # The order they are passed to the SideElement initialisation
            # function is important: must trace continuously around the
            # border of the SideElement, i.e. not jump between diagonal corners
            selements.append(SideElement(snode1_top, snode2_top, snode2_bottom, snode1_bottom))
            # Get ready for the next SideElement below
            snode1_top = snode1_bottom
            snode2_top = snode2_bottom

    return selements


# Process the special case where an entire side of a triangular Element lies
# on the line longitude=lon0. Convert the endpoint Nodes to SideNodes, and add
# to the snode_pairs list.
# Input:
# node1, node2 = endpoint Nodes from this Element
# lon0 = longitude for zonal slice
# data = FESOM output on original grid
# snode_pairs = list of SideNodePair objects to add to
def coincide_snode (node1, node2, data, snode_pairs):

    # Convert the Nodes into SideNodes
    snode1 = SideNode(node1.lon, node1.lat, node1.depth, data[node1.id])
    snode2 = SideNode(node2.lon, node2.lat, node2.depth, data[node2.id])
    # Save to SideNodePair list
    snode_pairs.append(SideNodePair(snode1, snode2))

    # Travel down the water column to similarly process the Nodes at each
    # depth level
    while True:
        # Find the Nodes directly below
        node1 = node1.below
        node2 = node2.below
        if node1 is None or node2 is None:
            # We've reached the bottom; stop
            break
        # Convert these to SideNodes
        snode1_below = SideNode(node1.lon, node1.lat, node1.depth, data[node1.id])
        snode2_below = SideNode(node2.lon, node2.lat, node2.depth, data[node2.id])
        # Save to linked list
        snode1.set_below(snode1_below)
        snode2.set_below(snode2_below)
        # Get ready for next iteration
        snode1 = snode1_below
        snode2 = snode2_below


# Given two Nodes where the straight line (in lon-lat space) between them
# intersects the line longitude=lon0, calculate the latitude of this
# intersection and linearly interpolate the model output at this intersection.
# Input:
# node1, node2 = Nodes at the endpoints of this line
# lon0 = longitude to interpolate to
# data = FESOM output on original grid
# Output:
# snode_sfc = SideNode object for the intersection at the surface, with all
#             SideNodes beneath it also interpolated and linked in
def interp_snode (node1, node2, lon0, data):

    # Calculate latitude at the intersection using basic equation of a line
    lat0 = node1.lat + (node2.lat - node1.lat)/(node2.lon - node1.lon)*(lon0 - node1.lon)
    # Calculate distances from intersection to node1 (d1) and to node2 (d2)
    d1 = sqrt((lon0 - node1.lon)**2 + (lat0 - node1.lat)**2)
    d2 = sqrt((lon0 - node2.lon)**2 + (lat0 - node2.lat)**2)
    # Save the values of the given variable at each Node
    var1 = data[node1.id]
    var2 = data[node2.id]
    # Linearly interpolate the variable at the intersection
    var0 = var1 + (var2 - var1)/(d2 + d1)*d1
    depth = 0.0
    # Create a surface SideNode at this intersection
    snode_sfc = SideNode(lon0, lat0, depth, var0)

    # Now travel down the water column to interpolate the intersection at
    # each depth level
    snode = snode_sfc
    while True:
        # Find the Nodes directly below
        node1 = node1.below
        node2 = node2.below
        if node1 is None or node2 is None:
            # We've reached the bottom; stop
            break
        # Latitude and distances will not change; just interpolate the new var0
        var1 = data[node1.id]
        var2 = data[node2.id]
        var0 = var1 + (var2 - var1)/(d2 + d1)*d1
        # Similarly interpolate depth
        depth = node1.depth + (node2.depth - node1.depth)/(d2 + d1)*d1
        # Create a SideNode at this depth
        snode_below = SideNode(lon0, lat0, depth, var0)
        # Add to the linked list
        snode.set_below(snode_below)
        # Get ready for next iteration
        snode = snode_below

    return snode_sfc

    
                

    

    
                        
                     
    
