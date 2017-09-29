from numpy import *
from netCDF4 import Dataset
from triangle_area_latdepth import *
from on_line import *
from side_of_line import *
from intersects_triangle import *

# Classes and routines to extract a depth slice along any horizontal line on
# the FESOM grid.

# Global variables
# Radius of the Earth in metres
r = 6.371e6
# Degrees to radians conversion factor
deg2rad = pi/180.0


# SideNode object containing longitude, latitude, depth, value of variable
# (chosen earlier by the user), and the node directly below it (if there is one)
# These SideNodes are interpolated between the original grid Nodes and represent
# the intersections of the original grid Elements with the specified horizontal
# line.
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


# SideElement object containing the four SideNodes making up the quadrilateral
# element (intersection of an original grid Element, and its 3D extension down
# through the water column, with the user-defined horizontal line).
class SideElement:

    # Initialise with four nodes (assumed to trace continuously around the
    # border of the SideElement, i.e. not jump between diagonal corners)
    # as well as the coordinates of the beginning of the line segment
    def __init__ (self, snode1, snode2, snode3, snode4, lon0, lat0):

        self.snodes = array([snode1, snode2, snode3, snode4])
        # Calculate great-circle distance from (lon0, lat0) in metres
        distance = []
        for snode in [snode1, snode2, snode3, snode4]:
            distance.append(r*arccos(sin(snode.lat)*sin(lat0) + cos(snode.lat)*cos(lat0)*cos(abs(snode.lon-lon))))
        self.y = array(distance)
        depth = array([snode1.depth, snode2.depth, snode3.depth, snode4.depth])
        # Make depth negative
        self.z = -depth
        # Set the value of the user-defined variable to be the mean of the
        # values at each corner (not quite mathematically correct but this
        # is just a visualisation, and much easier than integrating around a
        # quadrilateral!!)
        self.var = (snode1.var + snode2.var + snode3.var + snode4.var)/4.0

    # Return the area of the quadrilateral making up this Element
    def area (self):

        # Divide into two triangles and add the areas
        lat1 = array([self.y[0], self.y[1], self.y[2]])
        depth1 = array([self.z[0], self.z[1], self.z[2]])
        area1 = triangle_area_latdepth(lat1, depth1)
        
        lat2 = array([self.y[0], self.y[2], self.y[3]])
        depth2 = array([self.z[0], self.z[2], self.z[3]])
        area2 = triangle_area_latdepth(lat2, depth2)

        return area1 + area2    


# Function to build SideElement mesh
# Input:
# elm2D = elements from regular FESOM grid
# data = FESOM output at each node; can be a single time index or a timeseries
# lon0, lat0 = coordinates of beginning of line segment
# lon1, lat1 = coordinates of end of line segment
def fesom_sidegrid_general (elm2D, data, lon0, lat0, lon1, lat1):

    snode_pairs = []
    for elm in elm2D:
        # Only consider elements within the rectangle enclosing the line segment
        if any(elm.lon >= min(lon0, lon1)) and any(elm.lon <= max(lon0, lon1)) and any(elm.lat >= min(lat0, lat)) and any(elm.lat <= max(lat0, lat1)):
            # Select elements which intersect, overlap, or touch the line segment
            if intersects_triangle(elm, lon0, lat0, lon1, lat1):
                # Check if any of the nodes (corners) of the element fall on
                # the line
                nodes_touching = 0
                for n in range(3):
                    if on_line(elm.lon[n], elm.lat[n], lon0, lat0, lon1, lat1):
                        nodes_touching += 1
                if nodes_touching == 1:
                    # If exactly one of the nodes is on the line, check if the
                    # opposite edge of the triangle intersects the line segment
                    if on_line(elm.lon[0], elm.lat[0], lon0, lat0, lon1, lat1) and side_of_line(elm.lon[1], elm.lat[1], elm.lon[2], elm.lat[2], lon0, lat0, lon1, lat1) <= 0:
                        # Get one SideNode directly from the node on the line,
                        # and interpolate another SideNode on the edge between
                        # the other 2 nodes
                        snode_pairs.append(SideNodePair(single_coincide_snode(elm.nodes[0], data), interp_snode(elm.nodes[1], elm.nodes[2], lon0, lat0, lon1, lat1, data)))
                    elif on_line(elm.lon[1], elm.lat[1], lon0, lat0, lon1, lat1) and side_of_line(elm.lon[0], elm.lat[0], elm.lon[2], elm.lat[2], lon0, lat0, lon1, lat1) <= 0:
                        snode_pairs.append(SideNodePair(single_coincide_snode(elm.nodes[1], data), interp_snode(elm.nodes[0], elm.nodes[2], lon0, lat0, lon1, lat1, data)))
                    elif on_line(elm.lon[2], elm.lat[2], lon0, lat0, lon1, lat1) and side_of_line(elm.lon[0], elm.lat[0], elm.lon[1], elm.lat[1], lon0, lat0, lon1, lat1) <= 0:
                        snode_pairs.append(SideNodePair(single_coincide_snode(elm.nodes[2], data), interp_snode(elm.nodes[0], elm.nodes[1], lon0, lat0, lon1, lat1, data)))
                    else:
                        # The element skims across the line at exactly one
                        # point. Ignore it, because that point will be dealt
                        # with inside another element.
                        pass
                if nodes_touching == 2:
                    # If two of the nodes are on the line, an entire side of
                    # the element lies along the line.
                    # Convert these two nodes to SideNodes and add them to
                    # snode_pairs.
                    if on_line(elm.lon[0], elm.lat[0], lon0, lat0, lon1, lat1) and on_line(elm.lon[1], elm.lat[1], lon0, lat0, lon1, lat1):
                        double_coincide_snode(elm.nodes[0], elm.nodes[1], data, snode_pairs)
                    elif on_line(elm.lon[0], elm.lat[0], lon0, lat0, lon1, lat1) and on_line(elm.lon[2], elm.lat[2], lon0, lat0, lon1, lat1):
                        double_coincide_snode(elm.nodes[0], elm.nodes[2], data, snode_pairs)
                    elif on_line(elm.lon[1], elm.lat[1], lon0, lat0, lon1, lat1) and on_line(elm.lon[2], elm.lat[2], lon0, lat0, lon1, lat1):
                        double_coincide_snode(elm.nodes[1], elm.nodes[2], data, snode_pairs)
                # Impossible for all three corners to be on the line
                else:
                    # Regular case
                    snodes_curr = []
                    # Find the two edges of the triangular element which
                    # intersect (or touch) the line segment
                    # For each such edge, interpolate a SideNode between the
                    # two endpoint Nodes.
                    if side_of_line(elm.lon[0], elm.lat[0], elm.lon[1], elm.lat[1], lon0, lat0, lon1, lat1) <= 0:
                        snodes_curr.append(interp_snode(elm.nodes[0], elm.nodes[1], lon0, lat0, lon1, lat1, data))
                    elif side_of_line(elm.lon[0], elm.lat[0], elm.lon[2], elm.lat[2], lon0, lat0, lon1, lat1) <= 0:
                        snodes_curr.append(interp_snode(elm.nodes[0], elm.nodes[2], lon0, lat0, lon1, lat1, data))
                    elif side_of_line(elm.lon[1], elm.lat[1], elm.lon[2], elm.lat[2], lon0, lat0, lon1, lat1) <= 0:
                        snodes_curr.append(interp_snode(elm.nodes[1], elm.nodes[2], lon0, lat0, lon1, lat1, data))
                    # Add the two resulting SideNodes to snode_pairs
                    snode_pairs.append(SideNodePair(snodes_curr[0], snodes_curr[1]))

    selements = []
    # Build the quadrilateral SideElements
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
            selements.append(SideElement(snode1_top, snode2_top, snode2_bottom, snode1_bottom, lon0, lat0))
            # Get ready for the next SideElement below
            snode1_top = snode1_bottom
            snode2_top = snode2_bottom

    return selements


# Process the special case where a single node lies on the line (and the
# opposite edge of the triangle intersects the line, so we care).
# Convert this node to a SideNode.
# Input:
# node = Node on the line
# data = FESOM output at each node; can be a single time index or a timeseries
# Output:
# snode_sfc = SideNode object for the intersection at the surface, with all
#             SideNodes beneath it also interpolated and linked in
def single_coincide_snode(node, data):

    if len(data.shape) == 2:
        # Timeseries
        snode_sfc = SideNode(node.lon, node.lat, node.depth, data[:,node.id])
    else:
        # Single time index
        snode_sfc = SideNode(node.lon, node.lat, node.depth, data[node.id])
    # Travel down the water column to similarly process the node at each depth level
    snode = snode_sfc
    while True:
        # Find the node directly below
        node = node.below
        if node is None:
            # We've reached the bottom; stop
            break
        # Convert to SideNode
        if len(data.shape) == 2:
            # Timeseries
            snode_below = SideNode(node.lon, node.lat, node.depth, data[:,node.id])
        else:
            # Single time index
            snode_below = SideNode(node.lon, node.lat, node.depth, data[node.id])
        # Save to linked list
        snode.set_below(snode_below)
        # Get ready for next iteration
        snode = snode_below


# Process the special case where an entire side of a triangular Element lies
# on the line. Convert the endpoint Nodes to SideNodes, and add to the
# snode_pairs list.
# Input:
# node1, node2 = endpoint Nodes from this Element
# data = FESOM output at each node; can be a single time index or a timeseries
# snode_pairs = list of SideNodePair objects to add to
def double_coincide_snode(node1, node2, data, snode_pairs):

    # Convert the Nodes into SideNodes
    if len(data.shape) == 2:
        # Timeseries
        snode1 = SideNode(node1.lon, node1.lat, node1.depth, data[:,node1.id])
        snode2 = SideNode(node2.lon, node2.lat, node2.depth, data[:,node2.id])
    else:
        # Single time index
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
        if len(data.shape) == 2:
            # Timeseries
            snode1_below = SideNode(node1.lon, node1.lat, node1.depth, data[:,node1.id])
            snode2_below = SideNode(node2.lon, node2.lat, node2.depth, data[:,node2.id])
        else:
            # Single time index
            snode1_below = SideNode(node1.lon, node1.lat, node1.depth, data[node1.id])
            snode2_below = SideNode(node2.lon, node2.lat, node2.depth, data[node2.id])
        # Save to linked list
        snode1.set_below(snode1_below)
        snode2.set_below(snode2_below)
        # Get ready for next iteration
        snode1 = snode1_below
        snode2 = snode2_below


# Given two Nodes where the straight line (in lon-lat space) between them
# intersects the line we are interpolating to, calculate the location of the
# intersection and linearly interpolate the model output to that point.
# Input:
# node1, node2 = Nodes at the endpoints of the given edge of the triangle
# lon0, lat0 = coordinates of beginning of line segment
# lon1, lat1 = coordinates of end of line segment
# data = FESOM output on regular grid; can be a single time index or a 
#        timeseries
# Output:
# snode_sfc = SideNode object for the intersection at the surface, with all
#             SideNodes beneath it also interpolated and linked in
def interp_snode(node1, node2, lon0, lat0, lon1, lat1, data):

    # Calculate slope of the line containing (lon0, lat0) and (lon1, lat1)
    slope_line = (lat1 - lat0)/(lon1 - lon0)
    # Calculate slope of the triangle edge we care about
    slope_edge = (node2.lat - node1.lat)/(node2.lon - node1.lon)
    # Calculate the location of the intersection
    lon_int = ((slope_line*lon0 - lat0) - (slope_edge*node1.lon - node1.lat))/(slope_line - slope_edge)
    lat_int = slope_line*(lon_int - lon0) + lat0
    # Calculate distances from intersection to node1 (d1) and to node2 (d2)
    d1 = sqrt((lon_int - node1.lon)**2 + (lat_int - node1.lat)**2)
    d2 = sqrt((lon_int - node2.lon)**2 + (lat_int - node2.lat)**2)
    # Save the values of the given variable at each Node
    if len(data.shape) == 2:
        # Timeseries
        var1 = data[:,node1.id]
        var2 = data[:,node2.id]
    else:
        # Single time index
        var1 = data[node1.id]
        var2 = data[node2.id]
    # Linearly interpolate the variable at the intersection
    var0 = var1 + (var2 - var1)/(d2 + d1)*d1
    # Also interpolate depth
    depth = node1.depth + (node2.depth - node1.depth)/(d2 + d1)*d1
    # Create a surface SideNode at this intersection
    snode_sfc = SideNode(lon_int, lat_int, depth, var0)

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
        if len(data.shape) == 2:
            # Timeseries
            var1 = data[:,node1.id]
            var2 = data[:,node2.id]
        else:
            # Single time index
            var1 = data[node1.id]
            var2 = data[node2.id]
        var0 = var1 + (var2 - var1)/(d2 + d1)*d1
        # Similarly interpolate depth
        depth = node1.depth + (node2.depth - node1.depth)/(d2 + d1)*d1
        # Create a SideNode at this depth
        snode_below = SideNode(lon_int, lat_int, depth, var0)
        # Add to the linked list
        snode.set_below(snode_below)
        # Get ready for next iteration
        snode = snode_below

    return snode_sfc
    
                    

    

    
