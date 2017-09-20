from numpy import *
from triangle_area import *
from unrotate_grid import *

# Classes and routines to build an object-oriented 2D FESOM grid data structure

# Node object containing node id number, longitude, latitude
class Node:

    # Initialise with location data
    def __init__ (self, id, lon, lat):

        self.id = id
        self.lon = lon
        self.lat = lat


# Element object containing the three nodes making up the triangular element.
class Element:

    # Initialise
    # Input: node1, node2, node3 = component Nodes
    # cavity_nodes = array of booleans indicating whether each node is in an
    #                ice shelf cavity
    # circumpolar = boolean flag indicating whether the user's plot will be
    #               circumpolar Antarctic (True) or global (False)
    # is_repeat = optional boolean flag indicating if this Element has already
    #             been processed once (Elements crossing the line 180W=180E are
    #             processed twice)
    def __init__ (self, node1, node2, node3, cavity_nodes, circumpolar, is_repeat=False):

        self.nodes = array([node1, node2, node3])
        lon = array([node1.lon, node2.lon, node3.lon])
        lat = array([node1.lat, node2.lat, node3.lat])
        self.lon = lon
        self.lat = lat
        self.repeat_next = False
        self.cavity_nodes = cavity_nodes
        self.cavity = all(cavity_nodes)

        # Check for elements which cross longitude 180W = 180E
        # We want a copy of these elements on both sides of this line
        if (abs(lon[2]-lon[1]) > 170) or (abs(lon[1]-lon[0]) > 170) or \
                (abs(lon[2]-lon[0]) > 170):
          if is_repeat:
            # Already processed this element in the eastern hemisphere
            # Process it in the western hemisphere now
            index = nonzero(lon > 0)
            self.lon[index] = lon[index] - 360
          else:
            # Haven't seen this element before
            # Process it in the eastern hemisphere
            index = nonzero(lon < 0)
            lon[index] = lon[index] + 360
            # Flag to process it again in the western hemisphere
            self.repeat_next = True
        if circumpolar:
          # Convert to Cartesian coordinates for plotting
          deg2rad = pi/180
          self.x = -(lat + 90)*cos(lon*deg2rad+pi/2)
          self.y = (lat + 90)*sin(lon*deg2rad+pi/2)
        else:
          self.x = lon
          self.y = lat

    # Return the area of the triangle making up this Element
    def area (self):
        return triangle_area(self.lon, self.lat)


# Function to read FESOM grid files and build Element mesh
# Input: 
# mesh_path = path to directory containing grid files
# circumpolar = optional boolean flag indicating if the user's plot
#               will be circumpolar Antarctic (otherwise global)
# cross_180 = optional boolean flag indicating that elements which cross the
#             line 180W=180E should be copied to both sides of the line
#             (this is desirable for plotting but not for calculations like
#             global integrals). Default True.
# Output:
# elements = array of Element objects
def fesom_grid_2d (mesh_path, circumpolar=False, cross_180=True):

    # Northern boundary of circumpolar Antarctic domain
    nbdry = -30

    # Read 2D node information
    file = open(mesh_path + 'nod2d.out', 'r')
    file.readline()
    rlon = []
    rlat = []
    for line in file:
        tmp = line.split()
        lon_tmp = float(tmp[1])
        lat_tmp = float(tmp[2])
        # Make sure longitude is in the range [-180, 180]
        if lon_tmp < -180:
            lon_tmp = lon_tmp + 360
        elif lon_tmp > 180:
            lon_tmp = lon_tmp - 360
        rlon.append(lon_tmp)
        rlat.append(lat_tmp)
    file.close()
    rlon = array(rlon)
    rlat = array(rlat)
    # Unrotate grid
    lon, lat = unrotate_grid(rlon, rlat)

    # Create Nodes from location data
    nodes = []
    for i in range(size(lon)):
        nodes.append(Node(i, lon[i], lat[i]))

    # Read 2D cavity flag
    file = open(mesh_path + 'cavity_flag_nod2d.out', 'r')
    cavity = []
    for line in file:
        tmp = int(line)
        if tmp == 1:
            cavity.append(True)
        else:
            cavity.append(False)
    file.close()

    # Read 2D elements (triangles of connecting nodes)
    file = open(mesh_path + 'elem2d.out', 'r')
    file.readline()
    elements = []
    for line in file:
        tmp = line.split()
        id1 = int(tmp[0])-1
        id2 = int(tmp[1])-1
        id3 = int(tmp[2])-1
        cavity_nodes = [cavity[id1], cavity[id2], cavity[id3]]
        # Initialise the Element
        elm = Element(nodes[id1], nodes[id2], nodes[id3], cavity_nodes, circumpolar)
        if circumpolar:
            # Only save elm if it's within the circumpolar Antarctic domain
            if elm.nodes[0].lat < nbdry or elm.nodes[1].lat < nbdry or elm.nodes[2].lat < nbdry:
                elements.append(elm)
        else:
            elements.append(elm)

        if cross_180 and elm.repeat_next:
            # This element crosses the line of longitude 180W = 180E
            # Process it a second time so it shows up on both sides of the line
            elm_rep = Element(nodes[id1], nodes[id2], nodes[id3], cavity_nodes, circumpolar, True)
            if circumpolar:
                # Only save elm if it's within the circumpolar Antarctic domain
                if elm.nodes[0].lat < nbdry or elm.nodes[1].lat < nbdry or elm.nodes[2].lat < nbdry:
                    elements.append(elm_rep)
            else:
                elements.append(elm_rep)
    file.close()

    return elements

