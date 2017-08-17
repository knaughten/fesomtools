from numpy import *
from triangle_area import *
from unrotate_grid import *

# Classes and routines to build an object-oriented FESOM grid data structure


# Node object containing node id number, longitude, latitude, depth,
# and the node directly below it (if there is one).
class Node:

  # Initialise with location data
  def __init__ (self, id, lon, lat, depth):

    self.id = id
    self.lon = lon
    self.lat = lat
    self.depth = depth
    self.below = None


  # Save the node directly below the current node
  def set_below (self, node_below):

    self.below = node_below
  

  # Travel straight down the water column to find the node at the seafloor.
  def find_bottom (self):

    node = self
    id = node.id

    while node.below is not None:
      # Go down another level
      node = node.below

    return node


  # Given a depth z, travel straight down the water column to find nodes
  # A and B such that depth(A) <= z <= depth (B) and coefficients c and d
  # such that c*depth(A) + d*depth(B) = z. Return id(A), id(B), c, and d.
  # In practice this function will only be called for nodes at the surface.
  def find_depth (self, z):

    node_above = self
    node_below = node_above.below
    id1 = id2 = NaN
    coeff1 = coeff2 = NaN

    while True:
      if node_above.depth > z:
        # We're already too deep
        # This might happen if there is an ice shelf in this location and
        # z is less than the ice shelf draft
        # Return NaNs for everything
        break
      elif node_below.below is None:
        # We've reached the seafloor
        # Return NaNs for everything
        break
      elif node_below.depth >= z:
        # z is between node_above and node_below
        # Save the ids and calculate the coefficients
        id1 = node_above.id
        id2 = node_below.id
        coeff1 = (node_below.depth - z)/(node_below.depth - node_above.depth)
        coeff2 = 1 - coeff1
        break
      else:
        # We're still too shallow
        # Go one level down for the next iteration of the loop
        node_above = node_below
        node_below = node_above.below  
    
    return id1, id2, coeff1, coeff2    
           


# Element object containing the three nodes making up the triangular element.
# In practice this will only be called for surface nodes, but you can travel
# straight down to any depth by calling Node.find_depth or Node.find_bottom.
class Element:

  # Initialise
  # Input: node1, node2, node3 = component Nodes
  # cavity_nodes = array of booleans indicating whether each node is in an
  #                ice shelf cavity
  # coast_nodes = array of booleans indicating whether each node is coastal
  # circumpolar = boolean flag indicating whether the user's plot will be
  #               circumpolar Antarctic (True) or global (False)
  # is_repeat = optional boolean flag indicating if this Element has already
  #             been processed once (Elements crossing the line 180W=180E are
  #             processed twice)
  def __init__ (self, node1, node2, node3, cavity_nodes, coast_nodes, circumpolar, is_repeat=False):
    
    self.nodes = array([node1, node2, node3])
    lon = array([node1.lon, node2.lon, node3.lon])
    lat = array([node1.lat, node2.lat, node3.lat])
    self.lon = lon
    self.lat = lat
    self.repeat_next = False
    self.cavity_nodes = cavity_nodes
    self.coast_nodes = coast_nodes
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
def fesom_grid (mesh_path, circumpolar=False, cross_180=True):

  # Northern boundary of circumpolar Antarctic domain
  nbdry = -30

  # Read 3d node information
  file = open(mesh_path + 'nod3d.out', 'r')
  file.readline()
  rlon3d = []
  rlat3d = []
  depth3d = []

  for line in file:
      tmp = line.split()
      lon = float(tmp[1])
      lat = float(tmp[2])
      depth = float(tmp[3])
      rank = int(tmp[4])
      # Make sure longitude is in the range [-180, 180]
      if lon < -180:
          lon = lon + 360
      elif lon > 180:
          lon = lon - 360
      rlon3d.append(lon)
      rlat3d.append(lat)
      depth3d.append(depth)
  file.close()
  rlon3d = array(rlon3d)
  rlat3d = array(rlat3d)

  # Unrotate grid
  lon3d, lat3d = unrotate_grid(rlon3d, rlat3d)

  # Create Nodes from location data
  nodes = []
  for i in range(size(lon3d)):
    nodes.append(Node(i, lon3d[i], lat3d[i], -depth3d[i]))

  # Read lists of which nodes are directly below which
  file = open(mesh_path + 'aux3d.out', 'r')
  file.readline()
  # ID of the surface node in the current column
  sfc_id = int(file.readline())
  # Current node in the current column
  curr_node = nodes[sfc_id-1]

  for line in file:
    tmp = int(line)
    # -999 means there is no node below
    if tmp != -999:
      if tmp == sfc_id+1:
        # We are onto the next water column
        # Save the new sfc_id
        sfc_id = tmp
      elif tmp == 1:
        # We are finished all of the nodes
        break
      else:
        # Save the node below
        curr_node.set_below(nodes[tmp-1])
      # Save the new curr_node
      curr_node = nodes[tmp-1]

  # Read 2D cavity flag
  file = open(mesh_path + 'cavity_flag_nod2d.out', 'r')
  cavity = []

  for line in file:
    tmp = int(line)
    if tmp == 1:
      cavity.append(True)
    elif tmp == 0:
      cavity.append(False)
    else:
      print 'Problem'
  file.close()

  # Read coast flag
  file = open(mesh_path + 'nod2d.out', 'r')
  file.readline()
  coast = []

  for line in file:
    tmp = line.split()
    coast_tmp = int(tmp[3])
    if coast_tmp == 1:
      coast.append(True)
    elif coast_tmp == 0:
      coast.append(False)
    else:
      print 'Problem'
  file.close()

  # Read 2D elements (triangles of 3 connecting nodes)
  file = open(mesh_path + 'elem2d.out', 'r')
  file.readline()
  elements = []

  for line in file:

    tmp = line.split()
    id1 = int(tmp[0])-1
    id2 = int(tmp[1])-1
    id3 = int(tmp[2])-1
    # Make arrays of cavity and coast flags
    cavity_nodes = [cavity[id1], cavity[id2], cavity[id3]]
    coast_nodes = [coast[id1], coast[id2], coast[id3]]
    # Initialise the Element
    elm = Element(nodes[id1], nodes[id2], nodes[id3], cavity_nodes, coast_nodes, circumpolar)

    if circumpolar:
      # Only save elm if it's within the circumpolar Antarctic domain
      if elm.nodes[0].lat < nbdry or elm.nodes[1].lat < nbdry or elm.nodes[2].lat < nbdry:
        elements.append(elm)
    else:
      elements.append(elm)

    if cross_180 and elm.repeat_next:
      # This element crosses the line of longitude 180W = 180E
      # Process it a second time so it shows up on both sides of the line
      elm_rep = Element(nodes[id1], nodes[id2], nodes[id3], cavity_nodes, coast_nodes, circumpolar, True)
      if circumpolar:
        # Only save elm if it's within the circumpolar Antarctic domain
        if elm.nodes[0].lat < nbdry or elm.nodes[1].lat < nbdry or elm.nodes[2].lat < nbdry:
          elements.append(elm_rep)
      else:
        elements.append(elm_rep)
  file.close()

  return elements

  
