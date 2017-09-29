from numpy import *
from matplotlib.patches import Polygon
from fesom_sidegrid import *

# Linearly interpolate FESOM data to the specified longitude.
# Input:
# elements = array of 2D Elements created by fesom_grid.py
# lat_max = maximum latitude to consider
# lon0 = longitude to interpolate to, from -180 to 180
# data = array of FESOM data on original mesh
def side_patches (elements, lat_max, lon0, data):

    # Get SideElements interpolated to lon0
    selements = fesom_sidegrid(elements, data, lon0, lat_max)
    # Build an array of quadrilateral patches for the plot, and of data
    # values corresponding to each SideElement
    patches = []
    values = []
    lat_min = lat_max
    for selm in selements:
        # Make patch
        coord = transpose(vstack((selm.y,selm.z)))
        patches.append(Polygon(coord, True, linewidth=0.))
        # Save data value
        values.append(selm.var)
        lat_min = min(lat_min, amin(selm.y))
    # Show a little bit of the land mask
    lat_min = lat_min-0.5

    return patches, values, lat_min
