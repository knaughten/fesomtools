from numpy import *
from matplotlib.patches import Polygon
from fesom_grid import *


# Create the FESOM grid Elements and convert them to triangular patches for
# the user's plot.
# Input:
# mesh_path = path to directory containing grid files
# circumpolar = optional boolean flag indicating if the user's plot will be
#               circumpolar Antarctic (otherwise global)
# mask_cavities = optional boolean flag indicating if Elements in ice shelf
#                 cavities should be excluded
# Output:
# elements = array of Element objects which make up the FESOM 2D mesh
# patches = array of triangular Polygon objects to be used on the plot
def make_patches (mesh_path, circumpolar=False, mask_cavities=False):

    # Read the grid and build the Element array
    elements = fesom_grid(mesh_path, circumpolar)

    patches = []
    for elm in elements:
        if mask_cavities:
            # Only build patches for Elements not in ice shelf cavities
            if elm.cavity == False:
                coord = transpose(vstack((elm.x, elm.y)))
                patches.append(Polygon(coord, True, linewidth=0.))
        else:
            coord = transpose(vstack((elm.x, elm.y)))
            patches.append(Polygon(coord, True, linewidth=0.))

    return elements, patches


# Given the FESOM grid elements, select only the ones in ice shelf cavities
# and build a separate array of patches for use in masking the ice shelves.
def iceshelf_mask (elements):

    mask_patches = []
    for elm in elements:
        if elm.cavity:
            coord = transpose(vstack((elm.x, elm.y)))
            mask_patches.append(Polygon(coord, True, linewidth=0.))

    return mask_patches
