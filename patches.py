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
def make_patches (mesh_path, circumpolar=False, mask_cavities=False, only_major=False):

    if only_major:
        lon_min = [-62.67, -65.5, -79.17, -85, -104.17, -102.5, -108.33, -114.5, -135.67, -149.17, -155, 144, 115, 94.17, 80.83, 65, 33.83, 19, 12.9, 9.33, -10.05, -28.33, -181, 158.33]
        lon_max = [-59.33, -60, -66.67, -28.33, -88.83, -99.17, -103.33, -111.5, -114.33, -140, -145, 146.62, 123.33, 102.5, 89.17, 75, 37.67, 33.33, 16.17, 12.88, 7.6, -10.33, -146.67, 181]
        lat_min = [-73.03, -69.35, -74.17, -83.5, -73.28, -75.5, -75.5, -75.33, -74.9, -76.42, -78, -67.83, -67.17, -66.67, -67.83, -73.67, -69.83, -71.67, -70.5, -70.75, -71.83, -76.33, -85, -84.5]
        lat_max = [-69.37, -66.13, -69.5, -74.67, -71.67, -74.17, -74.67, -73.67, -73, -75.17, -76.41, -66.67, -66.5, -64.83, -66.17, -68.33, -68.67, -68.33, -69.33, -69.83, -69.33, -71.5, -77.77, -77]
        num_shelves = len(lon_min)-1

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
                if only_major:
                    # Only exclude major ice shelves
                    keep = True
                    for index in range(num_shelves):
                        if all(elm.lon >= lon_min[index]) and all(elm.lon <= lon_max[index]) and all(elm.lat >= lat_min[index]) and all(elm.lat <= lat_max[index]):
                            keep = False
                        if index == num_shelves-1:
                            if all(elm.lon >= lon_min[index+1]) and all(elm.lon <= lon_max[index+1]) and all(elm.lat >= lat_min[index+1]) and all(elm.lat <= lat_max[index+1]):
                                keep = False
                    if keep:
                        coord = transpose(vstack((elm.x, elm.y)))
                        patches.append(Polygon(coord, True, linewidth=0.))
        else:
            coord = transpose(vstack((elm.x, elm.y)))
            patches.append(Polygon(coord, True, linewidth=0.))

    return elements, patches


# Given the FESOM grid elements, select only the ones in ice shelf cavities
# and build a separate array of patches for use in masking the ice shelves.
def iceshelf_mask (elements, only_major=False):

    if only_major:
        lon_min = [-62.67, -65.5, -79.17, -85, -104.17, -102.5, -108.33, -114.5, -135.67, -149.17, -155, 144, 115, 94.17, 80.83, 65, 33.83, 19, 12.9, 9.33, -10.05, -28.33, -181, 158.33]
        lon_max = [-59.33, -60, -66.67, -28.33, -88.83, -99.17, -103.33, -111.5, -114.33, -140, -145, 146.62, 123.33, 102.5, 89.17, 75, 37.67, 33.33, 16.17, 12.88, 7.6, -10.33, -146.67, 181]
        lat_min = [-73.03, -69.35, -74.17, -83.5, -73.28, -75.5, -75.5, -75.33, -74.9, -76.42, -78, -67.83, -67.17, -66.67, -67.83, -73.67, -69.83, -71.67, -70.5, -70.75, -71.83, -76.33, -85, -84.5]
        lat_max = [-69.37, -66.13, -69.5, -74.67, -71.67, -74.17, -74.67, -73.67, -73, -75.17, -76.41, -66.67, -66.5, -64.83, -66.17, -68.33, -68.67, -68.33, -69.33, -69.83, -69.33, -71.5, -77.77, -77]
        num_shelves = len(lon_min)-1

    mask_patches = []
    for elm in elements:
        if elm.cavity:
            if only_major:
                # Only include major ice shelves
                keep = False
                for index in range(num_shelves):
                    if all(elm.lon >= lon_min[index]) and all(elm.lon <= lon_max[index]) and all(elm.lat >= lat_min[index]) and all(elm.lat <= lat_max[index]):
                        keep = True
                    if index == num_shelves-1:
                        if all(elm.lon >= lon_min[index+1]) and all(elm.lon <= lon_max[index+1]) and all(elm.lat >= lat_min[index+1]) and all(elm.lat <= lat_max[index+1]):
                            keep = True
                if keep:
                    coord = transpose(vstack((elm.x, elm.y)))
                    mask_patches.append(Polygon(coord, True, linewidth=0.))   
            else:
                coord = transpose(vstack((elm.x, elm.y)))
                mask_patches.append(Polygon(coord, True, linewidth=0.))

    return mask_patches
