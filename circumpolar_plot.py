from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from patches import *


# Create a plot of a specified variable on a circumpolar Antarctic domain
# Input:
# file_path = string containing path to FESOM output file
# var_name = string containing name of variable in file_path
# tstep = int specifying index of time axis in file_path
# elements = array of Elements for the global grid (created using fesom_grid)
# patches = array of Polygon patches corresponding to elements (created using
#           make_patches)
# mask_cavities = optional boolean flag indicating whether to mask ice shelf
#                 cavities in grey
# save = optional boolean flag indicating whether to save plot to file
#        (otherwise will display on screen)
# fig_name = optional string containing name of figure file, if save = True
def circumpolar_plot (file_path, var_name, tstep, elements, patches, mask_cavities=False, save=False, fig_name=None):

    # Northern boundary 60S
    lat_max = -60+90
    # Set font sizes
    font_sizes = [18, 16, 12]
    # Seconds per year, for conversion of ice shelf melt rate
    sec_per_year = 365.25*24*3600

    # Read data
    file = Dataset(file_path, 'r')
    varid = file.variables[var_name]
    data = varid[tstep-1,:]
    # Set descriptive variable name and units for title
    if var_name == 'area':
        name = 'ice concentration'
        units = 'fraction'
    elif var_name == 'wnet':
        # Convert from m/s to m/y
        data = data*sec_per_year
        name = 'ice shelf melt rate'
        units = 'm/y'
    else:
        name = varid.getncattr('description')
        units = varid.getncattr('units')

    # Build an array of data values corresponding to each Element
    values = []
    for elm in elements:
        if mask_cavities:
            # Only include Elements not beneath ice shelf cavities
            if elm.cavity == False:
                # Average data for each of the three component nodes
                values.append(mean([data[elm.nodes[0].id], data[elm.nodes[1].id], data[elm.nodes[2].id]]))
        else:
            values.append(mean([data[elm.nodes[0].id], data[elm.nodes[1].id], data[elm.nodes[2].id]]))

    if mask_cavities:
        # Get mask array of patches for ice shelf cavity elements
        mask_patches = iceshelf_mask(elements)
        if var_name == 'wnet':
            # Swap with regular patches so that open ocean elements are masked,
            # ice shelf cavity nodes are not
            tmp = patches
            patches = mask_patches
            mask_patches = tmp            

    # Set up plot
    fig = figure(figsize=(16,12))
    ax = fig.add_subplot(1,1,1, aspect='equal')
    # Set colourmap for patches, and refer it to the values array
    img = PatchCollection(patches, cmap=jet)
    img.set_array(array(values))
    img.set_edgecolor('face')
    # Add patches to plot
    ax.add_collection(img)
    if mask_cavities:
        # Set colour to light grey for patches in mask
        overlay = PatchCollection(mask_patches, facecolor=(0.6, 0.6, 0.6))
        overlay.set_edgecolor('face')
        # Add mask to plot
        ax.add_collection(overlay)

    # Configure plot
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    title(name + ' (' + units + ')', fontsize=font_sizes[0])
    ax.get_xaxis().set_ticks([]
)
    ax.get_yaxis().set_ticks([])
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=font_sizes[2])

    if save:
        savefig(name)
    else:
        show()
