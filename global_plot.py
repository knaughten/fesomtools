from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from patches import *


# Create a plot of a specified variable on a global domain
# Input:
# file_path = string containing path to FESOM output file
# var_name = string containing name of variable in file_path
# depth_key = int specifying whether to plot surface nodes (0), bottom nodes
#             (1), or some specified depth (2)
# depth = if depth_key==2, specified depth in m
# tstep = int specifying index of time axis in file_path
# elements = array of Elements for the global grid (created using fesom_grid)
# patches = array of Polygon patches corresponding to elements (created using
#           make_patches)
# mask_cavities = optional boolean flag indicating whether to mask ice shelf
#                 cavities in grey
# save = optional boolean flag indicating whether to save plot to file
#        (otherwise will display on screen)
# fig_name = optional string containing name of figure file, if save = True
def global_plot (file_path, var_name, depth_key, depth, tstep, elements, patches, mask_cavities=False, save=False, fig_name=None):

    # Set bounds for domain
    lon_min = -180
    lon_max = 180
    lat_min = -90
    lat_max = 90
    # Configure position of latitude and longitude labels
    lon_ticks = arange(-120, 120+1, 60)
    lat_ticks = arange(-60, 60+1, 30)
    # Set font sizes
    font_sizes = [18, 16, 12]

    # Read data
    file = Dataset(file_path, 'r')
    varid = file.variables[var_name]
    data = varid[tstep-1,:]
    # Set descriptive variable name and units for title
    if var_name == 'area':
        name = 'ice concentration'
        units = 'fraction'
    else:
        name = varid.getncattr('description')
        units = varid.getncattr('units')
    if depth_key == 0:
        if var_name in ['temp', 'salt', 'u', 'v']:
            depth_string = 'at surface'
        else:
            depth_string = ''
    elif depth_key == 1:
        depth_string = 'at bottom'
    elif depth_key == 2:
        depth_string = 'at ' + str(depth) + ' m'

    # Build an array of data values corresponding to each Element
    values = []
    plot_patches = []
    for elm in elements:
        # If mask_cavities is true, only include elements which are not inside
        # an ice shelf cavity; otherwise, include all elements
        if (mask_cavities and not elm.cavity) or (not mask_cavities):
            if depth_key == 0:
                # Surface nodes; this is easy
                # Average the data value for each of the three component nodes
                values.append(mean([data[elm.nodes[0].id], data[elm.nodes[1].id], data[elm.nodes[2].id]]))
            elif depth_key == 1:
                # Bottom nodes
                values_tmp = []
                # For each of the three component nodes, find the id of the
                # bottom node beneath it
                for i in range(3):
                    id = elm.nodes[i].find_bottom()
                    values_tmp.append(data[id])
                # Average over these three values
                values.append(mean(values_tmp))
            elif depth_key == 2:
                # Specified depth
                values_tmp = []
                # For each of the three component nodes, linearly interpolate
                # to the correct depth
                for i in range(3):
                    # Find the ids of the nodes above and below this depth,
                    # and the coefficients for the linear interpolation
                    id1, id2, coeff1, coeff2 = elm.nodes[i].find_depth(depth)
                    if any(isnan(array([id1, id2, coeff1, coeff2]))):
                        # No nodes at this depth, so save NaN as placeholder
                        values_tmp.append(NaN)
                    else:
                        values_tmp.append(coeff1*data[id1] + coeff2*data[id2])
                if any (isnan(array(values_tmp))):
                    pass
                else:
                    values.append(mean(values_tmp))
                    coord = transpose(vstack((elm.x, elm.y)))
                    # Make new patches for elements which exist at this depth
                    plot_patches.append(Polygon(coord, True, linewidth=0.))

    if depth_key != 2:
        # Use all patches
        plot_patches = patches[:]

    # Set up plot
    fig = figure(figsize=(16,8))
    ax = fig.add_subplot(1,1,1)
    # Set colourmap for patches, and refer it to the values array
    img = PatchCollection(plot_patches, cmap=jet)
    img.set_array(array(values))
    img.set_edgecolor('face')
    # Add patches to plot
    ax.add_collection(img)
    if mask_cavities:
        # Get mask array of patches for ice shelf cavity elements
        mask_patches = iceshelf_mask(elements)
        # Set colour to light grey for all patches in the mask
        overlay = PatchCollection(mask_patches, facecolor=(0.6, 0.6, 0.6))
        overlay.set_edgecolor('face')
        # Add mask to plot
        ax.add_collection(overlay)

    # Configure plot
    xlim(lon_min, lon_max)
    ylim(lat_min, lat_max)
    xticks(lon_ticks)
    yticks(lat_ticks)
    title(name + ' (' + units + ') ' + depth_string, fontsize=font_sizes[0])    
    xlabel('Longitude', fontsize=font_sizes[1])
    ylabel('Latitude', fontsize=font_sizes[1])
    setp(ax.get_xticklabels(), fontsize=font_sizes[2])
    setp(ax.get_yticklabels(), fontsize=font_sizes[2])
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=font_sizes[2])

    if save:
        savefig(fig_name)
    else:
        show()
