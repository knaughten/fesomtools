from numpy import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from fesom_sidegrid import *

# Create a plot of a specified variable at a specified zonal slice, i.e. depth
# vs latitude.
# Input:
# mesh_path = path to directory containing grid files
# file_path = path to FESOM output file
# var_name = string containing name of variable in file_path
# tstep = int specifying index of time axis in file_path
# lon0 = longitude for zonal slice
# depth_min = deepest depth to plot (negative, in metres)
# save = optional boolean flag indicating whether to save plot to file
#        (otherwise will display on screen)
# fig_name = optional string containing name of figure file, if save=True
def zonal_slice_plot (mesh_path, file_path, var_name, tstep, lon0, depth_min, save=False, fig_name=None):

    # Set northern boundary and upper (surface) boundary
    lat_max = -50
    depth_max = 0
    # Font sizes for figure
    font_sizes = [18, 16, 12]

    # Read variable name and units for title
    id = Dataset(file_path, 'r')
    varid = id.variables[var_name]
    name = varid.getncattr('description')
    units = varid.getncattr('units')
    if lon0 < 0:
        lon_string = 'at ' + str(-lon0) + 'W'
    else:
        lon_string = 'at ' + str(lon0) + 'E'

    # Build the array of SideElements making up the zonal slice
    selements = fesom_sidegrid(mesh_path, file_path, var_name, tstep, lon0, lat_max)

    # Build an array of trapezoidal patches for the plot, and of data values
    # corresponding to each SideElement
    # Also find the minimum latitude of any SideElement
    patches = []
    values = []
    lat_min = lat_max
    for selm in selements:
        # Make patch
        coord = transpose(vstack((selm.y,selm.z)))
        patches.append(Polygon(coord, True, linewidth=0.))
        # Save data value
        values.append(selm.var)
        # Update minimum latitude if needed
        lat_min = min(lat_min, amin(selm.y))
    # Set southern boundary to be just south of the minimum latitude
    lat_min = lat_min-1

    # Set up plot
    fig = figure(figsize=(16,8))
    ax = fig.add_subplot(1,1,1)
    # Set colourmap for patches, and refer it to the values array
    img = PatchCollection(patches, cmap=jet)
    img.set_array(array(values))
    img.set_edgecolor('face')
    # Add patches to plot
    ax.add_collection(img)

    # Configure plot
    xlim(lat_min, lat_max)
    ylim(depth_min, depth_max)
    title(name + ' (' + units + ') ' + lon_string, fontsize=font_sizes[0])
    xlabel('Latitude', fontsize=font_sizes[1])
    ylabel('Depth (m)', fontsize=font_sizes[1])
    setp(ax.get_xticklabels(), fontsize=font_sizes[2])
    setp(ax.get_yticklabels(), fontsize=font_sizes[2])
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=font_sizes[2])

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()
