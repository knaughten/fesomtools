from numpy import *
from matplotlib.pyplot import *
from matplotlib.cm import *
from fesom_intersectgrid import *

# Create a plot of a specified variable zonally averaged between specified
# longitude bounds, i.e. depth vs. latitude.
# Input:
# mesh_path = path to directory containing grid files
# file_path = path to FESOM output file
# var_name = string containing name of variable in file_path
# tstep = int specifying index of time axis in file_path
# lon_min, lon_max = longitude bounds for zonal averages
# depth_min = deepest depth to plot (negative, in metres)
# save = optional boolean flag indicating whether to save plot to file
#        (otherwise will display on screen)
# fig_name = optional string containing name of figure file, if save=True
def zonal_avg_plot (mesh_path, file_path, var_name, tstep, lon_min, lon_max, depth_min, save=False, fig_name=None):

    # Set bounds on latitude for plot
    lat_min = -90
    lat_max = -50
    # Set upper (surface) boundary
    depth_max = 0
    # Set number of intervals for latitude and depth to interpolate to
    num_lat = 100
    num_depth = 50
    # Font sizes for figure
    font_sizes = [18, 16, 12]

    # Read variable name and units for title
    id = Dataset(file_path, 'r')
    varid = id.variables[var_name]
    name = varid.getncattr('description')
    units = varid.getncattr('units')
    if lon_min < 0:
        lon_min_string = 'averaged from ' + str(-lon_min) + 'W'
    else:
        lon_min_string = 'averaged from ' + str(lon_min) + 'E'
    if lon_max < 0:
        lon_max_string = ' to ' + str(-lon_max) + 'W'
    else:
        lon_max_string = ' to ' + str(lon_max) + 'E'

    # Build the array containing data on the regular grid lat_vals x depth_vals
    # Missing values will be NaN
    lat_vals, depth_vals, data_reg = fesom_intersectgrid(mesh_path, file_path, var_name, tstep, lon_min, lon_max, lat_min, lat_max, depth_min, depth_max, num_lat, num_depth)
    # Mask the NaNs
    data_reg = ma.masked_where(isnan(data_reg), data_reg)

    # Find the southernmost index with at least one unmasked value
    for j in range(num_lat):
        if sum(data_reg[:,j]) is not ma.masked:
            # This sum will only be masked if all entries are masked
            lat_min = lat_vals[j]
            break
    # Set southern boundary to be just south of the minimum latitude
    lat_min = lat_min-1

    # Set up plot
    fig = figure(figsize=(16,8))
    ax = fig.add_subplot(1,1,1)
    # Shade interpolated data
    img = ax.pcolorfast(lat_vals, depth_vals, data_reg)

    # Configure plot
    xlim(lat_min, lat_max)
    ylim(depth_min, depth_max)
    title(name + ' (' + units + '), ' + lon_min_string + lon_max_string, fontsize=font_sizes[0])
    xlabel('Latitude', fontsize=font_sizes[1])
    ylabel('Depth (m)', fontsize=font_sizes[1])
    setp(ax.get_xticklabels(), fontsize=font_sizes[2])
    setp(ax.get_yticklabels(), fontsize=font_sizes[2])
    cbar = colorbar(img, ax=ax)
    cbar.ax.tick_params(labelsize=font_sizes[2])

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()
