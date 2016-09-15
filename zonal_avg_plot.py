from numpy import *
from netCDF4 import Dataset
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
# set_limits = optional boolean flag indicating whether or not to set manual
#              limits on colourbar (otherwise limits determined automatically)
# limits = optional array containing min and max limits
def zonal_avg_plot (mesh_path, file_path, var_name, tstep, lon_min, lon_max, depth_min, save=False, fig_name=None, set_limits=False, limits=None):

    # Set bounds on latitude for plot
    lat_min = -90
    lat_max = -30
    # Set upper (surface) boundary
    depth_max = 0
    # Set number of intervals for latitude and depth to interpolate to
    num_lat = 100
    num_depth = 50
    # Font sizes for figure
    font_sizes = [30, 24, 20]

    # Read variable name
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

    # Choose colour bounds
    if set_limits:
        # User-specified bounds
        var_min = limits[0]
        var_max = limits[1]
        if var_min == -var_max:
            # Bounds are centered on zero, so choose a blue-to-red colourmap
            # centered on yellow
            colour_map = 'RdYlBu_r'
        else:
            colour_map = 'jet'
    else:
        # Determine bounds automatically
        if var_name in ['u', 'v', 'w']:
            # Center levels on 0 for certain variables, with a blue-to-red
            # colourmap
            max_val = amax(abs(data_reg))
            var_min = -max_val
            var_max = max_val
            colour_map = 'RdYlBu_r'
        else:
            var_min = amin(data_reg)
            var_max = amax(data_reg)
            colour_map = 'jet'

    # Set up plot
    fig = figure(figsize=(16,8))
    ax = fig.add_subplot(1,1,1)
    # Shade interpolated data
    img = ax.pcolorfast(lat_vals, depth_vals, data_reg, cmap=colour_map)

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
    img.set_clim(vmin=var_min, vmax=var_max)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()
