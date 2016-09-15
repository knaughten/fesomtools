from numpy import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from unrotate_vector import *
from fesom_grid import *
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
# set_limits = optional boolean flag indicating whether or not to set manual
#              limits on colourbar (otherwise limits determined automatically)
# limits = optional array containing min and max limits
def zonal_slice_plot (mesh_path, file_path, var_name, tstep, lon0, depth_min, save=False, fig_name=None, set_limits=False, limits=None):

    # Set northern boundary and upper (surface) boundary
    lat_max = -50
    depth_max = 0
    # Font sizes for figure
    font_sizes = [30, 24, 20]

    # Read variable name and units for title
    id = Dataset(file_path, 'r')
    varid = id.variables[var_name]
    name = varid.getncattr('description')
    units = varid.getncattr('units')
    if lon0 < 0:
        lon_string = 'at ' + str(-lon0) + 'W'
    else:
        lon_string = 'at ' + str(lon0) + 'E'

    # Read data
    data = id.variables[var_name][tstep-1,:]
    # Check for vector variables that need to be unrotated
    if var_name in ['u', 'v']:
        # Read the rotated lat and lon
        fid = open(mesh_path + 'nod3d.out', 'r')
        fid.readline()
        lon = []
        lat = []
        for line in fid:
            tmp = line.split()
            lon_tmp = float(tmp[1])
            lat_tmp = float(tmp[2])
            if lon_tmp < -180:
                lon_tmp += 360
            elif lon_tmp > 180:
                lon_tmp -= 360
            lon.append(lon_tmp)
            lat.append(lat_tmp)
        fid.close()
        lon = array(lon)
        lat = array(lat)
        if var_name == 'u':
            u_data = data[:]
            v_data = id.variables['v'][tstep-1,:]
            u_data_lonlat, v_data_lonlat = unrotate_vector(lon, lat, u_data, v_data)
            data = u_data_lonlat[:]
        elif var_name == 'v':
            v_data = data[:]
            u_data = id.variables['u'][tstep-1,:]
            u_data_lonlat, v_data_lonlat = unrotate_vector(lon, lat, u_data, v_data)
            data = v_data_lonlat[:]
    id.close()    

    # Build the regular FESOM grid
    elm2D = fesom_grid(mesh_path)

    # Build the array of SideElements making up the zonal slice
    selements = fesom_sidegrid(elm2D, data, lon0, lat_max)

    # Build an array of quadrilateral patches for the plot, and of data values
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
            max_val = amax(abs(array(values)))
            var_min = -max_val
            var_max = max_val
            colour_map = 'RdYlBu_r'
        else:
            var_min = amin(array(values))
            var_max = amax(array(values))
            colour_map = 'jet'

    # Set up plot
    fig = figure(figsize=(16,8))
    ax = fig.add_subplot(1,1,1)
    # Set colourmap for patches, and refer it to the values array
    img = PatchCollection(patches, cmap=colour_map)
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
    img.set_clim(vmin=var_min, vmax=var_max)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()
