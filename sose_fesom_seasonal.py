from netCDF4 import Dataset
from numpy import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from fesom_grid import *
from seasonal_avg import *
from side_patches import *

# Make a 4x2 plot comparing lat vs. depth slices of seasonally averaged
# temperature or salinity at the given longitude, between FESOM (given year
# of simulation) and SOSE (2005-2010 climatology).
# Input:
# elements = array of 2D Elements created using fesom_grid.py
# file_path1 = path to a FESOM output file containing one year of 5-day
#              averages for ocean variables (we will just use December)
# file_path2 = path to a FESOM output file containing the following year of
#              5-day averages for ocean variables (we will use January
#              through November)
# var_name = 'temp' for temperature or 'salt' for salinity
# lon0 = the specific longitude to plot (between -180 and 180)
# depth_min = deepest depth to plot (negative, in m)
# save = optional boolean flag; if True, the figure will be saved with file name
#        fig_name; if False, the figure will display on the screen
# fig_name = optional string containing filename for figure, if save=True
def sose_fesom_seasonal (elements, file_path1, file_path2, var_name, lon0, depth_min, save=False, fig_name=None):

    # Path to SOSE seasonal climatology file
    sose_file = '/short/m68/kaa561/SOSE_seasonal_climatology.nc'
    lat_max = -60 #-30
    season_names = ['DJF', 'MAM', 'JJA', 'SON']

    # Bounds on colour scale
    if var_name == 'temp':
        var_min = -2.5
        var_max = 3.5 #7.5
        var_ticks = 1
    elif var_name == 'salt':
        var_min = 33.8
        var_max = 34.8
        var_ticks = 0.2
    else:
        print 'Unknown variable ' + var_name
        return

    # Choose what to write on the title about the variable
    if var_name == 'temp':
        var_string = r'Temperature ($^{\circ}$C)'
    elif var_name == 'salt':
        var_string = 'Salinity (psu)'
    # Choose what to write on the title about longitude
    if lon0 < 0:
        lon_string = ' at ' + str(int(round(-lon0))) + r'$^{\circ}$W'
    else:
        lon_string = ' at ' + str(int(round(lon0))) + r'$^{\circ}$E'

    print 'Processing SOSE data'
    # Read grid and 3D data (already seasonally averaged)
    id = Dataset(sose_file, 'r')
    lon_sose = id.variables['longitude'][0,:]
    lat_sose = id.variables['latitude'][:,0]
    z_sose = id.variables['depth'][:]
    var_3d_sose = id.variables[var_name][:,:,:,:]

    # Calculate zonal slices for each season
    var_sose = ma.empty([4, size(z_sose), size(lat_sose,0)])
    var_sose[:,:,:] = 0.0
    for season in range(4):
        print 'Calculating zonal slices for ' + season_names[season]
        var_sose[season,:,:] = interp_lon_sose(var_3d_sose[season,:,:,:], lon_sose, lon0)

    # Get seasonal averages of the FESOM output
    fesom_data = seasonal_avg(file_path1, file_path2, var_name)

    # Set colour levels
    lev = linspace(var_min, var_max, num=50)

    # Choose southern boundary based on extent of SOSE grid
    lat_min = amin(lat_sose)

    # Plot
    fig = figure(figsize=(20,9))
    for season in range(4):
        # FESOM
        print 'Calculating zonal slices for ' + season_names[season]
        patches, values, tmp = side_patches(elements, lat_max, lon0, fesom_data[season,:])
        ax = fig.add_subplot(2, 4, season+1)
        img = PatchCollection(patches, cmap=jet)
        img.set_array(array(values))
        img.set_edgecolor('face')
        img.set_clim(vmin=var_min, vmax=var_max)
        ax.add_collection(img)
        xlim([lat_min, lat_max])
        ylim([depth_min, 0])
        title('FESOM (' + season_names[season] + ')', fontsize=24)
        if season == 0:
            ylabel('depth (m)', fontsize=18)
        # SOSE
        fig.add_subplot(2, 4, season+5)
        pcolormesh(lat_sose, z_sose, var_sose[season,:,:], vmin=var_min, vmax=var_max, cmap='jet')
        xlim([lat_min, lat_max])
        ylim([depth_min, 0])
        title('SOSE (' + season_names[season] + ')', fontsize=24)
        xlabel('Latitude', fontsize=18)
        if season == 0:
            ylabel('depth (m)', fontsize=18)
    # Add colorbar
    cbaxes = fig.add_axes([0.93, 0.2, 0.015, 0.6])
    cbar = colorbar(img, cax=cbaxes, ticks=arange(var_min, var_max+var_ticks, var_ticks))
    cbar.ax.tick_params(labelsize=16)
    # Add the main title
    suptitle(var_string + lon_string, fontsize=30)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()
    

# Linearly interpolate SOSE data to the specified longitude.
# Input:
# data_3d = arary of data, dimension depth x lat x lon
# lon = 1D array of longitude values (between 0 and 360)
# lon0 = longitude to interpolate to (between -180 and 180)
# Output:
# data = array of data interpolated to lon0, dimension depth x lat
def interp_lon_sose (data_3d, lon, lon0):

    # Convert lon0 to be between 0 and 360
    if lon0 < 0:
        lon0 += 360

    if lon0 < lon[0] or lon0 > lon[-1]:
        # Special case: lon0 on periodic boundary
        # Be careful with mod 360 here
        iw = size(lon)-1
        ie = 0
        # Calculate difference between lon0 and lon[iw], mod 360 if necessary
        dlon_num = lon0 - lon[iw]
        if dlon_num < -300:
            dlon_num += 360
        # Calculate difference between lon[ie] and lon[iw], mod 360
        dlon_den = lon[ie] - lon[iw] + 360
    else:
        # General case
        # Find the first index eastwards of lon0
        ie = nonzero(lon > lon0)[0][0]
        # The index before it will be the last index westward of lon0
        iw = ie - 1
        dlon_num = lon0 - lon[iw]
        dlon_den = lon[ie] - lon[iw]
    # Coefficients for interpolation
    coeff1 = dlon_num/dlon_den
    coeff2 = 1 - coeff1

    # Interpolate
    data = coeff1*data_3d[:,:,ie] + coeff2*data_3d[:,:,iw]
    return data


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    file_path1 = raw_input("Path to output oce.mean.nc containing one year of 5-day averages (December will be used): ")
    file_path2 = raw_input("Path to the following oce.mean.nc containing 5-day averages for the next year (January through November will be used): ")
    var_key = raw_input("Temperature (t) or salinity (s)? ")
    if var_key == 't':
        var_name = 'temp'
    elif var_key == 's':
        var_name = 'salt'
    lon0 = float(raw_input("Enter longitude (-180 to 180): "))
    depth_min = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None

    # Build the FESOM mesh ahead of time
    elements = fesom_grid(mesh_path)
    sose_fesom_seasonal (elements, file_path1, file_path2, var_name, lon0, depth_min, save, fig_name)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the input parameters; repeat until the user is finished
                changes = raw_input("Enter a parameter to change: (1) file paths, (2) temperature/salinity, (3) longitude, (4) deepest depth, (5) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New file paths
                        file_path1 = raw_input("Path to one year of 5-day averages for sea ice variables (December will be used): ")
                        file_path2 = raw_input("Path to the following year of 5-day averages for sea ice variables (January through November will be used): ")
                    elif int(changes) == 2:
                        # Switch from temperature to salinity or vice versa
                        if var_name == 'temp':
                            var_name = 'salt'
                        else:
                            var_name = 'temp'
                    elif int(changes) == 3:
                        # New longitude
                        lon0 = float(raw_input("Enter longitude (-180 to 180): "))
                    elif int(changes) == 4:
                        # New depth bound
                        depth_min = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
                    elif int(changes) == 5:
                        # Change from save to display, or vice versa
                        save = not save
            if save:
                # Get file name for figure
                fig_name = raw_input("File name for figure: ")
            # Make the plot
            sose_fesom_seasonal (elements, file_path1, file_path2, var_name, lon0, depth_min, save, fig_name)
        else:
            break

    

