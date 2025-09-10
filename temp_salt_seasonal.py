from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from fesom_grid import *
from seasonal_avg import *
from side_patches import *

# Make a 4x2 plot showing seasonally averaged temperature (top) and salinity
# (bottom) interpolated to a given longitude, i.e. latitude vs. depth slices.
# Input:
# elements = FESOM grid elements (from fesom_grid.py)
# file_path1, file_path2 = paths to two FESOM output oce.mean.nc files, each
#                          containing one year of 5-day averages. The script
#                          will use December from file_path1 and January
#                          through November from file_path2.
# lon0 = longitude to interpolate to (-180 to 180)
# depth_min = deepest depth to plot (negative, in metres)
# save = optional boolean indicating to save the figure, rather than display
# fig_name = if save=True, filename for figure
def temp_salt_seasonal (elements, file_path1, file_path2, lon0, depth_min, save=False, fig_name=None):

    # Northern boundary for plot
    lat_max = -30
    # Season names for titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']

    # Bounds on colour scales for temperature and salinity
    var_min = [-2.5, 33.8]
    var_max = [3.5, 34.8]
    var_ticks = [1, 0.2]

    # Choose what to write on the title about longitude
    if lon0 < 0:
        lon_string = str(int(round(-lon0))) + r'$^{\circ}$W'
    else:
        lon_string = str(int(round(lon0))) + r'$^{\circ}$E'

    # Get seasonal averages of temperature and salinity
    temp_data = seasonal_avg(file_path1, file_path2, 'temp')
    salt_data = seasonal_avg(file_path1, file_path2, 'salt')

    # Set colour levels
    lev1 = linspace(var_min[0], var_max[0], num=50)
    lev2 = linspace(var_min[1], var_max[1], num=50)

    # Plot
    fig = figure(figsize=(20,9))
    # Loop over seasons
    for season in range(4):
        print('Calculating zonal slices for ' + season_names[season])
        # Interpolate temperature to lon0 and get plotting patches
        patches, values, lat_min = side_patches(elements, lat_max, lon0, temp_data[season,:])
        ax = fig.add_subplot(2, 4, season+1)
        img = PatchCollection(patches, cmap=jet)
        img.set_array(array(values))
        img.set_edgecolor('face')
        img.set_clim(vmin=var_min[0], vmax=var_max[0])
        ax.add_collection(img)
        xlim([lat_min, lat_max])
        ylim([depth_min, 0])
        title('Temperature (' + season_names[season] + ')', fontsize=24)
        if season == 0:
            # Add depth label on the left
            ylabel('depth (m)', fontsize=18)
        if season == 3:
            # Add a colourbar on the right
            cbaxes = fig.add_axes([0.93, 0.6, 0.015, 0.3])
            cbar1 = colorbar(img, cax=cbaxes, ticks=arange(var_min[0], var_max[0]+var_ticks[0], var_ticks[0]))
            cbar1.ax.tick_params(labelsize=16)
        # Repeat for salinity
        patches, values, lat_min = side_patches(elements, lat_max, lon0, salt_data[season,:])
        ax = fig.add_subplot(2, 4, season+5)
        img = PatchCollection(patches, cmap=jet)
        img.set_array(array(values))
        img.set_edgecolor('face')
        img.set_clim(vmin=var_min[1], vmax=var_max[1])
        ax.add_collection(img)
        xlim([lat_min, lat_max])
        ylim([depth_min, 0])
        title('Salinity (' + season_names[season] + ')', fontsize=24)
        if season == 0:
            ylabel('Depth (m)', fontsize=18)
        xlabel('Latitude', fontsize=18)
        if season == 3:
            cbaxes = fig.add_axes([0.93, 0.1, 0.015, 0.3])
            cbar2 = colorbar(img, cax=cbaxes, ticks=arange(var_min[1], var_max[1]+var_ticks[1], var_ticks[1]))
            cbar2.ax.tick_params(labelsize=16)
    suptitle(lon_string, fontsize=30)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()    


# Command-line interface
if __name__ == "__main__":

    mesh_path = input("Path to FESOM mesh directory: ")
    file_path1 = input("Path to output oce.mean.nc containing one year of 5-day averages (December will be used): ")
    file_path2 = input("Path to the following oce.mean.nc containing 5-day averages for the next year (January through November will be used): ")
    lon0 = float(input("Enter longitude (-180 to 180): "))
    depth_min = -1*float(input("Deepest depth to plot (positive, metres): "))
    action = input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None

    # Build the FESOM mesh ahead of time
    elements = fesom_grid(mesh_path)
    temp_salt_seasonal(elements, file_path1, file_path2, lon0, depth_min, save, fig_name)

    # Repeat until the user wants to exit
    while True:
        repeat = input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the input parameters; repeat until the user is finished
                changes = input("Enter a parameter to change: (1) file paths, (2) longitude, (3) deepest depth, (4) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New file paths
                        file_path1 = input("Path to one year of 5-day averages for sea ice variables (December will be used): ")
                        file_path2 = input("Path to the following year of 5-day averages for sea ice variables (January through November will be used): ")
                    elif int(changes) == 2:
                        # New longitude
                        lon0 = float(input("Enter longitude (-180 to 180): "))
                    elif int(changes) == 3:
                        # New depth bound
                        depth_min = -1*float(input("Deepest depth to plot (positive, metres): "))
                    elif int(changes) == 4:
                        # Change from save to display, or vice versa
                        save = not save
            if save:
                # Get file name for figure
                fig_name = input("File name for figure: ")
            # Make the plot
            temp_salt_seasonal(elements, file_path1, file_path2, lon0, depth_min, save, fig_name)
        else:
            break

    

