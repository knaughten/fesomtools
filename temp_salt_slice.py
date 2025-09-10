from numpy import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from unrotate_vector import *
from fesom_grid import *
from fesom_sidegrid import *

# Make a 2x1 plot showing temperature and salinity interpolated to a given
# longitude (i.e. latitude vs. depth slices) at a single time index.
# Input:
# elm2D = FESOM grid elements (from fesom_grid.py)
# file_path = path to FESOM output oce.mean.nc file
# tstep = time index in file_path to plot (1-based)
# lon0 = longitude to interpolate to (-180 to 180)
# depth_min = deepest depth to plot (negative, in metres)
# save = optional boolean indicating to save the figure, rather than display
# fig_name = if save=True, filename for figure
def temp_salt_slice (elm2D, file_path, tstep, lon0, depth_min, save=False, fig_name=None):

    # Set northern boundary and upper (surface) boundary
    lat_max = -30
    depth_max = 0

    # Bounds on colour scales for temperature and salinity
    var_min = [-2, 33.8]
    var_max = [3, 34.8]
    var_tick = [1, 0.2]

    # Read temperature and salinity at each node
    id = Dataset(file_path, 'r')
    temp = id.variables['temp'][tstep-1,:]
    salt = id.variables['salt'][tstep-1,:]
    id.close()

    # Choose what to write on the title about longitude
    if lon0 < 0:
        lon_string = str(-lon0) + 'W'
    else:
        lon_string = str(lon0) + 'E'

    # Set up plots
    fig = figure(figsize=(24,6))
    # Make SideElements with temperature data
    selements_temp = fesom_sidegrid(elm2D, temp, lon0, lat_max)
    # Build an array of quadrilateral patches for the plot, and of data values
    # corresponding to each SideElement
    # Also find the minimum latitude of any SideElement
    patches = []
    values = []
    lat_min = lat_max
    for selm in selements_temp:
        # Make patch
        coord = transpose(vstack((selm.y,selm.z)))
        patches.append(Polygon(coord, True, linewidth=0.))
        # Save data value
        values.append(selm.var)
        # Update minimum latitude if needed
        lat_min = min(lat_min, amin(selm.y))
    # Set southern boundary to be just south of the minimum latitude
    lat_min = lat_min-1
    # Add to plot
    ax1 = fig.add_subplot(1,2,1)
    img1 = PatchCollection(patches)
    img1.set_array(array(values))
    img1.set_edgecolor('face')
    img1.set_clim(vmin=var_min[0], vmax=var_max[0])
    ax1.add_collection(img1)
    xlim(lat_min, lat_max)
    ylim(depth_min, depth_max)
    xlabel('Latitude')
    ylabel('Depth (m)')
    title(r'Temperature ($^{\circ}$C)', fontsize=20)
    cbar1 = colorbar(img1, ticks=arange(var_min[0], var_max[0]+var_tick[0], var_tick[0]), extend='both')
    cbar1.ax.tick_params(labelsize=16)
    # Repeat for salinity
    selements_salt = fesom_sidegrid(elm2D, salt, lon0, lat_max)
    patches = []
    values = []
    for selm in selements_salt:
        coord = transpose(vstack((selm.y,selm.z)))
        patches.append(Polygon(coord, True, linewidth=0.))
        values.append(selm.var)
    ax2 = fig.add_subplot(1,2,2)
    img2 = PatchCollection(patches)
    img2.set_array(array(values))
    img2.set_edgecolor('face')
    img2.set_clim(vmin=var_min[1], vmax=var_max[1])
    ax2.add_collection(img2)
    xlim(lat_min, lat_max)
    ylim(depth_min, depth_max)
    xlabel('Latitude')
    ylabel('Depth (m)')
    title('Salinity (psu)', fontsize=20)
    cbar2 = colorbar(img2, ticks=arange(var_min[1], var_max[1]+var_tick[1], var_tick[1]), extend='both')
    cbar2.ax.tick_params(labelsize=16)
    suptitle('Time index ' + str(tstep) + ', ' + lon_string, fontsize=24)
    subplots_adjust(wspace=0.025)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    mesh_path = input("Path to FESOM mesh directory: ")
    file_path = input("Path to FESOM oce.mean.nc file: ")
    tstep = int(input("Time index to plot (starting at 1): "))
    lon0 = float(input("Longitude in degrees (-180 to 180): "))
    depth_min = -1*float(input("Deepest depth to plot (positive, metres): "))
    action = input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None    
    elm2D = fesom_grid(mesh_path)
    temp_salt_slice(elm2D, file_path, tstep, lon0, depth_min, save, fig_name)

    # Repeat until the user is finished
    while True:
        repeat = input("Make another plot (y/n)? ")
        if repeat == 'y':
            update_mesh = False
            while True:                
                changes = input("Enter a parameter to change: (1) mesh path, (2) file path, (3) timestep, (4) longitude, (5) deepest depth, (6) save/display; or enter to continue: ")
                if len(changes) == 0:
                    break
                else:
                    if int(changes) == 1:
                        update_mesh = True
                        mesh_path = input("Path to FESOM mesh directory: ")
                        file_path = input("Path to FESOM oce.mean.nc file: ")
                    elif int(changes) == 2:
                        file_path = input("Path to FESOM oce.mean.nc file: ")
                    elif int(changes) == 3:
                        tstep = int(input("Time index to plot (starting at 1): "))
                    elif int(changes) == 4:
                        lon0 = float(input("Longitude in degrees (-180 to 180): "))
                    elif int(changes) == 5:
                        depth_min = -1*float(input("Deepest depth to plot (positive, metres): "))
                    elif int(changes) == 6:
                        save = not save
            if save:
                fig_name = input("File name for figure: ")
            if update_mesh:
                elm2D = fesom_grid(mesh_path)
            temp_salt_slice(elm2D, file_path, tstep, lon0, depth_min, save, fig_name)
        else:
            break
