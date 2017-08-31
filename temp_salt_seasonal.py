from netCDF4 import Dataset
from numpy import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from fesom_grid import *
from fesom_sidegrid import *
from seasonal_avg import *

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
        print 'Calculating zonal slices for ' + season_names[season]
        # Interpolate temperature to lon0
        patches, values, lat_min = interp_lon_fesom(elements, lat_max, lon0, temp_data[season,:])
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
        patches, values, lat_min = interp_lon_fesom(elements, lat_max, lon0, salt_data[season,:])
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


# Linearly interpolate FESOM data to the specified longitude.
# Input:
# elements = array of 2D Elements created by fesom_grid.py
# lat_max = maximum latitude to consider
# lon0 = longitude to interpolate to, from -180 to 180
# data = array of FESOM data on original mesh
def interp_lon_fesom (elements, lat_max, lon0, data): 

    snode_pairs = []
    for elm in elements:
        # Don't consider elements outside the Southern Ocean
        if any(elm.y <= lat_max):
            # Select elements which intersect lon0
            if any(elm.x <= lon0) and any(elm.x >= lon0):
                # Special case where nodes (corners) of the element are
                # exactly at longitude lon0
                if any(elm.x == lon0):
                    # If exactly one of the corners is at lon0, ignore it;
                    # this element only touches lon0 at one point
                    # If two of the corners are at lon0, an entire side of
                    # the element lies along the line lon0
                    if count_nonzero(elm.x == lon0) == 2:
                        # Select these two Nodes
                        index = nonzero(elm.x == lon0)
                        nodes = elm.nodes[index]
                        node1 = nodes[0]
                        node2 = nodes[1]
                        # Convert to SideNodes and add them to snode_pairs
                        coincide_snode(node1, node2, data, snode_pairs)
                    # Impossible for all three corners to be at lon0
                else:
                    # Regular case
                    snodes_curr = []
                    # Find the two sides of the triangular element which
                    # intersect longitude lon0
                    # For each such side, interpolate a SideNode between the
                    # two endpoint Nodes.
                    if any(array([elm.x[0], elm.x[1]]) < lon0) and any(array([elm.x[0], elm.x[1]]) > lon0):
                        snodes_curr.append(interp_snode(elm.nodes[0], elm.nodes[1], lon0, data))
                    if any(array([elm.x[1], elm.x[2]]) < lon0) and any(array([elm.x[1], elm.x[2]]) > lon0):
                        snodes_curr.append(interp_snode(elm.nodes[1], elm.nodes[2], lon0, data))
                    if any(array([elm.x[0], elm.x[2]]) < lon0) and any(array([elm.x[0], elm.x[2]]) > lon0):
                        snodes_curr.append(interp_snode(elm.nodes[0], elm.nodes[2], lon0, data))
                    # Add the two resulting SideNodes to snode_pairs
                    snode_pairs.append(SideNodePair(snodes_curr[0], snodes_curr[1]))
    selements = []
    # Build the quadrilateral SideElements
    for pair in snode_pairs:
        # Start at the surface
        snode1_top = pair.south
        snode2_top = pair.north
        while True:
            # Select the SideNodes directly below
            snode1_bottom = snode1_top.below
            snode2_bottom = snode2_top.below
            if snode1_bottom is None or snode2_bottom is None:
                # Reached the bottom, so stop
                break
            # Make a SideElement from these four SideNodes
            # The order they are passed to the SideElement initialisation
            # function is important: must trace continuously around the
            # border of the SideElement, i.e. not jump between diagonal
            # corners.
            selements.append(SideElement(snode1_top, snode2_top, snode2_bottom, snode1_bottom))
            # Get ready for the next SideElement below
            snode1_top = snode1_bottom
            snode2_top = snode2_bottom
    # Build an array of quadrilateral patches for the plot, and of data
    # values corresponding to each SideElement
    patches = []
    values = []
    lat_min = lat_max
    for selm in selements:
        # Make patch
        coord = transpose(vstack((selm.y,selm.z)))
        patches.append(Polygon(coord, True, linewidth=0.))
        # Save data value
        values.append(selm.var)
        lat_min = min(lat_min, amin(selm.y))
    # Show a little bit of the land mask
    lat_min = lat_min-0.5

    return patches, values, lat_min


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    file_path1 = raw_input("Path to output oce.mean.nc containing one year of 5-day averages (December will be used): ")
    file_path2 = raw_input("Path to the following oce.mean.nc containing 5-day averages for the next year (January through November will be used): ")
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
    temp_salt_seasonal(elements, file_path1, file_path2, lon0, depth_min, save, fig_name)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the input parameters; repeat until the user is finished
                changes = raw_input("Enter a parameter to change: (1) file paths, (2) longitude, (3) deepest depth, (4) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New file paths
                        file_path1 = raw_input("Path to one year of 5-day averages for sea ice variables (December will be used): ")
                        file_path2 = raw_input("Path to the following year of 5-day averages for sea ice variables (January through November will be used): ")
                    elif int(changes) == 2:
                        # New longitude
                        lon0 = float(raw_input("Enter longitude (-180 to 180): "))
                    elif int(changes) == 3:
                        # New depth bound
                        depth_min = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
                    elif int(changes) == 4:
                        # Change from save to display, or vice versa
                        save = not save
            if save:
                # Get file name for figure
                fig_name = raw_input("File name for figure: ")
            # Make the plot
            temp_salt_seasonal(elements, file_path1, file_path2, lon0, depth_min, save, fig_name)
        else:
            break

    

