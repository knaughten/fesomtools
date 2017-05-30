from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.pyplot import *
from matplotlib.colors import LinearSegmentedColormap
from patches import *

# Make a circumpolar plot of the percentage change in mass loss for each major
# ice shelf over the given RCP simulation (2091-2100 average minus 2006-2015).
# Input:
# mesh_path = path to FESOM mesh directory
# log_path = path to logfile from timeseries_massloss.py, containing 5-day
#            averages from 1992 to 2100 inclusive
# save = optional boolean indicating to save the figure rather than display
# fig_name = if save=True, filename for figure
def massloss_percent_change (mesh_path, log_path, save=False, fig_name=None):

    # Limits on longitude and latitude for each ice shelf
    # These depend on the source geometry, in this case RTopo 1.05
    # Note there is one extra index at the end of each array; this is because
    # the Ross region crosses the line 180W and therefore is split into two
    # We have -181 and 181 not -180 and 180 at this boundary so that
    # elements which cross the boundary are still counted
    lon_min = [-62.67, -65.5, -79.17, -85, -104.17, -102.5, -108.33, -114.5, -135.67, -149.17, -155, 144, 115, 94.17, 80.83, 65, 33.83, 19, 12.9, 9.33, -10.05, -28.33, -181, 158.33]
    lon_max = [-59.33, -60, -66.67, -28.33, -88.83, -99.17, -103.33, -111.5, -114.33, -140, -145, 146.62, 123.33, 102.5, 89.17, 75, 37.67, 33.33, 16.17, 12.88, 7.6, -10.33, -146.67, 181]
    lat_min = [-73.03, -69.35, -74.17, -83.5, -73.28, -75.5, -75.5, -75.33, -74.9, -76.42, -78, -67.83, -67.17, -66.67, -67.83, -73.67, -69.83, -71.67, -70.5, -70.75, -71.83, -76.33, -85, -84.5]
    lat_max = [-69.37, -66.13, -69.5, -74.67, -71.67, -74.17, -74.67, -73.67, -73, -75.17, -76.41, -66.67, -66.5, -64.83, -66.17, -68.33, -68.67, -68.33, -69.33, -69.83, -69.33, -71.5, -77.77, -77]
    num_shelves = len(lat_max)-1

    # Plotting parameters
    max_lat_plot = -63+90
    mask_cavities = True
    circumpolar = True
    # Assume timeseries are 5-day averages
    peryear = 365/5
    # First 14 years don't count as RCP
    skipyears=14

    # Build FESOM mesh
    # Get separate patches for the open ocean and minor ice shelf elements
    # so we can mask them out
    elements, mask_patches = make_patches(mesh_path, circumpolar, mask_cavities, only_major=True)
    patches = iceshelf_mask(elements, only_major=True)

    # Read log file
    f = open(log_path, 'r')
    # Skip the first line (header)
    f.readline()
    # Count the number of time indices for the first variable (total mass loss
    # for all ice shelves, which we don't care about)
    num_time = 0
    for line in f:
        try:
            tmp = float(line)
            num_time += 1
        except(ValueError):
            # Reached header for first individual ice shelf
            break
    # Set up array for mass loss percent change at each ice shelf
    massloss_change = empty(num_shelves)
    index = 0
    # Loop over ice shelves
    while index < num_shelves:
        massloss_tmp = []
        for line in f:
            try:
                massloss_tmp.append(float(line))
            except(ValueError):
                # Reached the header for the next ice shelf
                break
        massloss_start = mean(array(massloss_tmp[skipyears*peryear:(skipyears+10)*peryear]))
        massloss_end = mean(array(massloss_tmp[-10*peryear:]))
        massloss_change[index] = (massloss_end - massloss_start)/massloss_start*100
        index += 1

    # Set up colour map
    max_val = amax(massloss_change)
    cmap_vals = array([0, max_val/4, max_val/2, max_val])
    cmap_colours = [(1,1,1),(1,0.9,0.4),(0.99,0.59,0.18),(0.5,0.0,0.08)]
    cmap_vals_norm = cmap_vals/max_val
    cmap_list = []
    for i in range(size(cmap_vals)):
        cmap_list.append((cmap_vals_norm[i], cmap_colours[i]))
    my_cmap = LinearSegmentedColormap.from_list('my_cmap',cmap_list)

    # Save these values to the right Elements
    values = []
    for elm in elements:
        # Make sure we're actually in an ice shelf cavity
        if elm.cavity:
            keep = False
            # Loop over ice shelves
            for index in range(num_shelves):
                # Figure out whether or not this element is part of the given
                # ice shelf
                if all(elm.lon >= lon_min[index]) and all(elm.lon <= lon_max[index]) and all(elm.lat >= lat_min[index]) and all(elm.lat <= lat_max[index]):
                    keep = True
                    value_tmp = massloss_change[index]
                if index == num_shelves-1:
                    # Ross region is split into two
                    if all(elm.lon >= lon_min[index+1]) and all(elm.lon <= lon_max[index+1]) and all(elm.lat >= lat_min[index+1]) and all(elm.lat <= lat_max[index+1]):
                        keep = True
                        value_tmp = massloss_change[index]
            if keep:
                values.append(value_tmp)                

    # Set up a grey square covering the domain, anything that isn't covered
    # up later is land
    x_reg, y_reg = meshgrid(linspace(-max_lat_plot, max_lat_plot, num=100), linspace(-max_lat_plot, max_lat_plot, num=100))
    land_square = zeros(shape(x_reg))

    # Plot
    fig = figure(figsize=(16,12))
    ax = fig.add_subplot(1,1,1,aspect='equal')
    # Start with grey square background for land
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches, cmap=my_cmap)
    img.set_array(array(values))
    img.set_edgecolor('face')
    img.set_clim(0, max_val)
    ax.add_collection(img)
    # Mask out the open ocean in white
    overlay = PatchCollection(mask_patches, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)

    # Contour ice shelf front
    contour_lines = []
    for elm in elements:
        # Select elements where exactly 2 of the 3 nodes are in a cavity
        if count_nonzero(elm.cavity_nodes) == 2:
            # Save the coastal flags and x- and y- coordinates of these 2
            coast_tmp = []
            x_tmp = []
            y_tmp = []
            for i in range(3):
                if elm.cavity_nodes[i]:
                    coast_tmp.append(elm.coast_nodes[i])
                    x_tmp.append(elm.x[i])
                    y_tmp.append(elm.y[i])
            # Select elements where at most 1 of these 2 nodes are coastal
            if count_nonzero(coast_tmp) < 2:
                # Draw a line between the 2 nodes
                contour_lines.append([(x_tmp[0], y_tmp[0]), (x_tmp[1], y_tmp[1])])
    # Add all the lines to the plot
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax.add_collection(contours)

    # Configure plot
    xlim([-max_lat_plot, max_lat_plot])
    ylim([-max_lat_plot, max_lat_plot])
    axis('off')
    title('Percent Change in Ice Shelf Mass Loss\n2091-2100 vs 2006-2015', fontsize=30)
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=20)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    log_path = raw_input("Path to mass loss logfile for RCP: ")
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    massloss_percent_change(mesh_path, log_path, save, fig_name)
        
                
                

    
