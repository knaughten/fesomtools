from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.pyplot import *
from patches import *

# Make a map of unexplained percent error in annually averaged simulated melt
# rate from each ice shelf that is over 5,000 km^2 in Rignot et al., 2013.
# Input:
# mesh_path = path to FESOM mesh directory
# log_path = path to log file created by timeseries_massloss.py
# save = optional boolean to save the figure to a file, rather than displaying
#        it on the screen
# fig_name = if save=True, path to the desired filename for the figure
def ismr_map (mesh_path, log_path, save=False, fig_name=None):

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
    # Area of each ice shelf in m^2 (printed to screen during
    # timeseries_massloss.py, update if the mesh changes)
    area = [8913349811.85, 46584257183.0, 42808123452.1, 427584932349.0, 24722762717.7, 3781726582.37, 2621932896.61, 3962399239.18, 28888280291.9, 4135254416.01, 11700133115.8, 3582900716.02, 3789935863.14, 25475880917.6, 11521918479.5, 63540422738.5, 2734104115.99, 40342528417.8, 6779671751.73, 5666704758.2, 52230325985.9, 70048756435.4, 475369140579.0]
    # Observed melt rate (Rignot 2013) and uncertainty for each ice shelf, in Gt/y
    obs_ismr = [0.1, 0.4, 3.1, 0.3, 1.7, 16.2, 17.7, 7.8, 4.3, 0.6, 1.5, 1.4, 7.7, 2.8, 1.7, 0.6, -0.4, 0.4, 0.7, 0.5, 0.5, 0.1, 0.1]
    obs_ismr_error = [0.6, 1, 0.8, 0.1, 0.6, 1, 1, 0.6, 0.4, 0.3, 0.3, 0.6, 0.7, 0.6, 0.7, 0.4, 0.6, 0.4, 0.2, 0.2, 0.2, 0.2, 0.1]
    # Density of ice in kg/m^3
    rho_ice = 916

    # Plotting parameters
    max_lat_plot = -63+90
    mask_cavities = True
    circumpolar = True
    # Assume timeseries are 5-day averages
    days_per_output = 5
    output_per_year = 365.0/days_per_output

    # Build FESOM mesh
    # Get separate patches for the open ocean elements so we can mask them out
    elements, mask_patches = make_patches(mesh_path, circumpolar, mask_cavities)
    patches = iceshelf_mask(elements)

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
    # Set up array for melt rate values at each ice shelf
    ismr_ts = empty([len(obs_ismr), num_time])
    index = 0
    # Loop over ice shelves
    while index < len(obs_ismr):
        t = 0
        for line in f:
            try:
                # Convert from mass loss to area-averaged melt rate
                ismr_ts[index,t] = float(line)*1e12/(rho_ice*area[index])
                t += 1
            except(ValueError):
                # Reached the header for the next ice shelf
                break
        index += 1

    # Find unexplained error in annual average
    ismr = empty(len(obs_ismr))
    error_vals = empty(len(obs_ismr))
    for index in range(len(obs_ismr)):
        # Calculate annual average
        ismr[index] = mean(ismr_ts[index, -output_per_year:])
        # Calculate range of observations
        ismr_low = obs_ismr[index] - obs_ismr_error[index]
        ismr_high = obs_ismr[index] + obs_ismr_error[index]
        # Calculate unexplained percent error in melt rate
        if ismr[index] < ismr_low:
            # Simulated melt rate too low
            error_vals[index] = (ismr[index] - ismr_low)/ismr_low*100
        elif ismr[index] > ismr_high:
            # Simulated melt rate too high
            error_vals[index] = (ismr[index] - ismr_high)/ismr_high*100
        else:
            # Simulated melt rate within observational error estimates
            error_vals[index] = 0

    # Scale for plotting
    max_val = 100 #amax(abs(error_vals))

    # Build a field of ice shelf melt rate unexplained percent error
    values = []
    # Keep track of which elements are in ice shelves we care about
    location_flag = zeros(len(elements))
    for i in range(len(elements)):
        elm = elements[i]
        # Make sure we're actually in an ice shelf cavity
        if elm.cavity:
            error_tmp = 0.0
            # Loop over ice shelves
            for index in range(len(obs_ismr)):
                # Figure out whether or not this element is part of the given
                # ice shelf
                if all(elm.lon >= lon_min[index]) and all(elm.lon <= lon_max[index]) and all(elm.lat >= lat_min[index]) and all(elm.lat <= lat_max[index]):
                    location_flag[i] = 1
                    error_tmp = error_vals[index]
                if index == len(obs_ismr)-1:
                    # Ross region is split in two
                    if all(elm.lon >= lon_min[index+1]) and all(elm.lon <= lon_max[index+1]) and all(elm.lat >= lat_min[index+1]) and all(elm.lat <= lat_max[index+1]):
                        location_flag[i] = 1
                        error_tmp = error_vals[index]
            values.append(error_tmp)

    # Plot
    fig = figure(figsize=(16,12))
    ax = fig.add_subplot(1,1,1,aspect='equal')
    img = PatchCollection(patches, cmap='RdBu_r')
    img.set_array(array(values))
    img.set_edgecolor('face')
    img.set_clim(-max_val, max_val)
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
    title('Bias in Ice Shelf Melt Rate (%)', fontsize=30)
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
    log_path = raw_input("Path to mass loss logfile: ")
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    ismr_map(mesh_path, log_path, save, fig_name)
        
                
                

    
