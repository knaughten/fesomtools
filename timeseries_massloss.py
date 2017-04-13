from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from fesom_grid import *

# Calculate and plot timeseries of basal mass loss and area-averaged ice shelf
# melt rates from major ice shelves and from the entire continent during a 
# FESOM simulation.
# Takes 8 GB memory on raijin for Kaitlin's low_res mesh.
# Input:
# mesh_path = path to FESOM mesh directory
# diag_file = path to output forcing.diag.nc file that contains variable "wnet"
#             (surface freshwater flux), assumed to have 5-day averages
# log_file = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
# fig_dir = optional string containing directory to save figures into. Make
#           sure it ends with a "/". Default is an empty string.
def timeseries_massloss (mesh_path, diag_file, log_file, fig_dir=''):

    # Titles and figure names for each ice shelf
    names = ['All Ice Shelves', 'Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    fig_names = ['total_massloss.png', 'larsen_d.png', 'larsen_c.png', 'wilkins_georgevi_stange.png', 'ronne_filchner.png', 'abbot.png', 'pig.png', 'thwaites.png', 'dotson.png', 'getz.png', 'nickerson.png', 'sulzberger.png', 'mertz.png', 'totten_moscowuni.png', 'shackleton.png', 'west.png', 'amery.png', 'princeharald.png', 'baudouin_borchgrevink.png', 'lazarev.png', 'nivl.png', 'fimbul_jelbart_ekstrom.png', 'brunt_riiserlarsen.png', 'ross.png']
    # Limits on longitude and latitude for each ice shelf
    # These depend on the source geometry, in this case RTopo 1.05
    # Note there is one extra index at the end of each array; this is because
    # the Ross region crosses the line 180W and therefore is split into two
    # We have -181 and 181 not -180 and 180 at this boundary so that
    # elements which cross the boundary are still counted
    lon_min = [-181, -62.67, -65.5, -79.17, -85, -104.17, -102.5, -108.33, -114.5, -135.67, -149.17, -155, 144, 115, 94.17, 80.83, 65, 33.83, 19, 12.9, 9.33, -10.05, -28.33, -181, 158.33]
    lon_max = [181, -59.33, -60, -66.67, -28.33, -88.83, -99.17, -103.33, -111.5, -114.33, -140, -145, 146.62, 123.33, 102.5, 89.17, 75, 37.67, 33.33, 16.17, 12.88, 7.6, -10.33, -146.67, 181]
    lat_min = [-90, -73.03, -69.35, -74.17, -83.5, -73.28, -75.5, -75.5, -75.33, -74.9, -76.42, -78, -67.83, -67.17, -66.67, -67.83, -73.67, -69.83, -71.67, -70.5, -70.75, -71.83, -76.33, -85, -84.5]
    lat_max = [-30, -69.37, -66.13, -69.5, -74.67, -71.67, -74.17, -74.67, -73.67, -73, -75.17, -76.41, -66.67, -66.5, -64.83, -66.17, -68.33, -68.67, -68.33, -69.33, -69.83, -69.33, -71.5, -77.77, -77]
    # Observed mass loss (Rignot 2013) and uncertainty for each ice shelf, in Gt/y
    obs_massloss = [1325, 1.4, 20.7, 135.4, 155.4, 51.8, 101.2, 97.5, 45.2, 144.9, 4.2, 18.2, 7.9, 90.6, 72.6, 27.2, 35.5, -2, 21.6, 6.3, 3.9, 26.8, 9.7, 47.7]
    obs_massloss_error = [235, 14, 67, 40, 45, 19, 8, 7, 4, 14, 2, 3, 3, 8, 15, 10, 23, 3, 18, 2, 2, 14, 16, 34]
    # Observed ice shelf melt rates and uncertainty
    obs_ismr = [0.85, 0.1, 0.4, 3.1, 0.3, 1.7, 16.2, 17.7, 7.8, 4.3, 0.6, 1.5, 1.4, 7.7, 2.8, 1.7, 0.6, -0.4, 0.4, 0.7, 0.5, 0.5, 0.1, 0.1]
    obs_ismr_error = [0.1, 0.6, 1, 0.8, 0.1, 0.6, 1, 1, 0.6, 0.4, 0.3, 0.3, 0.6, 0.7, 0.6, 0.7, 0.4, 0.6, 0.4, 0.2, 0.2, 0.2, 0.2, 0.1]
    # Density of ice in kg/m^3
    rho_ice = 916

    circumpolar = True   # Only consider elements south of 30S
    cross_180 = False    # Don't make second copies of elements that cross 180E
    days_per_output = 5  # Number of days for each output step

    tmp_massloss = []
    # Check if the log file exists
    if exists(log_file):
        print 'Reading previously calculated values'
        f = open(log_file, 'r')
        # Skip the first line (header)
        f.readline()
        for line in f:
            try:
                tmp_massloss.append(float(line))
            except(ValueError):
                # Reached the header for the next variable
                break
        start_t = len(tmp_massloss)
        # Set up array for mass loss values at each ice shelf
        old_massloss = empty([len(names), start_t])
        # Fill in the first timeseries (entire continent)
        old_massloss[0,:] = tmp_massloss[:]
        index = 1
        # Loop over the individual ice shelves
        while index < len(names):
            t = 0
            for line in f:
                try:
                    old_massloss[index, t] = float(line)
                    t += 1
                except(ValueError):
                    # Reached the header for the next ice shelf
                    break
            index +=1
    else:
        start_t = 0

    print 'Building grid'
    elements = fesom_grid(mesh_path, circumpolar, cross_180)

    print 'Reading data'
    id = Dataset(diag_file, 'r')
    num_time = id.variables['time'].shape[0]
    # Set up array of mass loss values
    massloss = empty([len(names), start_t+num_time])
    if exists(log_file):
        # Fill first start_t timesteps with existing values
        massloss[:,0:start_t] = old_massloss[:,:]
    # Read melt rate and convert from m/s to m/y
    ismr = id.variables['wnet'][:,:]*365.25*24*60*60
    id.close()

    print 'Setting up arrays'
    # Melt rate timeseries at each element
    ismr_elm = zeros([num_time, len(elements)])
    # Area of each element
    area_elm = zeros(len(elements))
    # Flag to indicate which ice shelves the element is part of
    location_flag = zeros([len(names), len(elements)])
    # Loop over each element to fill these in
    for i in range(len(elements)):
        elm = elements[i]
        # Make sure we're actually in an ice shelf cavity
        if elm.cavity:
            # Average ice shelf melt rate timeseries over 3 component nodes
            ismr_elm[:,i] = (ismr[:,elm.nodes[0].id] + ismr[:,elm.nodes[1].id] +ismr[:,elm.nodes[2].id])/3
            # Call area function
            area_elm[i] = elm.area()
            # Loop over ice shelves
            for index in range(len(names)):
                # Figure out whether or not this element is part of the given
                # ice shelf
                if all(elm.lon >= lon_min[index]) and all(elm.lon <= lon_max[index]) and all(elm.lat >= lat_min[index]) and all(elm.lat <= lat_max[index]):
                    location_flag[index,i] = 1
                if index == len(names)-1:
                    # Ross region is split into two
                    if all(elm.lon >= lon_min[index+1]) and all(elm.lon <= lon_max[index+1]) and all(elm.lat >= lat_min[index+1]) and all(elm.lat <= lat_max[index+1]):
                        location_flag[index,i] = 1

    # Calculate conversion factors from mass loss to area-averaged melt rate
    # for each ice shelf
    factors = empty(len(names))
    for index in range(len(names)):
        # Calculate area of the ice shelf
        tmp_area = sum(area_elm*location_flag[index,:])
        factors[index] = 1e12/(rho_ice*tmp_area)
        print 'Area of ' + names[index] + ': ' + str(tmp_area) + ' m^2'

    # Build timeseries
    for t in range(num_time):
        # Loop over ice shelves
        for index in range(len(names)):
            # Integrate ice shelf melt rate over area to get volume loss
            volumeloss = sum(ismr_elm[t,:]*area_elm*location_flag[index,:])
            # Convert to mass loss in Gt/y
            massloss[index,start_t+t] = 1e-12*rho_ice*volumeloss

    # Calculate time values
    time = arange(size(massloss,1))*days_per_output/365.

    print 'Plotting'
    for index in range(len(names)):
        # Calculate the bounds on observed mass loss and melt rate
        massloss_low = obs_massloss[index] - obs_massloss_error[index]
        massloss_high = obs_massloss[index] + obs_massloss_error[index]
        ismr_low = obs_ismr[index] - obs_ismr_error[index]
        ismr_high = obs_ismr[index] + obs_ismr_error[index]
        # Set up plot: mass loss and melt rate are directly proportional (with
        # a different constant of proportionality for each ice shelf depending
        # on its area) so plot one line with two y-axes
        fig, ax1 = subplots()
        ax1.plot(time, massloss[index,:], color='black')
        # In blue, add dashed lines for observed mass loss
        ax1.axhline(massloss_low, color='b', linestyle='dashed')
        ax1.axhline(massloss_high, color='b', linestyle='dashed')
        # Make sure y-limits won't cut off observed melt rate
        ymin = amin([ismr_low/factors[index], massloss_low, amin(massloss[index,:])])
        ymax = amax([ismr_high/factors[index], massloss_high, amax(massloss[index,:])])
        # Adjust y-limits to line up with ticks
        ticks = ax1.get_yticks()
        min_tick = ticks[0]
        max_tick = ticks[-1]
        dtick = ticks[1]-ticks[0]
        while min_tick >= ymin:
            min_tick -= dtick
        while max_tick <= ymax:
            max_tick += dtick
        ax1.set_ylim([min_tick, max_tick])
        # Title and ticks in blue for this side of the plot
        ax1.set_ylabel('Basal Mass Loss (Gt/y)', color='b')
        for t1 in ax1.get_yticklabels():
            t1.set_color('b')
        ax1.set_xlabel('Years')
        ax1.grid(True)
        # Twin axis for melt rates
        ax2 = ax1.twinx()
        # Make sure the scales line up
        limits = ax1.get_ylim()
        ax2.set_ylim([limits[0]*factors[index], limits[1]*factors[index]])
        # In red, add dashed lines for observed ice shelf melt rates
        ax2.axhline(ismr_low, color='r', linestyle='dashed')
        ax2.axhline(ismr_high, color='r', linestyle='dashed')
        # Title and ticks in red for this side of the plot
        ax2.set_ylabel('Area-Averaged Ice Shelf Melt Rate (m/y)', color='r')
        for t2 in ax2.get_yticklabels():
            t2.set_color('r')
        # Name of the ice shelf for the main title
        title(names[index])
        fig.savefig(fig_dir + fig_names[index])     

    print 'Saving results to log file'
    f = open(log_file, 'w')
    for index in range(len(names)):
        f.write(names[index] + ' Basal Mass Loss\n')
        for t in range(size(time)):
            f.write(str(massloss[index, t]) + '\n')
    f.close()       


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    diag_file = raw_input("Path to FESOM forcing.diag.nc output file: ")
    log_file = raw_input("Path to logfile to save values and/or read previously calculated values: ")

    timeseries_massloss(mesh_path, diag_file, log_file)
            
                
                
                
        
