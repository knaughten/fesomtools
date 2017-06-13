from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from fesom_grid import *

# Plot timeseries of total basal mass loss and area-averaged ice shelf melt
# rates split up into 3 different depth classes for the ice shelf draft. 
# Input:
# mesh_path = path to FESOM mesh directory
# diag_file = path to FESOM output forcing.diag.nc file
# log_file = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
# fig_dir = optional string containing directory to save figures into. Make
#           sure it ends with a "/". Default is an empty string.
def timeseries_massloss_depth (mesh_path, diag_file, log_file, fig_dir=''):

    # Bounds on depth classes
    draft_min = array([0, 250, 500])
    draft_max = array([250, 500, 3000])
    num_classes = size(draft_min)
    # Labels for legend
    labels = ['<'+str(draft_max[0])+' m']
    for n in range(1, num_classes-1):
        labels.append(str(draft_min[n])+'-'+str(draft_max[n])+' m')
    labels.append('>'+str(draft_min[-1])+' m')

    circumpolar = True   # Only consider elements south of 30S
    cross_180 = False    # Don't make second copies of elements that cross 180E
    days_per_output = 5  # Number of days for each output step
    rho_ice = 916        # Density of ice in kg/m^3
    start_year = 1992    # Assumes 1 repetition of present-day forcing, possibly followed by RCP

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
        # Set up array for mass loss values for each depth class
        old_massloss = empty([num_classes, start_t])
        # Fill in the first depth class
        old_massloss[0,:] = tmp_massloss[:]
        n = 1
        # Loop over the other depth classes
        while n < num_classes:
            t = 0
            for line in f:
                try:
                    old_massloss[n,t] = float(line)
                    t += 1
                except(ValueError):
                    # Reached the header for the next depth class
                    break
            n +=1
    else:
        start_t = 0

    print 'Building grid'
    elements = fesom_grid(mesh_path, circumpolar, cross_180)

    print 'Reading data'
    id = Dataset(diag_file, 'r')
    num_time = id.variables['time'].shape[0]
    # Set up array of mass loss values
    massloss = empty([num_classes, start_t+num_time])
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
    # Flag to indicate which depth class the element is part of
    class_flag = zeros([num_classes, len(elements)])
    # Loop over each element to fill these in
    for i in range(len(elements)):
        elm = elements[i]
        # Make sure we're actually in an ice shelf cavity
        if elm.cavity:
            # Average ice shelf melt rate timeseries over 3 component nodes
            ismr_elm[:,i] = (ismr[:,elm.nodes[0].id] + ismr[:,elm.nodes[1].id] +ismr[:,elm.nodes[2].id])/3
            # Call area function
            area_elm[i] = elm.area()
            # Get ice shelf draft (average depth of surface nodes)
            draft = mean(array([(elm.nodes[0]).depth, (elm.nodes[1]).depth, (elm.nodes[2]).depth]))
            # Loop over depth classes
            found = False
            for n in range(num_classes):
                # Figure out whether or not this element is part of the given depth class
                if draft > draft_min[n] and draft <= draft_max[n]:
                    found = True
                    class_flag[n,i] = 1
            if not found:
                print "Couldn't find a depth class for ice shelf draft " + str(draft)
                return

    # Calculate conversion factors from mass loss to area-averaged melt rate
    # for each depth class
    factors = empty(num_classes)
    for n in range(num_classes):
        # Calculate total ice shelf area in this class
        tmp_area = sum(area_elm*class_flag[n,:])
        print 'Area of ice shelf draft between '+str(draft_min[n])+' and '+str(draft_max[n])+'m: '+str(tmp_area)+' m^2'
        factors[n] = 1e12/(rho_ice*tmp_area)

    # Build timeseries
    for t in range(num_time):
        # Loop over depth classes
        for n in range(num_classes):
            # Integrate ice shelf melt rate over area to get volume loss
            volumeloss = sum(ismr_elm[t,:]*area_elm*class_flag[n,:])
            # Convert to massloss in Gt/y
            massloss[n,start_t+t] = 1e-12*rho_ice*volumeloss

    # Calculate time values
    time = arange(size(massloss,1))*days_per_output/365. + start_year

    print "Plotting"

    # Start with mass loss
    fig, ax = subplots(figsize=(10,6))
    # One line for each depth class
    for n in range(num_classes):
        ax.plot(time, massloss[n,:], label=labels[n], linewidth=2)
    # Configure plot
    title('Basal Mass Loss', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel('Gt/y', fontsize=14)
    xlim([time[0], time[-1]])
    grid(True)
    # Move the plot over to make room for legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    # Make legend
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig(fig_dir + 'massloss_depth.png')

    # Repeat for average melt rate
    fig, ax = subplots(figsize=(10,6))
    for n in range(num_classes):
        ax.plot(time, massloss[n,:]*factors[n], label=labels[n], linewidth=2)
    # Configure plot
    title('Area-Averaged Ice Shelf Melt Rate', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel('m/y', fontsize=14)
    grid(True)
    # Move the plot over to make room for legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    # Make legend
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig(fig_dir + 'ismr_depth.png')

    print 'Saving results to log file'
    f = open(log_file, 'w')
    f.write('Basal Mass Loss for ice shelf drafts <' + str(draft_max[0]) + ' m:\n')
    for t in range(size(time)):
        f.write(str(massloss[0,t]) + '\n')
    for n in range(1, num_classes-1):
        f.write('Basal Mass Loss for ice shelf drafts ' + str(draft_min[n]) + '-' + str(draft_max[n]) + ' m:\n')
        for t in range(size(time)):
            f.write(str(massloss[n,t]) + '\n')
    f.write('Basal Mass Loss for ice shelf drafts >' + str(draft_min[-1]) + 'm:\n')
    for t in range(size(time)):
        f.write(str(massloss[-1,t]) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    diag_file = raw_input("Path to FESOM forcing.diag.nc output file: ")
    log_file = raw_input("Path to logfile to save values and/or read previously calculated values: ")

    timeseries_massloss_depth(mesh_path, diag_file, log_file)
        
        
    
            
        

    
