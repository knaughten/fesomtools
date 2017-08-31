from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from fesom_grid import *

def timeseries_massloss_sectors (mesh_path, diag_file, log_file, fig_dir=''):

    # Titles and figure names for each sector
    names = ['Filchner-Ronne Ice Shelf', 'Eastern Weddell Region', 'Amery Ice Shelf', 'Australian Sector', 'Ross Sea', 'Amundsen Sea', 'Bellingshausen Sea', 'Larsen Ice Shelves', 'Total Antarctica']
    fig_names = ['filchner_ronne.png', 'eweddell.png', 'amery.png', 'australian.png', 'ross.png', 'amundsen.png', 'bellingshausen.png', 'larsen.png', 'total_antarctica.png']
    num_sectors = len(fig_names)
    # Seconds per year
    sec_per_year = 365.25*24*3600
    # Density of ice in kg/m^3
    rho_ice = 916
    # Only consider elements south of 30S
    circumpolar = True
    # Don't make second copies of elements that cross 180E
    cross_180 = False
    # Number of days for each output step
    days_per_output = 5

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
        # Set up array for mass loss values for each sector
        old_massloss = empty([num_sectors, start_t])
        # Fill in the first sector
        old_massloss[0,:] = tmp_massloss[:]
        index = 1
        # Loop over the rest of the sectors
        while index < num_sectors:
            t = 0
            for line in f:
                try:
                    old_massloss[index, t] = float(line)
                    t += 1
                except(ValueError):
                    # Reached the header for the next sector
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
    massloss = empty([num_sectors, start_t+num_time])
    if exists(log_file):
        # Fill first start_t timesteps with existing values
        massloss[:,0:start_t] = old_massloss[:,:]
    # Read melt rate and convert from m/s to m/y
    ismr = id.variables['wnet'][:,:]*sec_per_year
    id.close()

    print 'Setting up arrays'
    # Melt rate timeseries at each element
    ismr_elm = zeros([num_time, len(elements)])
    # Area of each element
    area_elm = zeros(len(elements))
    # Flag to indicate which ice shelves the element is part of
    location_flag = zeros([num_sectors, len(elements)])
    for i in range(len(elements)):
        elm = elements[i]
        # Make sure we're actually in an ice shelf cavity
        if elm.cavity:
            # Average ice shelf melt rate timeseries over 3 component nodes
            ismr_elm[:,i] = (ismr[:,elm.nodes[0].id] + ismr[:,elm.nodes[1].id] +ismr[:,elm.nodes[2].id])/3
            # Call area function
            area_elm[i] = elm.area()            
            # Get average lon and lat across 3 Nodes
            lon = mean(elm.lon)
            lat = mean(elm.lat)
            # Figure out which sector this ice shelf element falls into
            if lon >= -85 and lon < -30 and lat < -74:
                # Filchner-Ronne
                location_flag[0,i] = 1           
            elif lon >= -30 and lon < 65:
                # Eastern Weddell region
                location_flag[1,i] = 1       
            elif lon >= 65 and lon < 76:
                # Amery
                location_flag[2,i] = 1       
            elif lon >= 76 and lon < 165 and lat >= -74:
                # Australian sector
                location_flag[3,i] = 1       
            elif (lon >= 155 and lon < 165 and lat < -74) or (lon >= 165) or (lon < -140):
                # Ross Sea
                location_flag[4,i] = 1         
            elif (lon >= -140 and lon < -105) or (lon >= -105 and lon < -98 and lat < -73.1):
                # Amundsen Sea
                location_flag[5,i] = 1          
            elif (lon >= -104 and lon < -98 and lat >= -73.1) or (lon >= -98 and lon < -66 and lat >= -75):
                # Bellingshausen Sea
                location_flag[6,i] = 1
            elif lon >= -66 and lon < -59 and lat >= -74:
                # Larsen Ice Shelves
                location_flag[7,i] = 1
            else:
                print 'No region found for lon=',str(lon),', lat=',str(lat)
                break #return
            # All ice shelf elements are in Total Antarctica
            location_flag[8,i] = 1

    # Calculate conversion factors from mass loss to area-averaged melt rate
    # for each ice shelf
    factors = empty(num_sectors)
    for index in range(num_sectors):
        # Calculate area of the ice shelf
        tmp_area = sum(area_elm*location_flag[index,:])
        factors[index] = 1e12/(rho_ice*tmp_area)
        print 'Area of ' + names[index] + ': ' + str(tmp_area) + ' m^2'

    # Build timeseries
    for t in range(num_time):
        # Loop over sectors
        for index in range(num_sectors):
            # Integrate ice shelf melt rate over area to get volume loss
            volumeloss = sum(ismr_elm[t,:]*area_elm*location_flag[index,:])
            # Convert to mass loss in Gt/y
            massloss[index,start_t+t] = 1e-12*rho_ice*volumeloss

    # Calculate time values
    time = arange(size(massloss,1))*days_per_output/365.

    print 'Plotting'
    for index in range(num_sectors):
        # Set up plot: mass loss and melt rate are directly proportional (with
        # a different constant of proportionality for each ice shelf depending
        # on its area) so plot one line with two y-axes
        fig, ax1 = subplots()
        ax1.plot(time, massloss[index,:], color='black')
        # Title and ticks in blue for left side of the plot
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
        # Title and ticks in red for this side of the plot
        ax2.set_ylabel('Area-Averaged Ice Shelf Melt Rate (m/y)', color='r')
        for t2 in ax2.get_yticklabels():
            t2.set_color('r')
        # Name of the sector for the main title
        title(names[index])
        fig.savefig(fig_dir + fig_names[index])    

    print 'Saving results to log file'
    f = open(log_file, 'w')
    for index in range(num_sectors):
        f.write(names[index] + ' Basal Mass Loss (Gt/y)\n')
        for t in range(size(time)):
            f.write(str(massloss[index, t]) + '\n')
    f.close()   


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    diag_file = raw_input("Path to FESOM forcing.diag.nc output file: ")
    log_file = raw_input("Path to logfile to save values and/or read previously calculated values: ")

    timeseries_massloss_sectors(mesh_path, diag_file, log_file)
    
