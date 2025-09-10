from netCDF4 import Dataset
from numpy import *
from fesom_grid import *

# Calculate timeseries of basal mass loss for the Wilkins, Stange, and
# George VI ice shelves separately. Output to a log file.
# Input:
# mesh_path = path to FESOM mesh directory
# directory = path to output directory for FESOM simulation
# start_year, end_year = integers containing first and last years to process
# log_file = path to desired output log file
def timeseries_massloss_bellingshausen (mesh_path, directory, start_year, end_year, log_file):

    # Name of each ice shelf
    names = ['Wilkins Ice Shelf', 'Stange Ice Shelf', 'George VI Ice Shelf']
    # Limits on longitude and latitude for each ice shelf
    # These depend on the source geometry, in this case RTopo 1.05
    # Note there is one extra index at the end of each array; this is because
    # the George VI region is awkwardly shaped and split into 2 (and hence why
    # all 3 ice shelves are lumped together in timeseries_massloss.py)
    lon_min = [-75, -79, -74.5, -69.5]
    lon_max = [-69, -74.5, -67, -66]
    lat_min = [-71.5, -73.8, -73.8, -72.6]
    lat_max = [-69, -72.6, -72.6, -70]
    # Density of ice in kg/m^3
    rho_ice = 916
    # Beginning of each file name
    expt_name = 'MK44005'

    circumpolar = True   # Only consider elements south of 30S
    cross_180 = False    # Don't make second copies of elements that cross 180E
    peryear = 365//5      # Number of output steps per year (5-day averages)

    print('Building grid')
    elements = fesom_grid(mesh_path, circumpolar, cross_180)

    print('Setting up arrays')
    # Area of each ice shelf element
    area_elm = zeros(len(elements))
    # Flag to indicate which ice shelves the element is part of
    location_flag = zeros([len(names), len(elements)])
    # Loop over each element to fill these in
    for i in range(len(elements)):
        elm = elements[i]
        # Make sure we're actually in an ice shelf cavity
        if elm.cavity:
            # Call area function
            area_elm[i] = elm.area()
            # Loop over ice shelves
            for shelf in range(len(names)):
                # Figure out whether or not this element is part of the given
                # ice shelf
                if all(elm.lon >= lon_min[shelf]) and all(elm.lon <= lon_max[shelf]) and all(elm.lat >= lat_min[shelf]) and all(elm.lat <= lat_max[shelf]):
                    location_flag[shelf,i] = 1
                if shelf == len(names)-1:
                    # George VI region is split into two
                    if all(elm.lon >= lon_min[shelf+1]) and all(elm.lon <= lon_max[shelf+1]) and all(elm.lat >= lat_min[shelf+1]) and all(elm.lat <= lat_max[shelf+1]):
                        location_flag[shelf,i] = 1
    # Set up array of mass loss values for each ice shelf
    num_time = (end_year-start_year+1)*peryear
    massloss = empty([len(names), num_time])

    t_posn = 0
    for year in range(start_year, end_year+1):
        file = directory + expt_name + '.' + str(year) + '.forcing.diag.nc'
        print('Processing ' + file)
        # Read melt rate and convert from m/s to m/y
        id = Dataset(file, 'r')        
        ismr = id.variables['wnet'][:,:]*365.25*24*60*60
        id.close()
        # Calculate melt rate timeseries at each element
        ismr_elm = zeros([peryear, len(elements)])
        for i in range(len(elements)):
            elm = elements[i]
            # Make sure we're actually in an ice shelf cavity
            if elm.cavity:
                # Average ice shelf melt rate timeseries over 3 component nodes
                ismr_elm[:,i] = (ismr[:,elm.nodes[0].id] + ismr[:,elm.nodes[1].id] +ismr[:,elm.nodes[2].id])/3
        # Build timeseries
        for t in range(peryear):
            for shelf in range(len(names)):
                # Integrate ice shelf melt rate over area to get volume loss
                volumeloss = sum(ismr_elm[t,:]*area_elm*location_flag[shelf,:])
                # Convert to mass loss in Gt/y
                massloss[shelf, t_posn+t] = 1e-12*rho_ice*volumeloss
        t_posn += peryear

    print('Saving results to log file')
    f = open(log_file, 'w')
    for index in range(len(names)):
        f.write(names[index] + ' Basal Mass Loss\n')
        for t in range(num_time):
            f.write(str(massloss[index, t]) + '\n')
    f.close()   


# Command-line interface
if __name__ == "__main__":

    mesh_path = input("Path to FESOM mesh directory: ")
    directory = input("Path to FESOM output directory: ")
    start_year = int(input("First year to process: "))
    end_year = int(input("Last year to process: "))
    log_file = input("Path to desired logfile: ")
    timeseries_massloss_bellingshausen(mesh_path, directory, start_year, end_year, log_file)

            
