from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from fesom_grid import *

# Calculate and plot timeseries of total sea ice area and volume during a
# FESOM simulation.
# Input:
# mesh_path = path to FESOM mesh directory
# ice_file = path to output ice.mean.nc, assumed to have 5-day averages
# log_file = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
def timeseries_seaice (mesh_path, ice_file, log_file):

    circumpolar = True   # Only consider elements south of 30S
    cross_180 = False    # Don't make second copies of elements that cross 180E
    days_per_output = 5  # Number of days for each output step

    total_area = []
    total_volume = []
    # Check if the log file exists
    if exists(log_file):
        print 'Reading previously calculated values'
        f = open(log_file, 'r')
        # Skip the first line (header)
        f.readline()
        for line in f:
            try:
                total_area.append(float(line))
            except(ValueError):
                # Reached the header for the next variable
                break
        for line in f:
            total_volume.append(float(line))
        f.close()

    print 'Building grid'
    elements = fesom_grid(mesh_path, circumpolar, cross_180)

    print 'Reading data'
    id = Dataset(ice_file, 'r')
    num_time = id.variables['time'].shape[0]
    aice = id.variables['area'][:,:]
    hice = id.variables['hice'][:,:]
    id.close()

    print 'Setting up arrays'
    # Sea ice concentration at each element
    aice_elm = zeros([num_time, len(elements)])
    # Sea ice height at each element
    hice_elm = zeros([num_time, len(elements)])
    # Area of each element
    area_elm = zeros(len(elements))
    # Loop over elements to fill these in
    for i in range(len(elements)):
        elm = elements[i]
        # Average aice and hi over 3 component nodes
        aice_elm[:,i] = (aice[:,elm.nodes[0].id] + aice[:,elm.nodes[1].id] + aice[:,elm.nodes[2].id])/3
        hice_elm[:,i] = (hice[:,elm.nodes[0].id] + hice[:,elm.nodes[1].id] + hice[:,elm.nodes[2].id])/3
        # Call area function
        area_elm[i] = elm.area()

    # Build timeseries
    for t in range(num_time):
        # Integrate area and convert to million km^2
        total_area.append(sum(aice_elm[t,:]*area_elm)*1e-12)
        # Integrate volume and convert to million km^3
        total_volume.append(sum(aice_elm[t,:]*hice_elm[t,:]*area_elm)*1e-12)

    # Calculate time values
    time = arange(len(total_area))*days_per_output/365.

    print 'Plotting total sea ice area'
    clf()
    plot(time, total_area)
    xlabel('Years')
    ylabel(r'Total Sea Ice Area (million km$^2$)')
    grid(True)
    savefig('seaice_area.png')

    print 'Plotting total sea ice volume'
    clf()
    plot(time, total_volume)
    xlabel('Years')
    ylabel(r'Total Sea Ice Volume (million km$^2$)')
    grid(True)
    savefig('seaice_volume.png')

    print 'Saving results to log file'
    f = open(log_file, 'w')
    f.write('Total Sea Ice Area (million km^2):\n')
    for elm in total_area:
        f.write(str(elm) + '\n')
    f.write('Total Sea Ice Volume (million km^3):\n')
    for elm in total_volume:
        f.write(str(elm) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    ice_file = raw_input("Path to FESOM ice.mean.nc output file: ")
    log_file = raw_input("Path to logfile to save values and/or read previously calculated values: ")

    timeseries_seaice(mesh_path, ice_file, log_file)
