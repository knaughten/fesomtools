from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from fesom_grid import *

# Calculate and plot timeseries of sea ice extent (area of ice with
# concentration >= 15%) during a FESOM simulation.
# mesh_path = path to FESOM mesh directory
# ice_file = path to output ice.mean.nc, assumed to have 5-day averages
# log_file = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
def timeseries_seaice_extent (mesh_path, ice_file, log_file, fig_dir=''):

    circumpolar = True   # Only consider elements south of 30S
    cross_180 = False    # Don't make second copies of elements that cross 180E
    days_per_output = 5  # Number of days for each output step

    extent = []
    # Check if the log file exists
    if exists(log_file):
        print('Reading previously calculated values')
        f = open(log_file, 'r')
        # Skip the first line (header)
        f.readline()
        for line in f:
            extent.append(float(line))
        f.close()

    print('Building grid')
    elements = fesom_grid(mesh_path, circumpolar, cross_180)

    print('Reading data')
    id = Dataset(ice_file, 'r')
    num_time = id.variables['time'].shape[0]
    aice = id.variables['area'][:,:]
    id.close()

    print('Setting up arrays')
    # Sea ice concentration at each element
    aice_elm = zeros([num_time, len(elements)])
    # Area of each element
    area_elm = zeros(len(elements))
    # Loop over elements to fill these in
    for i in range(len(elements)):
        elm = elements[i]
        # Average aice over 3 component nodes
        aice_elm[:,i] = (aice[:,elm.nodes[0].id] + aice[:,elm.nodes[1].id] + aice[:,elm.nodes[2].id])/3
        # Call area function
        area_elm[i] = elm.area()
    # Select elements with concentration >= 15%
    flag = aice_elm >= 0.15

    print('Building timeseries')
    for t in range(num_time):
        # Integrate extent and convert to million km^2
        extent.append(sum(flag[t,:]*area_elm)*1e-12)

    # Calculate time values
    time = arange(len(extent))*days_per_output/365.

    print('Plotting')
    clf()
    plot(time, extent)
    xlabel('Years')
    ylabel(r'Sea Ice Extent (million km$^2$)')
    grid(True)
    savefig(fig_dir+'seaice_extent.png')

    print('Saving results to log file')
    f = open(log_file, 'w')
    f.write('Sea Ice Extent (million km^2):\n')
    for elm in extent:
        f.write(str(elm) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    mesh_path = input("Path to FESOM mesh directory: ")
    ice_file = input("Path to FESOM ice.mean.nc output file: ")
    log_file = input("Path to logfile to save values and/or read previously calculated values: ")

    timeseries_seaice_extent(mesh_path, ice_file, log_file)
    
