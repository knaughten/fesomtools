from netCDF4 import Dataset
from numpy import *
from fesom_grid import *

def timeseries_seaice_extent_faster (mesh_path, output_path, start_year, end_year, log_file):

    circumpolar = True   # Only consider elements south of 30S
    cross_180 = False    # Don't make second copies of elements that cross 180E
    days_per_output = 5  # Number of days for each output step
    expt_name = 'MK44005'

    print 'Building grid'
    elements = fesom_grid(mesh_path, circumpolar, cross_180)

    extent = []
    for year in range(start_year, end_year+1):
        print year
        ice_file = output_path + expt_name + '.' + str(year) + '.ice.mean.nc'
        print 'Reading data'
        id = Dataset(ice_file, 'r')
        num_time = id.variables['time'].shape[0]
        aice = id.variables['area'][:,:]
        id.close()
        print 'Setting up arrays'
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
        print 'Building timeseries'
        for t in range(num_time):
            # Integrate extent and convert to million km^2
            extent.append(sum(flag[t,:]*area_elm)*1e-12)

    print 'Saving results to log file'
    f = open(log_file, 'w')
    f.write('Sea Ice Extent (million km^2):\n')
    for elm in extent:
        f.write(str(elm) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    output_path = raw_input("Path to FESOM output directory: ")
    start_year = int(raw_input("First year to process: "))
    end_year = int(raw_input("Last year to process: "))
    log_file = raw_input("Path to logfile to save values: ")
    timeseries_seaice_extent_faster(mesh_path, output_path, start_year, end_year, log_file)
    
