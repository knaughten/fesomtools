from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from fesom_grid import *

# Plot timeseries of the annually-averaged, volume-averaged temperature and
# salinity between the given latitude and depth bounds.
# Input:
# mesh_path = path to FESOM mesh directory
# output_path = path to FESOM experiment directory containing all oce.mean.nc
#               files (one for each year)
# lat_bounds = array of size 2 containing latitudes to average between, in the
#              form [south_bound, north_bound], between -90 and 90
# depth_bounds = array of size 2 containing depths to average between, in the
#                form [shallow_bound, deep_bound], positive, in metres
# start_year, end_year = integers containing range of years to process
# log_file = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
# fig_name = filename to save figure
def timeseries_drift (mesh_path, output_path, lat_bounds, depth_bounds, start_year, end_year, log_file, fig_name):

    if lat_bounds[0] > -30:
        circumpolar = False
    else:
        circumpolar = True
    cross_180 = False
    # Naming conventions for FESOM output files
    file_head = output_path + 'MK44005.'
    file_tail = '.oce.mean.nc'

    # Figure out what to say in the title about latitude
    if lat_bounds[0] < 0:
        lat_min_string = str(-lat_bounds[0]) + r'$^{\circ}$S'
    else:
        lat_min_string = str(lat_bounds[0]) + r'$^{\circ}$N'
    if lat_bounds[1] < 0:
        lat_max_string = str(-lat_bounds[1]) + r'$^{\circ}$S'
    else:
        lat_max_string = str(lat_bounds[1]) + r'$^{\circ}$N'

    temp_avg = []
    salt_avg = []
    # Check if the log file exists
    if exists(log_file):
        print('Reading previously calculated values')
        f = open(log_file, 'r')
        # Skip the first line (header)
        f.readline()
        for line in f:
            try:
                temp_avg.append(float(line))
            except(ValueError):
                # Reached the header for the next variable
                break
        for line in f:
            salt_avg.append(float(line))
        f.close()

    print('Building grid')
    elements = fesom_grid(mesh_path, circumpolar, cross_180)

    # Loop over years
    for year in range(start_year, end_year+1):
        print('Processing year ' + str(year))
        # Initialise integrals
        temp_int = 0.0
        salt_int = 0.0
        volume_int = 0.0
        # Read temperature and salinity for this year, annually average
        id = Dataset(file_head + str(year) + file_tail, 'r')
        temp = mean(id.variables['temp'][:,:], axis=0)
        salt = mean(id.variables['salt'][:,:], axis=0)
        id.close()
        # Loop over elements
        for elm in elements:
            # Check if between latitude bounds
            if all(elm.lat >= lat_bounds[0]) and all(elm.lat <= lat_bounds[1]):
                # Get area of 2D element
                area = elm.area()
                nodes = [elm.nodes[0], elm.nodes[1], elm.nodes[2]]
                # Loop downward
                while True:
                    if nodes[0].below is None or nodes[1].below is None or nodes[2].below is None:
                        # Reached the bottom
                        break
                    # Get depth for all 6 nodes of triangular prism
                    node_depths = array([nodes[0].depth, nodes[1].depth, nodes[2].depth, nodes[0].below.depth, nodes[1].below.depth, nodes[2].below.depth])
                    if amax(node_depths) > depth_bounds[1]:
                        # At least 1 node is too deep; exit the loop
                        break
                    # Check whether all nodes are deep enough yet
                    if amin(node_depths) >= depth_bounds[0]:
                        # Calculate average temperature, salinity, and layer
                        # thickness for this 3D triangular prism
                        temp_vals = []
                        salt_vals = []
                        dz_vals = []
                        for i in range(3):
                            temp_vals.append(temp[nodes[i].id])
                            salt_vals.append(salt[nodes[i].id])
                            dz_vals.append(abs(nodes[i].depth - nodes[i].below.depth))
                        # Calculate volume
                        volume = area*mean(array(dz_vals))
                        # Integrate temperature, salinity, volume
                        temp_int += mean(array(temp_vals))*volume
                        salt_int += mean(array(salt_vals))*volume
                        volume_int += volume
                    # Get ready for next iteration of loop
                    for i in range(3):
                        nodes[i] = nodes[i].below
        # Convert temperature and salinity from integrals to volume-averages,
        # append to timeseries
        temp_avg.append(temp_int/volume_int)
        salt_avg.append(salt_int/volume_int)

    time = arange(len(temp_avg))

    print('Plotting')
    fig, ax1 = subplots()
    # Temperature
    ax1.plot(time, temp_avg, color='b')
    ax1.set_ylabel(r'Average temperature ($^{\circ}$C)', color='b')
    for t1 in ax1.get_yticklabels():
        t1.set_color('b')
    ax1.set_xlabel('Years')
    ax1.grid(True, axis='x')
    ax2 = ax1.twinx()
    # Salinity
    ax2.plot(time, salt_avg, color='r')
    ax2.set_ylabel('Average salinity (psu)', color='r')
    for t2 in ax2.get_yticklabels():
        t2.set_color('r')
    title(lat_min_string+' to '+lat_max_string+', '+str(depth_bounds[0])+'m to '+str(depth_bounds[1])+'m')
    fig.savefig(fig_name)

    print('Saving results to log file')
    f = open(log_file, 'w')
    f.write('Average temperature (C):\n')
    for elm in temp_avg:
        f.write(str(elm) + '\n')
    f.write('Average salinity (psu):\n')
    for elm in salt_avg:
        f.write(str(elm) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    mesh_path = input("Path to FESOM mesh directory: ")
    output_path = input("Path to FESOM output directory: ")
    lat_min = float(input("Southern bound on latitude (-90 to 90): "))
    lat_max = float(input("Northern bound on latitude (-90 to 90): "))
    lat_bounds = [lat_min, lat_max]
    depth_min = float(input("Shallow bound on depth (positive, in metres): "))
    depth_max = float(input("Deep bound on depth (positive, in metres): "))
    depth_bounds = [depth_min, depth_max]
    start_year = int(input("First year to process: "))
    end_year = int(input("Last year to process: "))
    log_file = input("Path to logfile to save values and/or read previously calculated values: ")
    fig_name = input("Filename for figure: ")

    timeseries_drift (mesh_path, output_path, lat_bounds, depth_bounds, start_year, end_year, log_file, fig_name)

    
    
        

    

    
