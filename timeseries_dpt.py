from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from fesom_grid import *
from fesom_sidegrid import *
from unrotate_vector import *

# Calculate and plot timeseries of the Drake Passage transport during a FESOM
# simulation.
# Takes 8 GB memory on raijin for Kaitlin's low_res mesh.
# Input:
# mesh_path = path to FESOM mesh directory
# ocn_file = path to output oce.mean.nc, assumed to have 5-day averages
# log_file = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
# fig_dir = optional string containing directory to save figures into. Make
#           sure it ends with a "/". Default is an empty string.
def timeseries_dpt (mesh_path, ocn_file, log_file, fig_dir=''):

    circumpolar = False  # Needs to be global for SideElements
    cross_180 = False    # Don't make second copies of elements that cross 180E
    days_per_output = 5  # Number of days for each output step

    # Longitude of Drake Passage zonal slice
    lon0 = -67
    # Latitude bounds on Drake Passage zonal slice
    lat_min = -68
    lat_max = -54.5

    dpt = []
    # Check if the log file exists
    if exists(log_file):
        print('Reading previously calculated values')
        f = open(log_file, 'r')
        # Skip the first line (header)
        f.readline()
        for line in f:
            dpt.append(float(line))
        f.close()

    print('Building grid')
    # First get regular 2D elements
    elm2D = fesom_grid(mesh_path, circumpolar, cross_180)
    # Read longitude and latitude of each node in order (needed for rotation)
    fid = open(mesh_path + 'nod3d.out', 'r')
    fid.readline()
    lon = []
    lat = []
    for line in fid:
        tmp = line.split()
        lon_tmp = float(tmp[1])
        lat_tmp = float(tmp[2])
        if lon_tmp < -180:
            lon_tmp += 360
        elif lon_tmp > 180:
            lon_tmp -= 360
        lon.append(lon_tmp)
        lat.append(lat_tmp)
    fid.close()
    lon = array(lon)
    lat = array(lat)

    print('Reading data')
    id = Dataset(ocn_file, 'r')
    num_time = id.variables['time'].shape[0]
    # Read both u and v so we can rotate to get the real u
    u_r = id.variables['u'][:,:]
    v_r = id.variables['v'][:,:]
    id.close()

    print('Unrotating velocity vector')
    u = zeros(shape(u_r))
    # Rotate one time index at a time
    for t in range(num_time):
        u_tmp, v_tmp = unrotate_vector(lon, lat, u_r[t,:], v_r[t,:])
        u[t,:] = u_tmp

    print('Extracting zonal slice through Drake Passage')
    # Get quadrilateral elements in the latitude vs depth slice
    selements = fesom_sidegrid(elm2D, u, lon0, lat_max, lat_min)

    print('Setting up arrays')
    # Eastward velocity at each element
    u_selm = zeros([num_time, len(selements)])
    # Area of each element
    area_selm = zeros(len(selements))
    # Loop over elements to fill these in
    for i in range(len(selements)):
        selm = selements[i]
        u_selm[:,i] = selm.var
        area_selm[i] = selm.area()
    # Build timeseries
    for t in range(num_time):
        # Integrate u*area and convert to Sv
        dpt.append(sum(u_selm[t,:]*area_selm)*1e-6)

    # Calculate time values
    time = arange(len(dpt))*days_per_output/365.

    print('Plotting')
    clf()
    plot(time, dpt)
    xlabel('Years')
    ylabel('Drake Passage Transport (Sv)')
    grid(True)
    savefig(fig_dir + 'drakepsgtrans.png')

    print('Saving results to log file')
    f = open(log_file, 'w')
    f.write('Drake Passage Transport (Sv):\n')
    for elm in dpt:
        f.write(str(elm) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    mesh_path = input("Path to FESOM mesh directory: ")
    ocn_file = input("Path to FESOM oce.mean.nc output file: ")
    log_file = input("Path to logfile to save values and/or read previously calculated values: ")
    timeseries_dpt(mesh_path, ocn_file, log_file)
    
        

    
