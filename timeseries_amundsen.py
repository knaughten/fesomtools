from netCDF4 import Dataset
from numpy import *
from os.path import *
from fesom_grid import *

def timeseries_amundsen (mesh_path, ice_diag_file, log_file):

    # Mesh parameters
    circumpolar = True
    cross_180 = False
    # Number of days for each output step
    days_per_output = 5
    # Bounds on Amundsen Sea box
    lon_min = -115
    lon_max = -100
    lat_max = -71

    avg_ice2ocn = []
    # Check if the log file exists
    if exists(log_file):
        print('Reading previously calculated values')
        f = open(log_file, 'r')
        # Skip the first line (header)
        f.readline()
        for line in f:
            avg_ice2ocn.append(float(line))
        f.close()    

    print('Building grid')
    elements = fesom_grid(mesh_path, circumpolar, cross_180)
    num_elm2D = len(elements)

    print('Reading data')
    # Change sign on ice growth rate in m/s, multiply by 1e8 for visibility
    id = Dataset(ice_diag_file, 'r')
    ice2ocn = -1e8*id.variables['thdgr'][:,:]
    id.close()
    num_time = size(ice2ocn,0)

    print('Setting up arrays')
    # Location flag for non-cavity elements in Amundsen Sea
    location_flag = zeros(num_elm2D)
    # Area of each element in Amundsen Sea
    area_elm = zeros(num_elm2D)
    # Ice to ocean freshwater flux timeseries at each element
    ice2ocn_elm = zeros([num_time, num_elm2D])
    # Loop over each element to fill these in
    for i in range(num_elm2D):
        elm = elements[i]
        # Ignore ice shelf cavities
        if not elm.cavity:
            # Check if we're within the given lon and lat bounds
            if all(elm.lon >= lon_min) and all(elm.lon <= lon_max) and all(elm.lat <= lat_max):
                # Save area
                area_elm[i] = elm.area()
                # Set location flag
                location_flag[i] = 1
                # Average ice-ocean freshwater flux timeseries over 3 components
                ice2ocn_elm[:,i] = (ice2ocn[:,elm.nodes[0].id] + ice2ocn[:,elm.nodes[1].id] + ice2ocn[:,elm.nodes[2].id])/3.0

    print('Building timeseries')
    for t in range(num_time):
        # Average over area of the correct elements
        avg_ice2ocn.append(sum(ice2ocn_elm[t,:]*area_elm*location_flag)/sum(area_elm*location_flag))

    print('Saving results to log file')
    f = open(log_file, 'w')
    f.write('Average ice-to-ocean freshwater flux (1e-8 m/s): \n')
    for val in avg_ice2ocn:
        f.write(str(val) + '\n')


# Command-line interface
if __name__ == "__main__":

    mesh_path = input("Path to FESOM mesh directory: ")
    forcing_file = input("Path to FESOM ice.diag.nc file: ")
    log_file = input("Path to logfile to save values and/or read previously calculated values: ")
    timeseries_amundsen(mesh_path, ice_diag_file, log_file)
        
    
    

    
