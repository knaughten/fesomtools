from netCDF4 import Dataset
from numpy import *
from unesco import *

# Given an ocean output file with temperature and salinity data, calculate
# density fields at each timestep using the 1980 UNESCO seawater equation of
# state. Save in a new file.
# Input:
# mesh_path = path to FESOM mesh directory
# oce_file = path to ocean output file
# rho_file = desired path to new density file
def make_density_file (mesh_path, oce_file, rho_file):

    # Read depth values for 3D nodes
    file = open(mesh_path + 'nod3d.out', 'r')
    file.readline()
    depth = []
    for line in file:
        tmp = line.split()
        depth.append(float(tmp[3]))
    file.close()

    # Divide by 10 to convert to pressure
    depth = array(depth)
    press = abs(depth)/10.0
    num_nodes = size(press)

    # Set up output file
    out_id = Dataset(rho_file, 'w')
    # Define dimensions
    out_id.createDimension('nodes_3d', num_nodes)
    out_id.createDimension('T', None)
    # Define variables
    out_id.createVariable('time', 'f8', ('T'))
    out_id.variables['time'].long_name = 'model time'
    out_id.variables['time'].units = 's'
    out_id.createVariable('rho', 'f8', ('T', 'nodes_3d'))
    out_id.variables['rho'].long_name = 'density'
    out_id.variables['rho'].units = 'kg/m^3'

    # Read time values from input file
    in_id = Dataset(oce_file, 'r')
    time = in_id.variables['time'][:]

    # Process each timestep individually to conserve memory
    for t in range(size(time)):
        print 'Processing timestep '+str(t+1)+' of '+str(size(time))
        # Set a new time value in the output file
        out_id.variables['time'][t] = time[t]
        # Read temperature and salinity (convert to float128 to prevent
        # overflow in UNESCO calculations)
        temp = ma.asarray(in_id.variables['temp'][t,:], dtype=float128)
        salt = ma.asarray(in_id.variables['salt'][t,:], dtype=float128)
        # Magic happens here
        rho = unesco(temp, salt, press)
        # Save the results for this timestep
        out_id.variables['rho'][t,:] = rho

    in_id.close()
    out_id.close()


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input('Path to FESOM mesh directory: ')
    oce_file = raw_input('Path to ocean output file: ')
    rho_file = raw_input('Desired path to new density file: ')
    make_density_file(mesh_path, oce_file, rho_file)

    
    

    
        
