from netCDF4 import Dataset
from numpy import *

# Calculate a climatology of the given type of FESOM output file (eg oce.mean)
# over the given years. It will have the same frequency as the existing files
# (eg 5-day averages). Save to a new NetCDF file.
# Input:
# directory = path to directory containing FESOM output files
# file_type = string containing type of FESOM output file to process:
#             "forcing.diag", "ice.diag", "ice.mean", or "oce.mean"
# start_year, end_year = integers containing years to calculate climatology over
# out_file = path to desired climatology file
def average_years (directory, file_type, start_year, end_year, out_file):

    # Filename head
    expt_name = 'MK44005'

    # Variables to process, depending on file type
    if file_type == 'forcing.diag':
        vars = ['tair', 'tdew', 'shum', 'uwind', 'vwind', 'rain', 'snow', 'runoff', 'evap', 'lwrd', 'swrd', 'qnet', 'olat', 'osen', 'olwout', 'wnet', 'virtual_salt', 'relax_salt', 'stress_x', 'stress_y']
    elif file_type == 'ice.diag':
        vars = ['thdgr', 'thdgrsn', 'uhice', 'vhice', 'uhsnow', 'vhsnow', 'flice']
    elif file_type == 'ice.mean':
        vars = ['area', 'hice', 'hsnow', 'uice', 'vice']
    elif file_type == 'oce.mean':
        vars = ['ssh','u', 'v', 'w', 'temp', 'salt', 'ptr1']
        # For oce.mean, some variables are 2D and some are 3D
        is_3d = [False, True, True, True, True, True, True]
    else:
        print('Invalid file type')
        return

    print('Setting up file and processing year ' + str(start_year))
    # Look at the first variable in the first file to figure out sizes of each
    # dimension
    id_in = Dataset(directory + expt_name + '.' + str(start_year) + '.' + file_type + '.nc', 'r')
    dimensions = id_in.variables[vars[0]].shape
    # First dimension is time
    num_time = dimensions[0]
    if file_type == 'oce.mean':
        if is_3d[0]:
            # First variable was a 3D variable
            # Save the number of 3D nodes
            num_3d = dimensions[1]
            # Find a 2D variable
            for i in range(1, len(is_3d)):
                if not is_3d[i]:
                    # Save the number of 2D nodes
                    dimensions2 = id_in.variables[vars[i]].shape
                    num_2d = dimensions2[1]
                    break
        else:
            # First variable was a 2D variable
            # Save the number of 2D nodes
            num_2d = dimensions[1]
            # Find a 3D variable
            for i in range(1, len(vars)):
                if is_3d[i]:
                    # Save the number of 3D nodes
                    dimensions3 = id_in.variables[vars[i]].shape
                    num_3d = dimensions3[1]
    else:
        # All 2D variables
        num_2d = dimensions[1]
    # Save time values
    time = id_in.variables['time'][:]

    # Set up output file
    id_out = Dataset(out_file, 'w')
    id_out.createDimension('nodes_2d', num_2d)
    if file_type == 'oce.mean':
        id_out.createDimension('nodes_3d', num_3d)
    id_out.createDimension('T', None)
    id_out.createVariable('time', 'f8', ('T'))
    id_out.variables['time'].long_name = 'model time'
    id_out.variables['time'].units = 's'
    id_out.variables['time'][:] = time
    # Loop over variables
    for i in range(len(vars)):
        # Define variable
        if file_type == 'oce.mean':
            # Figure out if 2D or 3D
            if is_3d[i]:
                id_out.createVariable(vars[i], 'f8', ('T', 'nodes_3d'))
            else:
                id_out.createVariable(vars[i], 'f8', ('T', 'nodes_2d'))
        else:
            # Definitely 2D
            id_out.createVariable(vars[i], 'f8', ('T', 'nodes_2d'))
        # Carry forward any attributes
        if 'description' in id_in.variables[vars[i]].ncattrs():
            id_out.variables[vars[i]].description = id_in.variables[vars[i]].description
        if 'units' in id_in.variables[vars[i]].ncattrs():
            id_out.variables[vars[i]].units = id_in.variables[vars[i]].units
        # Copy the data for the first year
        id_out.variables[vars[i]][:,:] = id_in.variables[vars[i]][:,:]
    id_in.close()

    # Loop over the other years
    for year in range(start_year+1, end_year+1):
        print('Processing year ' + str(year))
        id_in = Dataset(directory + expt_name + '.' + str(year) + '.' + file_type + '.nc', 'r')
        for var_name in vars:
            # Accumulate the data in the output file
            id_out.variables[var_name][:,:] = id_out.variables[var_name][:,:] + id_in.variables[var_name][:,:]
        id_in.close()

    # Convert from sums to averages (over each timestep in climatology)
    num_years = end_year - start_year + 1
    for var_name in vars:
        id_out.variables[var_name][:,:] = id_out.variables[var_name][:,:]/num_years
    id_out.close()
        
    
# Command-line interface
if __name__ == "__main__":

    directory = input("Path to model output directory: ")
    file_type = input("File type (forcing.diag, ice.diag, ice.mean, or oce.mean): ")
    start_year = int(input("Starting year for averages: "))
    end_year = int(input("Ending year for averages: "))
    out_file = input("Path to desired output file: ")
    average_years(directory, file_type, start_year, end_year, out_file)
