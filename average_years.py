from netCDF4 import Dataset
from numpy import *

def average_years (directory, file_type, start_year, end_year, out_file):

    expt_name = 'MK44005'

    if file_type == 'forcing.diag':
        vars = ['tair', 'tdew', 'shum', 'uwind', 'vwind', 'rain', 'snow', 'runoff', 'evap', 'lwrd', 'swrd', 'qnet', 'olat', 'osen', 'olwout', 'wnet', 'virtual_salt', 'relax_salt', 'stress_x', 'stress_y']
    elif file_type == 'ice.diag':
        vars = ['thdgr', 'thdgrsn', 'uhice', 'vhice', 'uhsnow', 'vhsnow', 'flice']
    elif file_type == 'ice.mean':
        vars = ['area', 'hice', 'hsnow', 'uice', 'vice']
    elif file_type == 'oce.mean':
        vars = ['ssh','u', 'v', 'w', 'temp', 'salt', 'ptr1']
        is_3d = [False, True, True, True, True, True, True]
    else:
        print 'Invalid file type'
        return

    print 'Setting up file and processing year ' + str(start_year)
    id_in = Dataset(directory + expt_name + '.' + str(start_year) + '.' + file_type + '.nc', 'r')
    dimensions = id_in.variables[vars[0]].shape
    num_time = dimensions[0]
    if file_type == 'oce.mean':
        if is_3d[0]:
            num_3d = dimensions[1]
            for i in range(1, len(is_3d)):
                if not is_3d[i]:
                    dimensions2 = id_in.variables[vars[i]].shape
                    num_2d = dimensions2[1]
                    break
        else:
            num_2d = dimensions[1]
            for i in range(1, len(vars)):
                if is_3d[i]:
                    dimensions3 = id_in.variables[vars[i]].shape
                    num_3d = dimensions3[1]
    else:
        num_2d = dimensions[1]
    time = id_in.variables['time'][:]

    id_out = Dataset(out_file, 'w')
    id_out.createDimension('nodes_2d', num_2d)
    if file_type == 'oce.mean':
        id_out.createDimension('nodes_3d', num_3d)
    id_out.createDimension('T', None)
    id_out.createVariable('time', 'f8', ('T'))
    id_out.variables['time'].long_name = 'model time'
    id_out.variables['time'].units = 's'
    id_out.variables['time'][:] = time
    for i in range(len(vars)):
        if file_type == 'oce.mean':
            if is_3d[i]:
                id_out.createVariable(vars[i], 'f8', ('T', 'nodes_3d'))
            else:
                id_out.createVariable(vars[i], 'f8', ('T', 'nodes_2d'))
        else:
            id_out.createVariable(vars[i], 'f8', ('T', 'nodes_2d'))
        if 'description' in id_in.variables[vars[i]].ncattrs():
            id_out.variables[vars[i]].description = id_in.variables[vars[i]].description
        if 'units' in id_in.variables[vars[i]].ncattrs():
            id_out.variables[vars[i]].units = id_in.variables[vars[i]].units
        id_out.variables[vars[i]][:,:] = id_in.variables[vars[i]][:,:]
    id_in.close()

    for year in range(start_year+1, end_year+1):
        print 'Processing year ' + str(year)
        id_in = Dataset(directory + expt_name + '.' + str(year) + '.' + file_type + '.nc', 'r')
        for var_name in vars:
            id_out.variables[var_name][:,:] = id_out.variables[var_name][:,:] + id_in.variables[var_name][:,:]
        id_in.close()

    num_years = end_year - start_year + 1
    for var_name in vars:
        id_out.variables[var_name][:,:] = id_out.variables[var_name][:,:]/num_years
    id_out.close()
        
    
if __name__ == "__main__":

    directory = raw_input("Path to model output directory: ")
    file_type = raw_input("File type (forcing.diag, ice.diag, ice.mean, or oce.mean): ")
    start_year = int(raw_input("Starting year for averages: "))
    end_year = int(raw_input("Ending year for averages: "))
    out_file = raw_input("Path to desired output file: ")
    average_years(directory, file_type, start_year, end_year, out_file)
