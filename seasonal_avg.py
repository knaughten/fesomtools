from numpy import *
from netCDF4 import Dataset

# Calculate seasonal averages (DJF, MAM, JJA, SON) of the given variable.
# This is hard-coded and ugly, because the FESOM output does not have a
# proper calendar attached to the time axis.
# Input:
# file1 = Path to FESOM output file containing 1 year of 5-day averages
#         (December will be used)
# file2 = same as file1, but for the next year (Jan-Nov will be used)
#         If this is a climatology, set file2=file1
# var = string containing variable name
# Output:
# seasonal_data = array of size 4xn containing seasonal averages of "var"
#                 (DJF, MAM, JJA, SON)
def seasonal_avg (file1, file2, var):

    id = Dataset(file1, 'r')
    num_pts = id.variables[var].shape[1]
    seasonal_data = zeros([4, num_pts])

    # DJF: 1/5 of index 67 (1-based) and indices 68-73 in file1; 
    #      indices 1-11 and 4/5 of index 12 in file2; 90 days in total
    seasonal_data[0,:] = id.variables[var][66,:] + sum(id.variables[var][67:73,:]*5, axis=0)
    id.close()
    id = Dataset(file2, 'r')
    seasonal_data[0,:] += sum(id.variables[var][0:11,:]*5, axis=0) + id.variables[var][11,:]*4
    seasonal_data[0,:] /= 90

    # MAM: 1/5 of index 12, indices 13-30, and 1/5 of index 31 in file2;
    #      92 days in total
    seasonal_data[1,:] = id.variables[var][11,:] + sum(id.variables[var][12:30,:]*5, axis=0) + id.variables[var][30,:]
    seasonal_data[1,:] /= 92

    # JJA: 4/5 of index 31, indices 32-48, and 3/5 of index 49 in file2;
    #      92 days in total
    seasonal_data[2,:] = id.variables[var][30,:]*4 + sum(id.variables[var][31:48]*5, axis=0) + id.variables[var][48,:]*3
    seasonal_data[2,:] /= 92

    # SON: 2/5 of index 49, indices 50-66, and 4/5 of index 67 in file2;
    #      91 days in total
    seasonal_data[3,:] = id.variables[var][48,:]*2 + sum(id.variables[var][49:66,:]*5, axis=0) + id.variables[var][66,:]*4
    seasonal_data[3,:] /= 91
    id.close()

    return seasonal_data

    
    

    
