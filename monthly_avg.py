from numpy import *
from netCDF4 import Dataset

# Calculate a monthly average of the given variable over the given month.
# This is hard-coded and ugly, because the FESOM output does not have a
# proper calendar attached to the time axis.
# Input:
# file_path = Path to FESOM output containing 1 year of 5-day averages
# var = string containing variable name
# month = month number (0 to 11)
# Output:
# monthly_data = array of size n containing monthly average of "var"
def monthly_avg (file_path, var, month):

    # Number of days per month
    ndays_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    id = Dataset(file_path, 'r')
    num_pts = id.variables[var].shape[1]

    if month == 0:
        # January: indices 1-6 (1-based) and 1/5 of index 7
        monthly_data = sum(id.variables[var][0:6,:]*5, axis=0) + id.variables[var][6,:]
    elif month == 1:
        # Feburary: 4/5 of index 7, indices 8-11, and 4/5 of index 12
        monthly_data = id.variables[var][6,:]*4 + sum(id.variables[var][7:11,:]*5, axis=0) + id.variables[var][11,:]*4
    elif month == 2:
        # March: 1/5 of index 12 and indices 13-18
        monthly_data = id.variables[var][12,:] + sum(id.variables[var][12:18,:]*5, axis=0)
    elif month == 3:
        # April: indices 19-24
        monthly_data = sum(id.variables[var][18:24,:]*5, axis=0)
    elif month == 4:
        # May: indices 25-30, 1/5 of index 31
        monthly_data = sum(id.variables[var][24:30,:]*5, axis=0) + id.variables[var][30,:]
    elif month == 5:
        # June: 4/5 of index 31, indices 32-36, 1/5 of index 37
        monthly_data = id.variables[var][30,:]*4 + sum(id.variables[var][31:36,:]*5, axis=0) + id.variables[var][36,:]
    elif month == 6:
        # July: 4/5 of index 37, indices 38-42, 2/5 of index 43
        monthly_data = id.variables[var][36,:]*4 + sum(id.variables[var][37:42,:]*5, axis=0) + id.variables[var][42,:]*2
    elif month == 7:
        # August: 3/5 of index 43, indices 44-48, 3/5 of index 49
        monthly_data = id.variables[var][42,:]*3 + sum(id.variables[var][43:48,:]*5, axis=0) + id.variables[var][48,:]*3
    elif month == 8:
        # September: 2/5 of index 49, indices 50-54, 3/5 of index 55
        monthly_data = id.variables[var][48,:]*2 + sum(id.variables[var][49:54,:]*5, axis=0) + id.variables[var][54,:]*3
    elif month == 9:
        # October: 2/5 of index 55, indices 56-60, 4/5 of index 61
        monthly_data = id.variables[var][54,:]*2 + sum(id.variables[var][55:60,:]*5, axis=0) + id.variables[var][60,:]*4
    elif month == 10:
        # November: 1/5 of index 61, indices 62-66, 4/5 of index 67
        monthly_data = id.variables[var][60,:] + sum(id.variables[var][61:66,:]*5, axis=0) + id.variables[var][66,:]*4
    elif month == 11:
        # December: 1/5 of index 67, indices 68-73
        monthly_data = id.variables[var][66,:] + sum(id.variables[var][67:73,:]*5, axis=0)
    id.close()
    # Convert from sum to average
    monthly_data /= ndays_month[month]

    return monthly_data
