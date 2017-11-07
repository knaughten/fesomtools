from numpy import *
from netCDF4 import Dataset

def monthly_climatology_temp (directory, start_year, end_year, out_file):

    # Filename head
    expt_name = 'MK44005'
    # Number of days per month
    # Note leap days are disregarded in output
    ndays_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    num_years = end_year-start_year+1

    # Read number of 3D nodes from first file
    id = Dataset(directory + expt_name + '.' + str(start_year) + '.oce.mean.nc', 'r')
    n3d = id.variables['temp'].shape[1]
    id.close()
    # Set up array to integrate monthly climatology of temperature
    monthly_temp = ma.empty([12, n3d])
    monthly_temp[:,:] = 0.0
    # Also integrate number of days in each month
    ndays = zeros(12)
    # Loop over years
    for year in range(start_year, end_year+1):
        print '...' + str(year)
        id = Dataset(directory + expt_name + '.' + str(year) + '.oce.mean.nc', 'r')
        # January: indices 1-6 (1-based) and 1/5 of index 7
        monthly_temp[0,:] += sum(id.variables['temp'][0:6,:]*5, axis=0) + id.variables['temp'][6,:]
        # Feburary: 4/5 of index 7, indices 8-11, and 4/5 of index 12
        monthly_temp[1,:] += id.variables['temp'][6,:]*4 + sum(id.variables['temp'][7:11,:]*5, axis=0) + id.variables['temp'][11,:]*4
        # March: 1/5 of index 12 and indices 13-18
        monthly_temp[2,:] += id.variables['temp'][11,:] + sum(id.variables['temp'][12:18,:]*5, axis=0)
        # April: indices 19-24
        monthly_temp[3,:] += sum(id.variables['temp'][18:24,:]*5, axis=0)
        # May: indices 25-30, 1/5 of index 31
        monthly_temp[4,:] += sum(id.variables['temp'][24:30,:]*5, axis=0) + id.variables['temp'][30,:]
        # June: 4/5 of index 31, indices 32-36, 1/5 of index 37
        monthly_temp[5,:] += id.variables['temp'][30,:]*4 + sum(id.variables['temp'][31:36,:]*5, axis=0) + id.variables['temp'][36,:]
        # July: 4/5 of index 37, indices 38-42, 2/5 of index 43
        monthly_temp[6,:] += id.variables['temp'][36,:]*4 + sum(id.variables['temp'][37:42,:]*5, axis=0) + id.variables['temp'][42,:]*2
        # August: 3/5 of index 43, indices 44-48, 3/5 of index 49
        monthly_temp[7,:] += id.variables['temp'][42,:]*3 + sum(id.variables['temp'][43:48,:]*5, axis=0) + id.variables['temp'][48,:]*3
        # September: 2/5 of index 49, indices 50-54, 3/5 of index 55
        monthly_temp[8,:] += id.variables['temp'][48,:]*2 + sum(id.variables['temp'][49:54,:]*5, axis=0) + id.variables['temp'][54,:]*3
        # October: 2/5 of index 55, indices 56-60, 4/5 of index 61
        monthly_temp[9,:] += id.variables['temp'][54,:]*2 + sum(id.variables['temp'][55:60,:]*5, axis=0) + id.variables['temp'][60,:]*4
        # November: 1/5 of index 61, indices 62-66, 4/5 of index 67
        monthly_temp[10,:] += id.variables['temp'][60,:] + sum(id.variables['temp'][61:66,:]*5, axis=0) + id.variables['temp'][66,:]*4
        # December: 1/5 of index 67, indices 68-73
        monthly_temp[11,:] += id.variables['temp'][66,:] + sum(id.variables['temp'][67:73,:]*5, axis=0)
    # Convert from sums to averages
    for month in range(12):
        monthly_temp[month,:] = monthly_temp[month,:]/(ndays_month[month]*num_years)
    # Write to file
    print 'Writing ' + out_file
    id = Dataset(out_file, 'w')
    id.createDimension('nodes_3d', n3d)
    id.createDimension('T', None)
    id.createVariable('month', 'f8', ('T'))
    id.variables['month'][:] = arange(1,12+1)
    id.createVariable('temp', 'f8', ('T', 'nodes_3d'))
    id.variables['temp'].description = 'mean potential temperature'
    id.variables['temp'].units = 'degC'
    id.variables['temp'][:,:] = monthly_temp
    id.close()


# Command-line interface
if __name__ == "__main__":

    directory = raw_input("Path to FESOM output directory: ")
    start_year = int(raw_input("First year to process: "))
    end_year = int(raw_input("Last year to process: "))
    out_file = raw_input("Path to desired output climatology file: ")
    monthly_climatology_temp(directory, start_year, end_year, out_file)
