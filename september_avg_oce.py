from netCDF4 import Dataset
from numpy import *

# Calculate a September climatology of temperature and salinity over the given
# years. Save to a new NetCDF file.
# Input:
# directory = path to directory containing FESOM output files
# start_year, end_year = integers containing years to calculate climatology over
# out_file = path to desired climatology file
def september_avg_oce (directory, start_year, end_year, out_file):

    # Filename head
    expt_name = 'MK44005'

    # Read number of 3D nodes from first file
    id = Dataset(directory + expt_name + '.' + str(start_year) + '.oce.mean.nc', 'r')
    n3d = id.variables['temp'].shape[1]
    id.close()
    # Set up arrays to integrate temperature and salinity
    sept_temp = zeros(n3d)
    sept_salt = zeros(n3d)
    # Also integrate the number of days
    ndays = 0
    # Loop over years
    for year in range(start_year, end_year+1):
        print '...' + str(year)
        id = Dataset(directory + expt_name + '.' + str(year) + '.oce.mean.nc', 'r')
        # September: 2/5 of index 49, indices 50-54, 3/5 of index 55
        sept_temp += id.variables['temp'][48,:]*2 + sum(id.variables['temp'][49:54,:]*5, axis=0) + id.variables['temp'][54,:]*3
        sept_salt += id.variables['salt'][48,:]*2 + sum(id.variables['salt'][49:54,:]*5, axis=0) + id.variables['salt'][54,:]*3
        ndays += 30
        id.close()
    # Convert from sums to averages
    sept_temp /= ndays
    sept_salt /= ndays
    # Write to file
    print 'Writing ' + out_file
    id = Dataset(out_file, 'w')
    id.createDimension('nodes_3d', n3d)
    id.createDimension('T', None)
    id.createVariable('temp', 'f8', ('T', 'nodes_3d'))
    id.variables['temp'].description = 'mean potential temperature'
    id.variables['temp'].units = 'degC'
    id.variables['temp'][0,:] = sept_temp
    id.createVariable('salt', 'f8', ('T', 'nodes_3d'))
    id.variables['salt'].description = 'mean salinity'
    id.variables['salt'].units = 'psu'
    id.variables['salt'][0,:] = sept_salt
    id.close()
    

# Command-line interface
if __name__ == "__main__":

    directory = raw_input("Path to model output directory: ")
    start_year = int(raw_input("Starting year for averages: "))
    end_year = int(raw_input("Ending year for averages: "))
    out_file = raw_input("Path to desired output file: ")
    september_avg_oce(directory, start_year, end_year, out_file)
