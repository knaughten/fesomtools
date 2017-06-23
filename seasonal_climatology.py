from netCDF4 import Dataset
from numpy import *

# Create a seasonal climatology (DJF, MAM, JJA, SON) of ocean variables
# (3D temperature and salinity) and sea ice variables (2D concentration and
# effective thickness) over the simulation. Output to NetCDF files. 
# NB: This script assumes the FESOM output is in 5-day averages.
# Input:
# directory = path to FESOM model output directory containing one oce.mean.nc
#             file and one ice.mean.nc file for each year
# start_year, end_year = integers containing range of years to process
# out_file_oce, out_file_ice = paths to desired NetCDF files to save the
#                              seasonal climatology of ocean and sea ice
#                              variables respectively
def seasonal_climatology (directory, start_year, end_year, out_file_oce, out_file_ice):

    # Filename head
    expt_name = 'MK44005'

    print 'Processing ocean'
    # Read number of 3D nodes from first file
    id = Dataset(directory + expt_name + '.' + str(start_year) + '.oce.mean.nc', 'r')
    n3d = id.variables['temp'].shape[1]
    id.close()
    # Set up arrays to integrate seasonal climatology of temp and salt
    seasonal_temp = ma.empty([4, n3d]) 
    seasonal_salt = ma.empty([4, n3d])
    seasonal_temp[:,:] = 0.0
    seasonal_salt[:,:] = 0.0
    # Also integrate number of days in each season
    ndays = zeros(4)
    # Loop over years
    for year in range(start_year, end_year+1):
        print '...' + str(year)
        id = Dataset(directory + expt_name + '.' + str(year) + '.oce.mean.nc', 'r')
        # Indices 1-11 and 4/5 of index 12 are DJF (59 days)
        seasonal_temp[0,:] += sum(id.variables['temp'][0:11,:]*5, axis=0) + id.variables['temp'][11,:]*4
        seasonal_salt[0,:] += sum(id.variables['salt'][0:11,:]*5, axis=0) + id.variables['salt'][11,:]*4
        ndays[0] += 59
        # 1/5 of index 12, indices 13-30, and 1/5 of index 31 are MAM (92 days)
        seasonal_temp[1,:] += id.variables['temp'][11,:] + sum(id.variables['temp'][12:30,:]*5, axis=0) + id.variables['temp'][30,:]
        seasonal_salt[1,:] += id.variables['salt'][11,:] + sum(id.variables['salt'][12:30,:]*5, axis=0) + id.variables['salt'][30,:]
        ndays[1] += 92
        # 4/5 of index 31, indices 32-48, and 3/5 of index 49 are JJA (92 days)
        seasonal_temp[2,:] += id.variables['temp'][30,:]*4 + sum(id.variables['temp'][31:48]*5, axis=0) + id.variables['temp'][48,:]*3
        seasonal_salt[2,:] += id.variables['salt'][30,:]*4 + sum(id.variables['salt'][31:48]*5, axis=0) + id.variables['salt'][48,:]*3
        ndays[2] += 92
        # 2/5 of index 49, indices 50-66, and 4/5 of index 67 are SON (91 days)
        seasonal_temp[3,:] += id.variables['temp'][48,:]*2 + sum(id.variables['temp'][49:66,:]*5, axis=0) + id.variables['temp'][66,:]*4
        seasonal_salt[3,:] += id.variables['salt'][48,:]*2 + sum(id.variables['salt'][49:66,:]*5, axis=0) + id.variables['salt'][66,:]*4
        ndays[3] += 91
        # 1/5 of index 67 and indices 68-73 are DJF again (31 days)
        seasonal_temp[0,:] += id.variables['temp'][66,:] + sum(id.variables['temp'][67:73,:]*5, axis=0)
        seasonal_salt[0,:] += id.variables['salt'][66,:] + sum(id.variables['salt'][67:73,:]*5, axis=0)
        ndays[0] += 31
        id.close()
    # Convert from sums to averages
    for season in range(4):
        seasonal_temp[season,:] = seasonal_temp[season,:]/ndays[season]
        seasonal_salt[season,:] = seasonal_salt[season,:]/ndays[season]
    # Write to file
    print 'Writing ' + out_file_oce
    id = Dataset(out_file_oce, 'w')
    id.createDimension('nodes_3d', n3d)
    id.createDimension('T', None)
    id.createVariable('season', 'f8', ('T'))
    id.variables['season'].long_name = 'DJF, MAM, JJA, SON'
    id.variables['season'][:] = arange(1,4+1)
    id.createVariable('temp', 'f8', ('T', 'nodes_3d'))
    id.variables['temp'].description = 'mean potential temperature'
    id.variables['temp'].units = 'degC'
    id.variables['temp'][:,:] = seasonal_temp
    id.createVariable('salt', 'f8', ('T', 'nodes_3d'))
    id.variables['salt'].description = 'mean salinity'
    id.variables['salt'].units = 'psu'
    id.variables['salt'][:,:] = seasonal_salt
    id.close()

    print 'Processing sea ice'
    # Similar, but 2D nodes not 3D, sea ice area and effective thickness
    id = Dataset(directory + expt_name + '.' + str(start_year) + '.ice.mean.nc', 'r')
    n2d = id.variables['area'].shape[1]
    id.close()
    seasonal_area = ma.empty([4, n2d]) 
    seasonal_hice = ma.empty([4, n2d])
    seasonal_area[:,:] = 0.0
    seasonal_hice[:,:] = 0.0
    ndays = zeros(4)
    for year in range(start_year, end_year+1):
        print '...' + str(year)
        id = Dataset(directory + expt_name + '.' + str(year) + '.ice.mean.nc', 'r')
        # DJF
        seasonal_area[0,:] += sum(id.variables['area'][0:11,:]*5, axis=0) + id.variables['area'][11,:]*4
        seasonal_hice[0,:] += sum(id.variables['hice'][0:11,:]*5, axis=0) + id.variables['hice'][11,:]*4
        ndays[0] += 59
        # MAM
        seasonal_area[1,:] += id.variables['area'][11,:] + sum(id.variables['area'][12:30,:]*5, axis=0) + id.variables['area'][30,:]
        seasonal_hice[1,:] += id.variables['hice'][11,:] + sum(id.variables['hice'][12:30,:]*5, axis=0) + id.variables['hice'][30,:]
        ndays[1] += 92
        # JJA
        seasonal_area[2,:] += id.variables['area'][30,:]*4 + sum(id.variables['area'][31:48]*5, axis=0) + id.variables['area'][48,:]*3
        seasonal_hice[2,:] += id.variables['hice'][30,:]*4 + sum(id.variables['hice'][31:48]*5, axis=0) + id.variables['hice'][48,:]*3
        ndays[2] += 92
        # SON
        seasonal_area[3,:] += id.variables['area'][48,:]*2 + sum(id.variables['area'][49:66,:]*5, axis=0) + id.variables['area'][66,:]*4
        seasonal_hice[3,:] += id.variables['hice'][48,:]*2 + sum(id.variables['hice'][49:66,:]*5, axis=0) + id.variables['hice'][66,:]*4
        ndays[3] += 91
        # DJF again
        seasonal_area[0,:] += id.variables['area'][66,:] + sum(id.variables['area'][67:73,:]*5, axis=0)
        seasonal_hice[0,:] += id.variables['hice'][66,:] + sum(id.variables['hice'][67:73,:]*5, axis=0)
        ndays[0] += 31
        id.close()
    for season in range(4):
        seasonal_area[season,:] = seasonal_area[season,:]/ndays[season]
        seasonal_hice[season,:] = seasonal_hice[season,:]/ndays[season]
    print 'Writing ' + out_file_ice
    id = Dataset(out_file_ice, 'w')
    id.createDimension('nodes_2d', n2d)
    id.createDimension('T', None)
    id.createVariable('season', 'f8', ('T'))
    id.variables['season'].long_name = 'DJF, MAM, JJA, SON'
    id.variables['season'][:] = arange(1,4+1)
    id.createVariable('area', 'f8', ('T', 'nodes_2d'))
    id.variables['area'].description = 'ice concentration [0 to 1]'
    id.variables['area'][:,:] = seasonal_area
    id.createVariable('hice', 'f8', ('T', 'nodes_2d'))
    id.variables['hice'].description = 'effective ice thickness'
    id.variables['hice'].units = 'm'
    id.variables['hice'][:,:] = seasonal_hice
    id.close()
    

# Command-line interface
if __name__ == "__main__":

    directory = raw_input("Path to FESOM output directory: ")
    start_year = int(raw_input("First year to process: "))
    end_year = int(raw_input("Last year to process: "))
    out_file_oce = raw_input("Path to desired output climatology file for ocean variables: ")
    out_file_ice = raw_input("Path to desired output climatology file for sea ice variables: ")
    seasonal_climatology(directory, start_year, end_year, out_file_oce, out_file_ice)
    
