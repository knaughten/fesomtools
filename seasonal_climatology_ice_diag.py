from netCDF4 import Dataset
from numpy import *

def seasonal_climatology_ice_diag (directory, start_year, end_year, out_file):

    # Filename head
    expt_name = 'MK44005'

    # Read number of 2D nodes from first file
    id = Dataset(directory + expt_name + '.' + str(start_year) + '.ice.diag.nc', 'r')
    n2d = id.variables['thdgr'].shape[1]
    id.close()
    # Set up arrays to integrate seasonal climatology of temp and salt
    seasonal_thdgr = ma.empty([4, n2d])
    seasonal_uhice = ma.empty([4, n2d])
    seasonal_vhice = ma.empty([4, n2d])
    seasonal_flice = ma.empty([4, n2d])
    seasonal_thdgr[:,:] = 0.0
    seasonal_uhice[:,:] = 0.0
    seasonal_vhice[:,:] = 0.0
    seasonal_flice[:,:] = 0.0
    # Also integrate number of days in each season
    ndays = zeros(4)
    # Loop over years
    for year in range(start_year, end_year+1):
        print '...' + str(year)
        id = Dataset(directory + expt_name + '.' + str(year) + '.ice.diag.nc', 'r')
        # Indices 1-11 and 4/5 of index 12 are DJF (59 days)
        seasonal_thdgr[0,:] += sum(id.variables['thdgr'][0:11,:]*5, axis=0) + id.variables['thdgr'][11,:]*4
        seasonal_uhice[0,:] += sum(id.variables['uhice'][0:11,:]*5, axis=0) + id.variables['uhice'][11,:]*4
        seasonal_vhice[0,:] += sum(id.variables['vhice'][0:11,:]*5, axis=0) + id.variables['vhice'][11,:]*4
        seasonal_flice[0,:] += sum(id.variables['flice'][0:11,:]*5, axis=0) + id.variables['flice'][11,:]*4
        ndays[0] += 59
        # 1/5 of index 12, indices 13-30, and 1/5 of index 31 are MAM (92 days)
        seasonal_thdgr[1,:] += id.variables['thdgr'][11,:] + sum(id.variables['thdgr'][12:30,:]*5, axis=0) + id.variables['thdgr'][30,:]
        seasonal_uhice[1,:] += id.variables['uhice'][11,:] + sum(id.variables['uhice'][12:30,:]*5, axis=0) + id.variables['uhice'][30,:]
        seasonal_vhice[1,:] += id.variables['vhice'][11,:] + sum(id.variables['vhice'][12:30,:]*5, axis=0) + id.variables['vhice'][30,:]
        seasonal_flice[1,:] += id.variables['flice'][11,:] + sum(id.variables['flice'][12:30,:]*5, axis=0) + id.variables['flice'][30,:]
        ndays[1] += 92
        # 4/5 of index 31, indices 32-48, and 3/5 of index 49 are JJA (92 days)
        seasonal_thdgr[2,:] += id.variables['thdgr'][30,:]*4 + sum(id.variables['thdgr'][31:48]*5, axis=0) + id.variables['thdgr'][48,:]*3
        seasonal_uhice[2,:] += id.variables['uhice'][30,:]*4 + sum(id.variables['uhice'][31:48]*5, axis=0) + id.variables['uhice'][48,:]*3
        seasonal_vhice[2,:] += id.variables['vhice'][30,:]*4 + sum(id.variables['vhice'][31:48]*5, axis=0) + id.variables['vhice'][48,:]*3
        seasonal_flice[2,:] += id.variables['flice'][30,:]*4 + sum(id.variables['flice'][31:48]*5, axis=0) + id.variables['flice'][48,:]*3
        ndays[2] += 92
        # 2/5 of index 49, indices 50-66, and 4/5 of index 67 are SON (91 days)
        seasonal_thdgr[3,:] += id.variables['thdgr'][48,:]*2 + sum(id.variables['thdgr'][49:66,:]*5, axis=0) + id.variables['thdgr'][66,:]*4
        seasonal_uhice[3,:] += id.variables['uhice'][48,:]*2 + sum(id.variables['uhice'][49:66,:]*5, axis=0) + id.variables['uhice'][66,:]*4
        seasonal_vhice[3,:] += id.variables['vhice'][48,:]*2 + sum(id.variables['vhice'][49:66,:]*5, axis=0) + id.variables['vhice'][66,:]*4
        seasonal_flice[3,:] += id.variables['flice'][48,:]*2 + sum(id.variables['flice'][49:66,:]*5, axis=0) + id.variables['flice'][66,:]*4
        ndays[3] += 91
        # 1/5 of index 67 and indices 68-73 are DJF again (31 days)
        seasonal_thdgr[0,:] += id.variables['thdgr'][66,:] + sum(id.variables['thdgr'][67:73,:]*5, axis=0)
        seasonal_uhice[0,:] += id.variables['uhice'][66,:] + sum(id.variables['uhice'][67:73,:]*5, axis=0)
        seasonal_vhice[0,:] += id.variables['vhice'][66,:] + sum(id.variables['vhice'][67:73,:]*5, axis=0)
        seasonal_flice[0,:] += id.variables['flice'][66,:] + sum(id.variables['flice'][67:73,:]*5, axis=0)
        ndays[0] += 31
        id.close()
    # Convert from sums to averages
    for season in range(4):
        seasonal_thdgr[season,:] = seasonal_thdgr[season,:]/ndays[season]
        seasonal_uhice[season,:] = seasonal_uhice[season,:]/ndays[season]
        seasonal_vhice[season,:] = seasonal_vhice[season,:]/ndays[season]
        seasonal_flice[season,:] = seasonal_flice[season,:]/ndays[season]
    # Write to file
    print 'Writing ' + out_file
    id = Dataset(out_file, 'w')
    id.createDimension('nodes_2d', n2d)
    id.createDimension('T', None)
    id.createVariable('season', 'f8', ('T'))
    id.variables['season'].long_name = 'DJF, MAM, JJA, SON'
    id.variables['season'][:] = arange(1,4+1)
    id.createVariable('thdgr', 'f8', ('T', 'nodes_2d'))
    id.variables['thdgr'].description = 'thermodynamic growth rate of eff. ice thickness'
    id.variables['thdgr'].units = 'm/s'
    id.variables['thdgr'][:,:] = seasonal_thdgr
    id.createVariable('uhice', 'f8', ('T', 'nodes_2d'))
    id.variables['uhice'].description = 'zonal advective flux of eff. ice thickness'
    id.variables['uhice'].units = 'm.m/s'
    id.variables['uhice'][:,:] = seasonal_uhice
    id.createVariable('vhice', 'f8', ('T', 'nodes_2d'))
    id.variables['vhice'].description = 'meridional advective flux of eff. ice thickness'
    id.variables['vhice'].units = 'm.m/s'
    id.variables['vhice'][:,:] = seasonal_vhice
    id.createVariable('flice', 'f8', ('T', 'nodes_2d'))
    id.variables['flice'].description = 'rate of flooding snow to ice'
    id.variables['flice'].units = 'm/s'
    id.variables['flice'][:,:] = seasonal_flice
    id.close()


# Command-line interface
if __name__ == "__main__":

    directory = raw_input("Path to FESOM output directory: ")
    start_year = int(raw_input("First year to process: "))
    end_year = int(raw_input("Last year to process: "))
    out_file = raw_input("Path to desired output climatology file: ")
    seasonal_climatology_ice_diag(directory, start_year, end_year, out_file)
