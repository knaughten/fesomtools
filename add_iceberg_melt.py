from netCDF4 import Dataset
from numpy import *
from scipy.interpolate import griddata

# Read Martin and Adcroft's monthly climatology of freshwater fluxes
# from iceberg melt, and add to the precipitation fields used by
# FESOM. I suppose I should technically use a separate field such
# as runoff, but this is easier and matches ROMS!
# Input:
# file = path to ERA-Interim FC forcing file, containing one year of
#        monthly averages for total precipitation ("tp") in m/12h
def add_iceberg_melt (file):

    # Naming conventions for iceberg files
    iceberg_head = '/short/m68/kaa561/ROMS-CICE-MCT/data/MartinAdcroft2010_iceberg_meltfluxes/icebergs.1861-1960.'
    iceberg_tail = '.melt.nc'
    # Iceberg grid information
    iceberg_grid = '/short/m68/kaa561/ROMS-CICE-MCT/data/MartinAdcroft2010_iceberg_meltfluxes/icebergs.static.nc'
    # File containing ERA-Interim grid
    era_grid = '/short/y99/kaa561/FESOM/ERA_Interim_monthly/FC_1995_monthly_orig.nc'
    # Density of freshwater
    rho_w = 1e3
    # Seconds in 12 hours
    seconds_per_12h = 60.*60.*12.

    # Read ERA-Interim grid
    id = Dataset(era_grid, 'r')
    lon_era_1d = id.variables['longitude'][:]
    lat_era_1d = id.variables['latitude'][:]
    id.close()
    # Get a 2D mesh of lon and lat
    lon_era = tile(lon_era_1d, (size(lat_era_1d),1))
    lat_era = transpose(tile(lat_era_1d, (size(lon_era_1d),1)))

    # Read the iceberg grid
    id = Dataset(iceberg_grid, 'r')
    lon_iceberg = id.variables['lon'][:,:]
    lat_iceberg = id.variables['lat'][:,:]
    id.close()

    # Make sure longitudes are between 0 and 360
    index = lon_iceberg < 0
    lon_iceberg[index] = lon_iceberg[index] + 360

    # Loop over months
    for month in range(12):
        print 'Processing month ' + str(month+1)
        # Reconstruct the filename of this month's iceberg data
        if month+1 < 10:
            month_str = '0' + str(month+1)
        else:
            month_str = str(month+1)
        iceberg_file = iceberg_head + month_str + iceberg_tail
        # Read iceberg freshwater flux in kg/m^2/s
        id = Dataset(iceberg_file, 'r')
        melt_iceberg = id.variables['melt'][0,:,:]
        id.close()
        # Interpolate to ERA-Interim grid
        melt_era = interp_iceberg2era(melt_iceberg, lon_iceberg, lat_iceberg, lon_era, lat_era)
        # Convert to m per 12 h
        melt_era = melt_era/rho_w*seconds_per_12h
        # Add to precipitation field for this month
        id = Dataset(file, 'a')
        id.variables['tp'][month,:,:] = id.variables['tp'][month,:,:] + melt_era
        id.close()


# Given a field A on the iceberg grid, linearly interpolate to the
# ERA-Interim grid.
# Input:
# A = 2D array (m x n) containing data on the iceberg grid
# lon_iceberg = 2D array (m x n) containing longitude values on the
#               iceberg grid, from 0 to 360
# lat_iceberg = 2D array (m x n) containing latitude values on the
#               iceberg grid
# lon_era = 2D array (p x q) containing longitude values on the ERA-Interim
#           grid, from 0 to 360
# lat_era = 2D array (p x q) containing latitude values on the ERA-Interim
#           grid
# Output:
# A_interp = 2D array (p x q) containing A interpolated to the ERA-Interim
#            grid
def interp_iceberg2era (A, lon_iceberg, lat_iceberg, lon_era, lat_era):

    # Set up an nx2 array containing the coordinates of each point in the
    # iceberg grid
    points = empty([size(lon_iceberg), 2])
    points[:,0] = ravel(lon_iceberg)
    points[:,1] = ravel(lat_iceberg)
    # Also flatten the data
    values = ravel(A)
    # Now set up an mx2 array containing the coordinates of each point
    # we want to interpolate to, in the ROMS grid
    xi = empty([size(lon_era), 2])
    xi[:,0] = ravel(lon_era)
    xi[:,1] = ravel(lat_era)
    # Now call griddata; fill out-of-bounds values (such as under ice shelves)
    # with zeros
    result = griddata(points, values, xi, method='linear', fill_value=0.0)
    # Un-flatten the result
    A_interp = reshape(result, shape(lon_era))

    return A_interp


# Command-line interface
if __name__ == "__main__":

    file = raw_input('Path to ERA-Interim FC file containing one year of monthly data: ')
    add_iceberg_melt(file)
