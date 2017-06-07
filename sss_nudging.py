from netCDF4 import Dataset
from numpy import *
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import KDTree

# Interpolate the World Ocean Atlas 2013 monthly climatology of sea surface
# salinity to the ERA-Interim grid so it can be used as a forcing file by
# FESOM, for sea surface salinity restoring.
def sss_nudging ():

    # Paths to reconstruct WOA 2013 monthly files
    woa_head = '/short/m68/kaa561/metroms_iceshelf/data/originals/WOA_2013/woa13_95A4_s'
    woa_tail = '_01v2.nc'
    # File containing ERA-Interim grid
    era_grid = '/short/y99/kaa561/FESOM/ERA_Interim_1992_2016/icebergs.nc'
    # Path to output file
    out_file = '../FESOM/ERA_Interim_1992_2016/sss_nudge.nc'

    # Read ERA-Interim grid
    id = Dataset(era_grid, 'r')
    lon_era_1d = id.variables['longitude'][:]
    lat_era_1d = id.variables['latitude'][:]
    id.close()
    num_lon = size(lon_era_1d)
    num_lat = size(lat_era_1d)
    # Get a 2D mesh of lon and lat
    lon_era = tile(lon_era_1d, (num_lat,1))
    lat_era = transpose(tile(lat_era_1d, (num_lon,1)))

    # Make sure longitude goes from -180 to 180
    index = lon_era > 180
    lon_era[index] -= 360

    # Read WOA grid
    print 'Reading WOA grid'
    id = Dataset(woa_head + '01' + woa_tail, 'r')
    lon_woa_raw = id.variables['lon'][:]
    lat_woa_raw = id.variables['lat'][:]
    id.close()

    # Make sure longitude goes from -180 to 180
    index = lon_woa_raw > 180
    lon_woa_raw[index] -= 360    

    # The WOA longitude axis doesn't wrap around; there is a gap between
    # almost-180W and almost-180E, and the ERA grid has points in this gap.
    # So copy the last longitude value (mod 360) to the beginning, and the
    # first longitude value (mod 360) to the end.
    lon_woa = zeros(size(lon_woa_raw)+2)
    lon_woa[0] = lon_woa_raw[-1]-360
    lon_woa[1:-1] = lon_woa_raw
    lon_woa[-1] = lon_woa_raw[0]+360

    # Also duplicate the southernmost and northernmost rows so they cover 90S-90N
    lat_woa = zeros(size(lat_woa_raw)+2)
    lat_woa[1:-1] = lat_woa_raw
    lat_woa[0] = lat_woa[1]-1
    lat_woa[-1] = lat_woa[-2]+1

    # Set up output file
    print "Setting up " + out_file
    out_id = Dataset(out_file, 'w')
    # Define dimensions
    out_id.createDimension('longitude', num_lon)
    out_id.createDimension('latitude', num_lat)
    out_id.createDimension('time', None)
    # Define variables
    out_id.createVariable('longitude', 'f8', ('longitude'))
    out_id.variables['longitude'].units = 'degrees_east'
    out_id.variables['longitude'].long_name = 'longitude'
    out_id.variables['longitude'][:] = lon_era_1d
    out_id.createVariable('latitude', 'f8', ('latitude'))
    out_id.variables['latitude'].units = 'degrees_north'
    out_id.variables['latitude'].long_name = 'latitude'
    out_id.variables['latitude'][:] = lat_era_1d
    out_id.createVariable('time', 'f8', ('time'))
    out_id.variables['time'].units = 'months'
    out_id.createVariable('SALT', 'f8', ('time', 'latitude', 'longitude'))
    out_id.variables['SALT'].long_name = 'sea surface salinity'
    out_id.variables['SALT'].units = 'psu'

    # Loop over months
    for month in range(12):
        print 'Processing month ' + str(month+1) + ' of 12'
        # Construct file path
        if month+1 < 10:
            filename = woa_head + '0' + str(month+1) + woa_tail
        else:
            filename = woa_head + str(month+1) + woa_tail
        # Read surface salinity
        id = Dataset(filename, 'r')
        sss_raw = transpose(id.variables['s_an'][0,0,:,:])
        id.close()
        # Wrap the longitude axis
        sss = ma.array(zeros((size(lon_woa), size(lat_woa))))
        sss[1:-1,1:-1] = ma.copy(sss_raw)
        sss[0,1:-1] = ma.copy(sss_raw[-1,:])
        sss[-1,1:-1] = ma.copy(sss_raw[0,:])
        # Copy the northernmost and southernmost rows
        sss[:,0] = sss[:,1]
        sss[:,-1] = sss[:,-2]
        # Find land mask with nearest neighbours
        j,i = mgrid[0:sss.shape[0], 0:sss.shape[1]]
        jigood = array((j[~sss.mask], i[~sss.mask])).T
        jibad = array((j[sss.mask], i[sss.mask])).T
        sss[sss.mask] = sss[~sss.mask][KDTree(jigood).query(jibad)[1]]
        # Build a function for linear interpolation
        interp_function = RegularGridInterpolator((lon_woa, lat_woa), sss, bounds_error=True)
        # Call this function for each point on the ERA-Interim grid
        sss_interp = interp_function((lon_era, lat_era))
        # Save to file
        out_id.variables['time'][month] = month+1
        out_id.variables['SALT'][month,:,:] = sss_interp
    out_id.close()


# Command-line interface
if __name__ == "__main__":

    sss_nudging()
        

    
    
