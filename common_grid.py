from numpy import *
from numpy.linalg import inv
from scipy.interpolate import griddata
from monthly_avg import *
from unrotate_vector import *

# For various 2D fields, calculate the monthly climatology of FESOM output and
# interpolate to a regular grid (quarter-degree, circumpolar to 50S) for easy
# comparison with ROMS. The fields are: sea surface temperature and salinity,
# surface heat and salt fluxes, sea ice concentration and thickness, ocean
# surface velocity vector, sea ice velocity vector, surface stress vector,
# and curl of the surface stress.
# Input:
# mesh_path = path to FESOM mesh directory
# output_dir = path to directory containing FESOM output: 5-day averages for
#              each year between start_year and end_year
# start_year, end_year = years to calculate climatology over
# common_file = file containing land mask on the common grid (interpolated from
#               ROMS, see common_grid.py in roms_tools GitHub repo)
# out_file = path to desired output file
def common_grid (mesh_path, output_dir, start_year, end_year, common_file, out_file):

    # Resolution of common grid (degrees, same for lat and lon)
    res = 0.25
    # Northern boundary to interpolate to
    nbdry = -50
    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Name stamped on FESOM output files
    expt_name = 'MK44005'

    print('Calculating grids')
    # Make the latitude and longitude arrays for the common grid
    lon_common = arange(-180, 180+res, res)
    lat_common = arange(-90, nbdry+res, res)
    # Get a 2D version of each to calculate dx and dy in metres
    lon_2d, lat_2d = meshgrid(lon_common, lat_common)
    # dx = r*cos(lat)*dlon where lat and dlon (i.e. res) are in radians    
    dx = r*cos(lat_2d*deg2rad)*res*deg2rad
    # dy = r*dlat where dlat (i.e. res) is in radians
    # This is constant so reshape to an array of the right dimensions
    dy = zeros(shape(dx)) + r*res*deg2rad

    # Read the land mask from the existing ROMS common grid file
    id = Dataset(common_file, 'r')
    mask_common = id.variables['mask'][:,:]
    id.close()

    # Read FESOM 2D grid
    file = open(mesh_path + 'nod2d.out', 'r')
    file.readline()
    rlon = []
    rlat = []
    for line in file:
        tmp = line.split()
        lon_tmp = float(tmp[1])
        lat_tmp = float(tmp[2])
        if lon_tmp < -180:
            lon_tmp += 360
        elif lon_tmp > 180:
            lon_tmp -= 360
        rlon.append(lon_tmp)
        rlat.append(lat_tmp)
    file.close()
    rlon = array(rlon)
    rlat = array(rlat)
    n2d = size(rlon)
    # Unrotate grid
    lon_fesom, lat_fesom = unrotate_grid(rlon, rlat)

    print('Setting up ' + out_file)
    id = Dataset(out_file, 'w')
    id.createDimension('longitude', size(lon_common))
    id.createDimension('latitude', size(lat_common))
    id.createDimension('time', None)
    id.createVariable('longitude', 'f8', ('longitude'))
    id.variables['longitude'].units = 'degrees'
    id.variables['longitude'][:] = lon_common
    id.createVariable('latitude', 'f8', ('latitude'))
    id.variables['latitude'].units = 'degrees'
    id.variables['latitude'][:] = lat_common
    id.createVariable('time', 'f8', ('time'))
    id.variables['time'].units = 'months'
    id.createVariable('mask', 'f8', ('latitude', 'longitude'))
    id.variables['mask'].units = '1'
    id.variables['mask'][:,:] = mask_common
    id.createVariable('sst', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['sst'].long_name = 'sea surface temperature'
    id.variables['sst'].units = 'C'
    id.createVariable('sss', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['sss'].long_name = 'sea surface salinity'
    id.variables['sss'].units = 'psu'
    id.createVariable('shflux', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['shflux'].long_name = 'surface heat flux into ocean'
    id.variables['shflux'].units = 'W/m^2'
    id.createVariable('ssflux', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['ssflux'].long_name = 'surface virtual salinity flux into ocean'
    id.variables['ssflux'].units = 'psu m/s'
    id.createVariable('aice', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['aice'].long_name = 'sea ice concentration'
    id.variables['aice'].units = '1'
    id.createVariable('hice', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['hice'].long_name = 'sea ice thickness'
    id.variables['hice'].units = 'm'
    id.createVariable('uocn', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['uocn'].long_name = 'ocean surface velocity eastward'
    id.variables['uocn'].units = 'm/s'
    id.createVariable('vocn', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['vocn'].long_name = 'ocean surface velocity northward'
    id.variables['vocn'].units = 'm/s'
    id.createVariable('uice', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['uice'].long_name = 'sea ice velocity eastward'
    id.variables['uice'].units = 'm/s'
    id.createVariable('vice', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['vice'].long_name = 'sea ice velocity northward'
    id.variables['vice'].units = 'm/s'
    id.createVariable('sustr', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['sustr'].long_name = 'zonal surface stress'
    id.variables['sustr'].units = 'N/m^2'
    id.createVariable('svstr', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['svstr'].long_name = 'meridional surface stress'
    id.variables['svstr'].units = 'N/m^2'
    id.createVariable('curl_str', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['curl_str'].long_name = 'curl of surface stress'
    id.variables['curl_str'].units = 'N/m^3'    

    for year in range(start_year, end_year+1):
        print('Processing year ' + str(year))
        for month in range(12):
            print('Processing month ' + str(month+1))
            curr_month = (year-start_year)*12 + month
            # Write time value for this month            
            id.variables['time'][curr_month] = curr_month+1

            # Construct file names
            oce_mean_file = output_dir + expt_name + '.' + str(year) + '.oce.mean.nc'
            forcing_diag_file = output_dir + expt_name + '.' + str(year) + '.forcing.diag.nc'
            ice_mean_file = output_dir + expt_name + '.' + str(year) + '.ice.mean.nc'

            print('...sea surface temperature')
            # Get monthly average of 3D variable
            temp_fesom = monthly_avg(oce_mean_file, 'temp', month)
            # Select surface nodes
            sst_fesom = temp_fesom[:n2d]
            # Interpolate to common grid
            sst_common = interp_fesom2common(lon_common, lat_common, lon_fesom, lat_fesom, sst_fesom)
            # Apply land mask
            sst = ma.masked_where(mask_common==0, sst_common)
            # Write to file
            id.variables['sst'][curr_month,:,:] = sst

            print('...sea surface salinity')
            # Get monthly average of 3D variable
            salt_fesom = monthly_avg(oce_mean_file, 'salt', month)
            # Select surface nodes
            sss_fesom = salt_fesom[:n2d]
            # Interpolate to common grid
            sss_common = interp_fesom2common(lon_common, lat_common, lon_fesom, lat_fesom, sss_fesom)
            # Apply land mask
            sss = ma.masked_where(mask_common==0, sss_common)
            # Write to file
            id.variables['sss'][curr_month,:,:] = sss            

            print('...surface heat flux')
            # Get monthly average
            shflux_fesom = monthly_avg(forcing_diag_file, 'qnet', month)
            # Interpolate to common grid
            shflux_common = interp_fesom2common(lon_common, lat_common, lon_fesom, lat_fesom, shflux_fesom)
            # Apply land mask
            shflux = ma.masked_where(mask_common==0, shflux_common)
            # Write to file
            id.variables['shflux'][curr_month,:,:] = shflux

            print('...surface salt flux')
            # Get monthly average
            ssflux_fesom = monthly_avg(forcing_diag_file, 'virtual_salt', month)
            # Interpolate to common grid
            ssflux_common = interp_fesom2common(lon_common, lat_common, lon_fesom, lat_fesom, ssflux_fesom)
            # Apply land mask
            ssflux = ma.masked_where(mask_common==0, ssflux_common)
            # Write to file
            id.variables['ssflux'][curr_month,:,:] = ssflux

            print('...sea ice concentration')
            # Get monthly average
            aice_fesom = monthly_avg(ice_mean_file, 'area', month)
            # Interpolate to common grid
            aice_common = interp_fesom2common(lon_common, lat_common, lon_fesom, lat_fesom, aice_fesom)
            # Apply land mask
            aice = ma.masked_where(mask_common==0, aice_common)
            # Write to file
            id.variables['aice'][curr_month,:,:] = aice

            print('...sea ice thickness')
            # Get monthly average
            hice_fesom = monthly_avg(ice_mean_file, 'hice', month)
            # Interpolate to common grid
            hice_common = interp_fesom2common(lon_common, lat_common, lon_fesom, lat_fesom, hice_fesom)
            # Apply land mask
            hice = ma.masked_where(mask_common==0, hice_common)
            # Write to file
            id.variables['hice'][curr_month,:,:] = hice

            print('...surface ocean velocity vector')
            # Get monthly averages of both vector components in 3D
            uocn_3d_tmp = monthly_avg(oce_mean_file, 'u', month)
            vocn_3d_tmp = monthly_avg(oce_mean_file, 'v', month)
            # Select surface nodes
            uocn_tmp = uocn_3d_tmp[:n2d]
            vocn_tmp = vocn_3d_tmp[:n2d]
            # Unrotate
            uocn_fesom, vocn_fesom = unrotate_vector(rlon, rlat, uocn_tmp, vocn_tmp)
            # Interpolate to common grid
            uocn_common = interp_fesom2common(lon_common, lat_common, lon_fesom, lat_fesom, uocn_fesom)
            vocn_common = interp_fesom2common(lon_common, lat_common, lon_fesom, lat_fesom, vocn_fesom)
            # Apply land mask
            uocn = ma.masked_where(mask_common==0, uocn_common)
            vocn = ma.masked_where(mask_common==0, vocn_common)
            # Write to file
            id.variables['uocn'][curr_month,:,:] = uocn
            id.variables['vocn'][curr_month,:,:] = vocn

            print('...sea ice velocity vector')
            # Get monthly averages of both vector components
            uice_tmp = monthly_avg(ice_mean_file, 'uice', month)
            vice_tmp = monthly_avg(ice_mean_file, 'vice', month)
            # Unrotate
            uice_fesom, vice_fesom = unrotate_vector(rlon, rlat, uice_tmp, vice_tmp)
            # Interpolate to common grid
            uice_common = interp_fesom2common(lon_common, lat_common, lon_fesom, lat_fesom, uice_fesom)
            vice_common = interp_fesom2common(lon_common, lat_common, lon_fesom, lat_fesom, vice_fesom)
            # Apply land mask
            uice = ma.masked_where(mask_common==0, uice_common)
            vice = ma.masked_where(mask_common==0, vice_common)
            # Write to file
            id.variables['uice'][curr_month,:,:] = uice
            id.variables['vice'][curr_month,:,:] = vice

            print('...surface stress vector')
            # Surface stresses
            # Get monthly averages of both vector components
            sustr_tmp = monthly_avg(forcing_diag_file, 'stress_x', month)
            svstr_tmp = monthly_avg(forcing_diag_file, 'stress_y', month)
            # Unrotate
            sustr_fesom, svstr_fesom = unrotate_vector(rlon, rlat, sustr_tmp, svstr_tmp)
            # Interpolate to common grid
            sustr_common = interp_fesom2common(lon_common, lat_common, lon_fesom, lat_fesom, sustr_fesom)
            svstr_common = interp_fesom2common(lon_common, lat_common, lon_fesom, lat_fesom, svstr_fesom)
            # Apply land mask
            sustr = ma.masked_where(mask_common==0, sustr_common)
            svstr = ma.masked_where(mask_common==0, svstr_common)
            # Write to file
            id.variables['sustr'][curr_month,:,:] = sustr
            id.variables['svstr'][curr_month,:,:] = svstr

            print('...curl of surface stress vector')
            # Curl of surface stress = d/dx (svstr) - d/dy (sustr)
            # First calculate the two derivatives
            dsvstr_dx = ma.empty(shape(svstr_common))
            # Forward difference approximation
            dsvstr_dx[:,:-1] = (svstr_common[:,1:] - svstr_common[:,:-1])/dx[:,:-1]
            # Backward difference for the last row
            dsvstr_dx[:,-1] = (svstr_common[:,-1] - svstr_common[:,-2])/dx[:,-1]
            dsustr_dy = ma.empty(shape(sustr_common))
            dsustr_dy[:-1,:] = (sustr_common[1:,:] - sustr_common[:-1,:])/dy[:-1,:]
            dsustr_dy[-1,:] = (sustr_common[-1,:] - sustr_common[-2,:])/dy[-1,:]
            curl_str = dsvstr_dx - dsustr_dy
            curl_str = ma.masked_where(mask_common==0, curl_str)
            # Write to file
            id.variables['curl_str'][curr_month,:,:] = curl_str

    print('Finished')
    id.close()


# Interpolate the given FESOM field to the regular grid.
# Input:
# lon_1d, lat_1d = 1D arrays containing regular longitude and latitude values
# lon_fesom, lat_fesom = 1D arrays containing (unrotated) longitude and latitude
#                        of each 2D node on the FESOM mesh
# data_fesom = 1D array of data at each 2D node on the FESOM mesh,
#              corresponding to lon_fesom and lat_fesom
# Output:
# data_common = data interpolated to regular grid, dimension lat x lon
def interp_fesom2common (lon_1d, lat_1d, lon_fesom, lat_fesom, data_fesom):

    # Get a 2D field of common latitude and longitude
    lon_2d, lat_2d = meshgrid(lon_1d, lat_1d)

    # Make an array of all the FESOM coordinates (already flattened)
    points = empty([size(lon_fesom), 2])
    points[:,0] = lon_fesom
    points[:,1] = lat_fesom
    # Now make an array of all the common grid coordinates, flattened
    xi = empty([size(lon_2d), 2])
    xi[:,0] = ravel(lon_2d)
    xi[:,1] = ravel(lat_2d)
    # Now call griddata
    result = griddata(points, data_fesom, xi)
    # Un-flatten the result
    data_common = reshape(result, shape(lon_2d))

    return data_common
    

# Command-line interface
if __name__ == "__main__":

    mesh_path = input("Path to FESOM mesh directory: ")
    output_dir = input("Path to FESOM output directory: ")
    start_year = int(input("First year to process: "))
    end_year = int(input("Last year to process: "))
    common_file = input("Path to common-grid file containing land mask as interpolated from ROMS: ")
    out_file = input("Path to desired output file: ")
    common_grid(mesh_path, output_dir, start_year, end_year, common_file, out_file)
    
            
           
        

    
    
  

    
