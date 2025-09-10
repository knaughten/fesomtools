from netCDF4 import Dataset
from numpy import *
from scipy.interpolate import griddata
from unrotate_grid import *

# Make an initial conditions file for the FESOM sea ice, based on NSIDC
# satellite data of monthly-averaged sea ice concentration in January 1992.
# Wherever NSIDC concentration exceeds 15%, set the initial sea ice area to 1,
# ice thickness to 1 m (in the Southern Hemisphere) or 2 m (in the Northern
# Hemisphere), and snow thickness to 0.2 m. Set initial sea ice velocity to 0.
def ice_ini (mesh_path, out_dir):

    # File paths
    nsidc_sh_file = '/short/m68/kaa561/nsidc_aice/seaice_conc_monthly_sh_f11_199201_v02r00.nc'
    nsidc_nh_file = '/short/m68/kaa561/nsidc_aice/seaice_conc_monthly_nh_f11_199201_v02r00.nc'
    out_file = out_dir + 'MK44005.initial_ice.nc'

    # Read rotated latitude and longitude of 2D nodes
    print("Reading FESOM grid")
    f = open(mesh_path + 'nod2d.out', 'r')
    f.readline()
    rlon = []
    rlat = []
    for line in f:
        tmp = line.split()
        lon_tmp = float(tmp[1])
        if lon_tmp < -180:
            lon_tmp += 360
        elif lon_tmp > 180:
            lon_tmp -= 360
        rlon.append(lon_tmp)
        rlat.append(float(tmp[2]))
    f.close()
    rlon = array(rlon)
    rlat = array(rlat)
    # Unrotate grid
    fesom_lon, fesom_lat = unrotate_grid(rlon, rlat)

    # Southern hemisphere
    print("Processing southern hemisphere")
    # Read the NSIDC grid and data
    print("...reading NSIDC")
    id = Dataset(nsidc_sh_file, 'r')
    nsidc_lon = id.variables['longitude'][:,:]
    nsidc_lat = id.variables['latitude'][:,:]
    nsidc_aice = id.variables['seaice_conc_monthly_cdr'][0,:,:]
    id.close()
    # Make sure longitude goes from -180 to 180
    index = nsidc_lon < -180
    nsidc_lon[index] = nsidc_lon[index] + 360
    index = nsidc_lon > 180
    nsidc_lon[index] = nsidc_lon[index] - 360
    # Make sure sea ice concentration values go from 0 to 1
    index = nsidc_aice < 0
    nsidc_aice[index] = 0.0
    index = nsidc_aice > 1
    nsidc_aice[index] = 1.0    
    # Make an nx2 array of the NSIDC coordinates, flattened
    points = empty([size(nsidc_lon), 2])
    points[:,0] = ravel(nsidc_lon)
    points[:,1] = ravel(nsidc_lat)
    # Also flatten the NSIDC data
    values = ravel(nsidc_aice)
    # Now make an mx2 array of the FESOM coordinates
    xi = empty([size(fesom_lon), 2])
    xi[:,0] = fesom_lon
    xi[:,1] = fesom_lat
    # Linear interpolation to the FESOM coordinates
    print("...interpolating to FESOM points")
    area_sh = griddata(points, values, xi, method='linear', fill_value=-999)
    # Also do a nearest-neighbour interpolation
    area_sh2 = griddata(points, values, xi, method='nearest')
    # Replace missing values (from land mask) in linear interpolation with
    # nearest-neighbour interpolation so there are no coastal artifacts
    index = area_sh==-999
    area_sh[index] = area_sh2[index]
    # Cutoff concentration of 15%
    index = area_sh > 0.15
    area_sh[index] = 1.0
    area_sh[~index] = 0.0
    # These cells get 1 m of ice
    hice_sh = zeros(size(area_sh))
    hice_sh[index] = 1.0
    # and 0.2 m of snow
    hsnow_sh = zeros(size(area_sh))
    hsnow_sh[index] = 0.2

    # Northern hemisphere
    print("Processing northern hemisphere")
    # Read the NSIDC grid and data
    print("...reading NSIDC")
    id = Dataset(nsidc_nh_file, 'r')
    nsidc_lon = id.variables['longitude'][:,:]
    nsidc_lat = id.variables['latitude'][:,:]
    nsidc_aice = id.variables['seaice_conc_monthly_cdr'][0,:,:]
    id.close()
    # Make sure longitude goes from -180 to 180
    index = nsidc_lon < -180
    nsidc_lon[index] = nsidc_lon[index] + 360
    index = nsidc_lon > 180
    nsidc_lon[index] = nsidc_lon[index] - 360
    # Make sure sea ice concentration values go from 0 to 1
    index = nsidc_aice < 0
    nsidc_aice[index] = 0.0
    index = nsidc_aice > 1
    nsidc_aice[index] = 1.0
    # Make an nx2 array of the NSIDC coordinates, flattened
    points = empty([size(nsidc_lon), 2])
    points[:,0] = ravel(nsidc_lon)
    points[:,1] = ravel(nsidc_lat)
    # Also flatten the NSIDC data
    values = ravel(nsidc_aice)
    # We already have the FESOM coordinates saved in xi
    # Linear interpolation to the FESOM coordinates
    print("...interpolating to FESOM points")
    area_nh = griddata(points, values, xi, method='linear', fill_value=-999)
    # Also do a nearest-neighbour interpolation
    area_nh2 = griddata(points, values, xi, method='nearest')
    # Replace missing values (from land mask) in linear interpolation with
    # nearest-neighbour interpolation so there are no coastal artifacts
    index = area_nh==-999
    area_nh[index] = area_nh2[index]
    # Cutoff concentration of 15%
    index = area_nh > 0.15
    area_nh[index] = 1.0
    area_nh[~index] = 0.0
    # These cells get 2 m of ice
    hice_nh = zeros(size(area_nh))
    hice_nh[index] = 2.0
    # and 0.2 m of snow
    hsnow_nh = zeros(size(area_nh))
    hsnow_nh[index] = 0.2

    # Add the NH and SH together to get the global arrays
    area = area_sh + area_nh
    hice = hice_sh + hice_nh
    hsnow = hsnow_sh + hsnow_nh
    # Zero ice velocity
    uice = zeros(size(area))
    vice = zeros(size(area))

    # Write to file
    print("Writing output")
    id = Dataset(out_file, 'w')
    id.createDimension('nodes_2d', size(area))
    id.createDimension('T', None)
    id.createVariable('area', 'f8', ('T', 'nodes_2d'))
    id.variables['area'].description = 'ice concentration [0 to 1]'
    id.variables['area'][0,:] = area
    id.createVariable('hice', 'f8', ('T', 'nodes_2d'))
    id.variables['hice'].description = 'effective ice thickness'
    id.variables['hice'].units = 'm'
    id.variables['hice'][0,:] = hice
    id.createVariable('hsnow', 'f8', ('T', 'nodes_2d'))
    id.variables['hsnow'].description = 'effective snow thickness'
    id.variables['hsnow'].units = 'm'
    id.variables['hsnow'][0,:] = hsnow
    id.createVariable('uice', 'f8', ('T', 'nodes_2d'))
    id.variables['uice'].description = 'zonal velocity'
    id.variables['uice'].units = 'm/s'
    id.variables['uice'][0,:] = uice
    id.createVariable('vice', 'f8', ('T', 'nodes_2d'))
    id.variables['vice'].description = 'meridional velocity'
    id.variables['vice'].units = 'm/s'
    id.variables['vice'][0,:] = vice
    id.close()


# Command-line interface
if __name__ == "__main__":

    mesh_path = input("Path to FESOM mesh directory: ")
    out_dir = input("Path to FESOM simulation's output directory: ")
    ice_ini(mesh_path, out_dir)
