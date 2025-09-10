from numpy import *
from netCDF4 import Dataset
from fesom_intersectgrid import *

def save_zonalavg (mesh_path, in_file, tstep, out_file):

    # Set bounds on latitude for plot
    lat_min = -90
    lat_max = -30
    # Set bounds on depth
    depth_min = -6000
    depth_max = 0
    # Set number of intervals for latitude and depth to interpolate to
    num_lat = 100
    num_depth = 50

    print('Processing temperature')
    lat_vals, depth_vals, temp = fesom_intersectgrid(mesh_path, in_file, 'temp', tstep, -180, 180, lat_min, lat_max, depth_min, depth_max, num_lat, num_depth)
    # Mask the NaNs
    temp = ma.masked_where(isnan(temp), temp)
    print('Processing salinity')
    lat_vals, depth_vals, salt = fesom_intersectgrid(mesh_path, in_file, 'salt', tstep, -180, 180, lat_min, lat_max, depth_min, depth_max, num_lat, num_depth)
    salt = ma.masked_where(isnan(salt), salt)

    print('Writing ' + out_file)
    id = Dataset(out_file, 'w')
    id.createDimension('latitude', num_lat)
    id.createDimension('depth', num_depth)
    id.createVariable('latitude', 'f8', ('latitude'))
    id.variables['latitude'].units = 'degrees_north'
    id.variables['latitude'][:] = lat_vals
    id.createVariable('depth', 'f8', ('depth'))
    id.variables['depth'].units = 'm'
    id.variables['depth'][:] = depth_vals
    id.createVariable('temp', 'f8', ('depth', 'latitude'))
    id.variables['temp'].units = 'degC'
    id.variables['temp'][:,:] = temp
    id.createVariable('salt', 'f8', ('depth', 'latitude'))
    id.variables['salt'].units = 'psu'
    id.variables['salt'][:,:] = salt
    id.close()


# Command-line interface
if __name__ == "__main__":

    mesh_path = input("Path to FESOM mesh directory: ")
    in_file = input("Path to oce.mean file to zonally average: ")
    tstep = int(input("Time index to zonally average (starting at 1): "))
    out_file = input("Path to desired output file: ")
    save_zonalavg (mesh_path, in_file, tstep, out_file)
    

    
