from netCDF4 import Dataset
from numpy import *
from fesom_grid import *
from in_triangle import *
from triangle_area import *

def interpolate_nick_climatology (melt_file, temp_file, out_file):

    nick_grid = '/short/y99/kaa561/nick_interpolation/lonlatPISM.nc'
    mesh_path = '/short/y99/kaa561/FESOM/mesh/meshB/'

    # Read Nick's lat and lon points
    id = Dataset(nick_grid, 'r')
    nick_lat = id.variables['lat'][:,:]
    nick_lon = id.variables['lon'][:,:]
    id.close()

    # Set up arrays for interpolated melt rate and surface temp
    melt_reg = ma.empty(shape(nick_lat))
    temp_reg = ma.empty(shape(nick_lat))
    # Fill with NaNs
    melt_reg[:,:] = NaN
    temp_reg[:,:] = NaN

    # Build FESOM mesh
    elements = fesom_grid(mesh_path, circumpolar=True, cross_180=True)
    # Read melt rate and temperature on FESOM mesh
    id = Dataset(melt_file, 'r')
    melt_nodes = mean(id.variables['wnet'][:,:], axis=0)
    id.close()
    id = Dataset(temp_file, 'r')
    temp_nodes = mean(id.variables['temp'][:,:], axis=0)
    id.close()

    # Loop over all cavity elements
    for elm in elements:
        if elm.cavity:
            # Find all grid points which may fall within this triangle
            tmp = where((nick_lat >= amin(elm.lat))*(nick_lat <= amax(elm.lat))*(nick_lon >= amin(elm.lon))*(nick_lon <= amax(elm.lon)))
            j_vals = tmp[0]
            i_vals = tmp[1]
            # Loop over each such grid point
            for point in range(len(j_vals)):
                j = j_vals[point]
                i = i_vals[point]
                lon0 = nick_lon[j,i]
                lat0 = nick_lat[j,i]
                if in_triangle(elm, lon0, lat0):
                    # This point does fall in the triangle
                    # Get area of entire triangle
                    area = triangle_area(elm.lon, elm.lat)
                    # Get area of each sub-triangle formed by (lon0, lat0)
                    area0 = triangle_area([lon0, elm.lon[1], elm.lon[2]], [lat0, elm.lat[1], elm.lat[2]])
                    area1 = triangle_area([lon0, elm.lon[0], elm.lon[2]], [lat0, elm.lat[0], elm.lat[2]])
                    area2 = triangle_area([lon0, elm.lon[0], elm.lon[1]], [lat0, elm.lat[0], elm.lat[1]])
                    # Find fractional area of each
                    cff = [area0/area, area1/area, area2/area]
                    # Find melt rate and temperature at each node
                    melt_vals = []
                    temp_vals = []
                    for n in range(3):
                        melt_vals.append(melt_nodes[elm.nodes[n].id])
                        # This is implicitly surface temp
                        temp_vals.append(temp_nodes[elm.nodes[n].id])
                    # Barycentric interpolation to lon0, lat0
                    melt_reg[j,i] = sum(array(cff)*array(melt_vals))
                    temp_reg[j,i] = sum(array(cff)*array(temp_vals))

    # Mask NaNs
    melt_reg = ma.masked_where(isnan(melt_reg), melt_reg)
    temp_reg = ma.masked_where(isnan(temp_reg), temp_reg)

    # Conversions
    # m/s to mm/s
    melt_reg *= 1e3
    # C to K
    temp_reg += 273.15

    # Output to NetCDF file
    id = Dataset(out_file, 'w')
    id.createDimension('y', size(nick_lat,0))
    id.createDimension('x', size(nick_lat,1))
    id.createDimension('time', None)
    id.createVariable('longitude', 'f8', ('y','x'))
    id.variables['longitude'][:,:] = nick_lon
    id.createVariable('latitude', 'f8', ('y', 'x'))
    id.variables['latitude'][:,:] = nick_lat
    id.createVariable('melt', 'f8', ('y', 'x'))
    id.variables['melt'].units = 'mm/s'
    id.variables['melt'][:,:] = melt_reg
    id.createVariable('temp', 'f8', ('y', 'x'))
    id.variables['temp'].units = 'K'
    id.variables['temp'][:,:] = temp_reg
    id.close()


# Command-line interface
if __name__ == "__main__":

    melt_file = raw_input("Path to FESOM melt rate file: ")
    temp_file = raw_input("Path to FESOM temperature file: ")
    out_file = raw_input("Path to desired output file: ")
    interpolate_nick_climatology(melt_file, temp_file, out_file)
