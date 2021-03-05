### Extract the variables Skip needs for the AAD project, interpolated to a regular grid

import numpy as np
import netCDF4 as nc
import datetime
from fesom_grid import *
from unesco import *
from in_triangle import *
from triangle_area import *


# Read, process, and interpolate a given variable, and save to output file
def process_var (var, output_dir, mesh_path, start_year, end_year, out_file):

    file_head = 'MK44005.'
    [xmin, xmax, ymin, ymax] = [40, 180, -80, -40]
    res = 1/8.
    dt = 5  # days
    days_per_year = 365
    num_time = days_per_year/dt
    sec_per_day = 24*60*60
    mps_to_mpy = days_per_year*sec_per_day
    deg2rad = np.pi/180.
    time_units = 'seconds since 1979-01-01 00:00:00'
    calendar = 'noleap'
    density_anom = 0.03

    # Set parameters for each possible variable
    depth = None
    double_var = False
    mask_outside_cavity = False
    factor = 1
    if var == 'bottom_temp':
        file_tail = '.oce.mean.nc'
        var_in = 'temp'
        depth = 'bottom'
        title = 'Bottom temperature'
        units = 'degC'
    elif var == 'sfc_temp':
        file_tail = '.oce.mean.nc'
        var_in = 'temp'
        depth = 'surface'
        title = 'Surface temperature'
        units = 'degC'
    if var == 'bottom_salt':
        file_tail = '.oce.mean.nc'
        var_in = 'salt'
        depth = 'bottom'
        title = 'Bottom salinity'
        units = 'psu'
    elif var == 'sfc_salt':
        file_tail = '.oce.mean.nc'
        var_in = 'salt'
        depth = 'surface'
        title = 'Surface salinity'
        units = 'psu'
    elif var == 'bottom_speed':
        file_tail = '.oce.mean.nc'
        double_var = True
        var_in_1 = 'u'
        var_in_2 = 'v'
        depth = 'bottom'
        title = 'Bottom current speed'
        units = 'm/s'
    elif var == 'sfc_speed':
        file_tail = '.oce.mean.nc'
        double_var = True
        var_in_1 = 'u'
        var_in_2 = 'v'
        depth = 'surface'
        title = 'Surface current speed'
        units = 'm/s'
    elif var == 'ssh':
        file_tail = '.oce.mean.nc'
        var_in = 'ssh'
        title = 'Sea surface height'
        units = 'm'
    elif var == 'ismr':
        file_tail = '.forcing.diag.nc'
        var_in = 'wnet'
        mask_outside_cavity = True
        factor = mps_to_mpy
        title = 'Ice shelf melt rate'
        units = 'm/y'
    elif var == 'seaice_conc':
        file_tail = '.ice.mean.nc'
        var_in = 'area'
        title = 'Sea ice concentration'
        units = 'fraction'
    elif var == 'seaice_thick':
        file_tail = '.ice.mean.nc'
        var_in = 'hice'
        title = 'Sea ice thickness'
        units = 'm'
    elif var == 'seaice_growth':
        file_tail = '.ice.diag.nc'
        var_in = 'thdgr'
        factor = mps_to_mpy
        title = 'Sea ice thermodynamic growth rate'
        units = 'm/y'
    elif var == 'sfc_stress':
        file_tail = '.forcing.diag.nc'
        double_var = True
        var_in_1 = 'stress_x'
        var_in_2 = 'stress_y'
        title = 'Surface stress'
        units = 'N/m^2'
    elif var == 'mixed_layer_depth':
        file_tail = '.oce.mean.nc'
        double_var = True
        var_in_1 = 'temp'
        var_in_2 = 'salt'
        title = 'Mixed layer depth'
        units = 'm'

    print 'Building FESOM mesh'
    nodes, elements = fesom_grid(mesh_path, return_nodes=True)
    # Count the number of 2D nodes
    f = open(mesh_path+'nod2d.out', 'r')
    n2d = int(f.readline())
    f.close()    
    if mask_outside_cavity:
        # Read the cavity flag
        f = open(mesh_path+'cavity_flag_nod2d.out', 'r')
        cavity = []
        for line in f:
            cavity.append(int(line))
        f.close()
        cavity = np.array(cavity)
    
    print 'Building regular grid'
    lon_reg = np.arange(xmin, xmax+res, res)
    # Iterative latitude axis scaled by cos of latitude, as in mitgcm_python
    lat_reg = [ymin]
    while lat_reg[-1] < ymax+res:
        lat_reg.append(lat_reg[-1] + res*np.cos(lat_reg[-1]*deg2rad))
    lat_reg = np.array(lat_reg)
    for i in range(10):
        lat_reg_c = 0.5*(lat_reg[:-1] + lat_reg[1:])
        j = 0
        lat_reg = [ymin]
        while lat_reg[-1] < ymax+res and j < lat_reg_c.size:
            lat_reg.append(lat_reg[-1] + res*np.cos(lat_reg_c[j]*deg2rad))
            j += 1
        lat_reg = np.array(lat_reg)
    ymax = lat_reg[-1]
    num_lon = lon_reg.size
    num_lat = lat_reg.size

    print 'Setting up '+out_file
    id_out = nc.Dataset(out_file, 'w')
    id_out.createDimension('time', None)
    id_out.createVariable('time', 'f8', ('time'))
    id_out.variables['time'].units = time_units
    id_out.variables['time'].calendar = calendar
    id_out.createDimension('lon', num_lon)
    id_out.createVariable('lon', 'f8', ('lon'))
    id_out.variables['lon'].long_name = 'longitude'
    id_out.variables['lon'].units = 'degrees'
    id_out.variables['lon'][:] = lon_reg
    id_out.createDimension('lat', num_lat)
    id_out.createVariable('lat', 'f8', ('lat'))
    id_out.variables['lat'].long_name = 'latitude'
    id_out.variables['lat'].units = 'degrees'
    id_out.variables['lat'][:] = lat_reg
    id_out.createVariable(var, 'f8', ('time', 'lat', 'lon'))
    id_out.variables[var].long_name = title
    id_out.variables[var].units = units

    t_start = 0
    for year in range(start_year, end_year+1):
        print 'Processing ' + str(year)

        # Set time axis
        time = [nc.date2num(datetime.datetime(year, 1, 1), time_units, calendar=calendar)]
        for t in range(num_time - 1):
            time.append(time[-1] + dt*sec_per_day)
        id_out.variables['time'][t_start:] = np.array(time)
        
        # Read data
        id_in = nc.Dataset(output_dir+file_head+str(year)+file_tail, 'r')
        if double_var:
            # Two variables to read
            data1 = id_in.variables[var_in_1][:]
            data2 = id_in.variables[var_in_2][:]
        else:
            data = id_in.variables[var_in][:]
        id_in.close()
        
        # Process data as needed
        if double_var and var.endswith('speed') or var.endswith('stress'):
            # Get magnitude of vector
            data = np.sqrt(data1**2 + data2**2)
        if var == 'mixed_layer_depth':
            # Calculate density of each node
            density = unesco(data1, data2, np.zeros(data1.shape))
            # Set up array for mixed layer depth
            data = np.zeros([num_time, n2d])
            # Loop over timesteps (I know this is gross)
            for t in range(num_time):
                # Loop over surface nodes
                for i in range(n2d):
                    node = nodes[i]
                    density_sfc = density[t,i]
                    depth_sfc = node.depth
                    depth_tmp = node.depth
                    # Now travel down through the water column until the density exceeds the surface density by density_anom
                    curr_node = node.below
                    while True:
                        if curr_node is None:
                            # Reached the bottom: mixed layer depth is full depth
                            data[t,i] = depth_tmp-depth_sfc
                            break
                        if density[t, curr_node.id] >= density_sfc + density_anom:
                            # Reached the critical density anomaly
                            data[t,i] = curr_node.depth-depth_sfc
                            break
                        depth_tmp = curr_node.depth
                        curr_node = curr_node.below
        if depth == 'surface':
            # Select only the surface nodes
            data = data[:,:n2d]
        elif depth == 'bottom':
            # Select the bottom of each water column
            data_bottom = np.zeros([num_time, n2d])
            for i in range(n2d):
                bottom_id = nodes[i].find_bottom().id
                data_bottom[:,i] = data[:,bottom_id]
            data = data_bottom
        if mask_outside_cavity:
            # Multiply by cavity flag (1 in cavity, 0 outside)
            data *= cavity
        data *= factor

        # Interpolate to regular grid
        data_reg = np.zeros([num_time, num_lat, num_lon])
        valid_mask = np.zeros([num_lat, num_lon])
        # For each element, check if a point on the regular lat-lon grid lies within. If so, do barycentric interpolation to that point.
        for elm in elements:
            # Check if we are within domain of regular grid
            if np.amax(elm.lon) < xmin or np.amin(elm.lon) > xmax or np.amax(elm.lat) < ymin or np.amin(elm.lat) > ymax:
                continue
            # Find largest regular longitude west of element
            tmp = np.nonzero(lon_reg > np.amin(elm.lon))[0]
            if len(tmp) == 0:
                # Element crosses the western boundary
                iW = 0
            else:
                iW = tmp[0] - 1
            # Find smallest regular longitude east of element
            tmp = np.nonzero(lon_reg > np.amax(elm.lon))[0]
            if len(tmp) == 0:
                # Element crosses the eastern boundary
                iE = num_lon
            else:
                iE = tmp[0]
            # Find largest regular latitude south of Element
            tmp = np.nonzero(lat_reg > np.amin(elm.lat))[0]
            if len(tmp) == 0:
                # Element crosses the southern boundary
                jS = 0
            else:
                jS = tmp[0] - 1
            # Find smallest regular latitude value north of Element
            tmp = np.nonzero(lat_reg > np.amax(elm.lat))[0]
            if len(tmp) == 0:
                jN = num_lat
            else:
                jN = tmp[0]
            for i in range(iW+1, iE):
                for j in range(jS+1, jN):
                    # There is a chance that the regular gridpoint at (i,j) lies within this element
                    lon0 = lon_reg[i]
                    lat0 = lat_reg[j]
                    if in_triangle(elm, lon0, lat0):
                        # Get area of entire triangle
                        area = triangle_area(elm.lon, elm.lat)
                        # Get area of each sub-triangle formed by (lon0, lat0)
                        area0 = triangle_area([lon0, elm.lon[1], elm.lon[2]], [lat0, elm.lat[1], elm.lat[2]])
                        area1 = triangle_area([lon0, elm.lon[0], elm.lon[2]], [lat0, elm.lat[0], elm.lat[2]])
                        area2 = triangle_area([lon0, elm.lon[0], elm.lon[1]], [lat0, elm.lat[0], elm.lat[1]])
                        # Find fractional area of each
                        cff = np.array([area0/area, area1/area, area2/area])
                        data_tmp = np.zeros([num_time,3])
                        for n in range(3):
                            data_tmp[:,n] = data[:,elm.nodes[n].id]
                        # Barycentric interpolation to lon0, lat0
                        data_reg[:,j,i] = np.sum(cff[None,:]*data_tmp, axis=1)
                        valid_mask[j,i] = 1
                        # Make sure it doesn't go outside the bounds given by the 3 nodes, at each timestep (truncation errors can cause this)
                        for t in range(num_time):
                            data_reg[t,j,i] = max(data_reg[t,j,i], np.amin(data_tmp[t,:]))
                            data_reg[t,j,i] = min(data_reg[t,j,i], np.amax(data_tmp[t,:]))
        # Mask out anywhere that had nothing to interpolate to
        valid_mask = np.tile(valid_mask, [num_time,1,1])
        data_reg = np.ma.masked_where(valid_mask==0, data_reg)
        # Append to output file
        id_out.variables[var][t_start:,:] = data_reg
        t_start += num_time

    id_out.close()


# Process all variables for the intercomparison high-res simulation
def process_all_intercomparison (base_dir='../FESOM/', out_file_dir='../FESOM/data_for_skip/intercomparison/'):

    output_dir = base_dir+'intercomparison_highres/output/'
    mesh_path = base_dir+'mesh/meshB/'
    start_year = 1997
    end_year = 2016

    for var in ['bottom_temp', 'sfc_temp', 'bottom_salt', 'sfc_salt', 'bottom_speed', 'sfc_speed', 'ssh', 'ismr', 'seaice_conc', 'seaice_thick', 'seaice_growth', 'sfc_stress', 'mixed_layer_depth']:
        print 'Processing variable ' + var
        process_var(var, output_dir, mesh_path, start_year, end_year, out_file_dir+var+'_test.nc')

    
