from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from patches import *
from unrotate_vector import *
from triangle_area import *
from unrotate_grid import *
from in_triangle import *

# var_name = ['hi', 'thdgr', 'sst', 'f/h', 'vel', 'div']
def peninsula_res (var_name):

    # File paths
    mesh_path_lr = '/short/y99/kaa561/FESOM/mesh/low_res/'
    mesh_path_hr = '/short/y99/kaa561/FESOM/mesh/high_res/'
    directory_lr = '/short/y99/kaa561/FESOM/intercomparison_lowres/output/'
    directory_hr = '/short/y99/kaa561/FESOM/intercomparison_highres/output/'
    if var_name == 'hi':
        seasonal_file = 'seasonal_climatology_ice.nc'
    elif var_name in ['thdgr', 'div']:
        seasonal_file = 'seasonal_climatology_ice_diag.nc'
    elif var_name == 'sst':
        seasonal_file = 'seasonal_climatology_oce.nc'
    elif var_name == 'vel':
        seasonal_file = 'seasonal_climatology_uv.nc'
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Bounds on plot (for polar coordinate transformation)
    x_min = -22.5
    x_max = -12
    y_min = 6
    y_max = 15
    if var_name == 'div':
        # Bounds on regular grid
        lon_min = -80
        lon_max = -30
        lat_min = -75
        lat_max = -55
        # Size of regular grid
        num_lon = 500
        num_lat = 500
        # Radius of the Earth in metres
        r = 6.371e6
    if var_name == 'vel':
        num_bins_x = 20
        num_bins_y = 20
    if var_name == 'f/h':
        omega = 7.2921
    # Plotting parameters
    circumpolar = True
    mask_cavities = True
    # Season names for plot titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    # Colour bounds and colour map
    if var_name == 'hi':
        bounds = [0, 1.5]
        cbar_ticks = arange(0, 1.5+0.5, 0.5)
        colour_map = 'jet'
    elif var_name == 'thdgr':
        bounds = [-4, 4]
        cbar_ticks = arange(-4, 4+2, 2)
        colour_map = 'RdBu_r'
    elif var_name == 'div':
        bounds = [-2, 2]
        cbar_ticks = arange(-2, 2+1, 1)
        colour_map = 'RdBu_r'
    elif var_name == 'sst':
        bounds = [-2, 0]
        cbar_ticks = arange(-2, 0+1, 1)
        colour_map = 'jet'
    elif var_name == 'vel':
        bounds = [0, 0.2]
        cbar_ticks = arange(0, 0.2+0.05, 0.05)
        colour_map = 'cool'
    elif var_name == 'f/h':
        bounds = [0.02, 0.06]
        cbar_ticks = arange(0.02, 0.06+0.01, 0.01)
        colour_map = 'jet'        
    if var_name == 'div':
        # Set up regular grid
        # Start with boundaries
        lon_reg_edges = linspace(lon_min, lon_max, num_lon+1)
        lat_reg_edges = linspace(lat_min, lat_max, num_lat+1)
        # Now get centres
        lon_reg = 0.5*(lon_reg_edges[:-1] + lon_reg_edges[1:])
        lat_reg = 0.5*(lat_reg_edges[:-1] + lat_reg_edges[1:])
        # Also get differentials in lon-lat space
        dlon = lon_reg_edges[1:] - lon_reg_edges[:-1]
        dlat = lat_reg_edges[1:] - lat_reg_edges[:-1]
        # Make 2D versions
        lon_reg_2d, lat_reg_2d = meshgrid(lon_reg, lat_reg)
        dlon_2d, dlat_2d = meshgrid(dlon, dlat)
        # Calculate polar coordinates transformation for plotting
        x_reg = -(lat_reg_2d+90)*cos(lon_reg_2d*deg2rad+pi/2)
        y_reg = (lat_reg_2d+90)*sin(lon_reg_2d*deg2rad+pi/2)
        # Calculate differentials in Cartesian space
        dx = r*cos(lat_reg_2d*deg2rad)*dlon_2d*deg2rad
        dy = r*dlat_2d*deg2rad
        # Set up arrays for uhice and vhice at each resolution
        uhice_reg_lr = zeros([4, num_lat, num_lon])
        vhice_reg_lr = zeros([4, num_lat, num_lon])
        uhice_reg_hr = zeros([4, num_lat, num_lon])
        vhice_reg_hr = zeros([4, num_lat, num_lon])
        # Set up colour levels
        lev = linspace(bounds[0], bounds[1], 50)
    if var_name == 'vel':
        # Set up bins for vectors
        x_bins = linspace(x_min, x_max, num=num_bins_x+1)
        y_bins = linspace(y_min, y_max, num=num_bins_y+1)
        # Calculate centres of bins (for plotting)
        x_centres = 0.5*(x_bins[:-1] + x_bins[1:])
        y_centres = 0.5*(y_bins[:-1] + y_bins[1:])

    print 'Processing low-res FESOM'
    # Build mesh
    elements_lr, patches_lr = make_patches(mesh_path_lr, circumpolar, mask_cavities)
    if var_name in ['div', 'vel']:
        # Read rotated latitude and longitude at each node
        file = open(mesh_path_lr + 'nod2d.out', 'r')
        file.readline()
        rlon_lr = []
        rlat_lr = []
        for line in file:
            tmp = line.split()
            lon_tmp = float(tmp[1])
            lat_tmp = float(tmp[2])
            if lon_tmp < -180:
                lon_tmp += 360
            elif lon_tmp > 180:
                lon_tmp -= 360
            rlon_lr.append(lon_tmp)
            rlat_lr.append(lat_tmp)
        file.close()
        rlon_lr = array(rlon_lr)
        rlat_lr = array(rlat_lr)
        if var_name == 'vel':
            # Unrotate
            lon_lr, lat_lr = unrotate_grid(rlon_lr, rlat_lr)
            # Calculate polar coordinates for vector plotting
            x_lr = -(lat_lr+90)*cos(lon_lr*deg2rad+pi/2)
            y_lr = (lat_lr+90)*sin(lon_lr*deg2rad+pi/2)
            # Read cavity flag for each 2D node
            cavity_lr = []
            f = open(mesh_path_lr + 'cavity_flag_nod2d.out', 'r')
            for line in f:
                tmp = int(line)
                if tmp == 1:
                    cavity_lr.append(True)
                elif tmp == 0:
                    cavity_lr.append(False)
                else:
                    print 'Problem'
            f.close()
    # Read data
    if var_name != 'f/h':
        id = Dataset(directory_lr + seasonal_file, 'r')
        if var_name == 'hi':
            data_nodes_lr = id.variables['hice'][:,:]
        elif var_name == 'thdgr':
            data_nodes_lr = id.variables['thdgr'][:,:]*1e7
        elif var_name == 'sst':
            # Fine to read all 3D nodes
            data_nodes_lr = id.variables['temp'][:,:]
        elif var_name == 'vel':
            # Only read the surface nodes
            n2d_lr = size(rlon_lr)
            ur_nodes_lr = id.variables['u'][:,:n2d_lr]
            vr_nodes_lr = id.variables['v'][:,:n2d_lr]
            # Unrotate vectors
            u_nodes_lr, v_nodes_lr = unrotate_vector(rlon_lr, rlat_lr, ur_nodes_lr, vr_nodes_lr)
            # Save speed
            data_nodes_lr = sqrt(u_nodes_lr**2 + v_nodes_lr**2)
        elif var_name == 'div':
            urhice_nodes_lr = id.variables['uhice'][:,:]
            vrhice_nodes_lr = id.variables['vhice'][:,:]
            # Unrotate vectors
            uhice_nodes_lr, vhice_nodes_lr = unrotate_vector(rlon_lr, rlat_lr, urhice_nodes_lr, vrhice_nodes_lr)
        id.close()
    if var_name != 'div':
        # Count the number of elements not in ice shelf cavities
        num_elm_lr = 0
        for elm in elements_lr:
            if not elm.cavity:
                num_elm_lr += 1
        # Set up array for element-averages for each season
        if var_name == 'f/h':
            # No seasonal variation
            data_lr = zeros(num_elm_lr)
        else:
            data_lr = zeros([4, num_elm_lr])
        # Loop over elements to fill this in
        i = 0
        for elm in elements_lr:
            if not elm.cavity:
                # Average over 3 component nodes
                if var_name== 'f/h':
                    data_lr[i] = (abs(2*omega*sin(elm.nodes[0].lat*deg2rad)/elm.nodes[0].find_bottom().depth) + abs(2*omega*sin(elm.nodes[1].lat*deg2rad)/elm.nodes[1].find_bottom().depth) + abs(2*omega*sin(elm.nodes[2].lat*deg2rad)/elm.nodes[2].find_bottom().depth))/3
                else:
                    data_lr[:,i] = (data_nodes_lr[:,elm.nodes[0].id] + data_nodes_lr[:,elm.nodes[1].id] + data_nodes_lr[:,elm.nodes[2].id])/3
                i += 1
    if var_name == 'vel':
        # Make vectors for overlay
        # First set up arrays to integrate velocity in each bin
        # Simple averaging of all the points inside each bin
        ubin_lr = zeros([4, size(y_centres), size(x_centres)])
        vbin_lr = zeros([4, size(y_centres), size(x_centres)])
        num_pts_lr = zeros([4, size(y_centres), size(x_centres)])
        # First convert to polar coordinates, rotate to account for
        # longitude in circumpolar projection, and convert back to vector
        # components
        theta_lr = arctan2(v_nodes_lr, u_nodes_lr)
        theta_circ_lr = theta_lr - lon_lr*deg2rad
        u_circ_lr = data_nodes_lr*cos(theta_circ_lr)
        v_circ_lr = data_nodes_lr*sin(theta_circ_lr)
        # Loop over 2D nodes
        for n in range(n2d_lr):
            # Ignore cavities
            if not cavity_lr[n]:
                if x_lr[n] > x_min and x_lr[n] < x_max and y_lr[n] > y_min and y_lr[n] < y_max:
                    # Figure out which bins this falls into
                    x_index = nonzero(x_bins > x_lr[n])[0][0]-1
                    y_index = nonzero(y_bins > y_lr[n])[0][0]-1
                    # Integrate
                    ubin_lr[:,y_index, x_index] += u_circ_lr[:,n]
                    vbin_lr[:,y_index, x_index] += v_circ_lr[:,n]
                    num_pts_lr[:, y_index, x_index] += 1
        # Convert from sums to averages                    
        # First mask out points with no data
        ubin_lr = ma.masked_where(num_pts_lr==0, ubin_lr)
        vbin_lr = ma.masked_where(num_pts_lr==0, vbin_lr)
        # Divide everything else by number of points
        flag = num_pts_lr > 0
        ubin_lr[flag] = ubin_lr[flag]/num_pts_lr[flag]
        vbin_lr[flag] = vbin_lr[flag]/num_pts_lr[flag]
    if var_name == 'div':
        # Interpolate to regular grid
        # For each element, check if a point on the regular lat-lon grid lies
        # within. If so, do barycentric interpolation to that point.
        for elm in elements_lr:
            # Don't care about ice shelf cavities
            if not elm.cavity:
                # Check if we are within domain of regular grid
                if amax(elm.lon) > lon_min and amin(elm.lon) < lon_max and amax(elm.lat) > lat_min and amin(elm.lat) < lat_max:
                    # Find largest regular longitude value west of Element
                    tmp = nonzero(lon_reg > amin(elm.lon))[0]
                    if len(tmp) == 0:
                        # Element crosses the western boundary
                        iW = 0
                    else:
                        iW = tmp[0] - 1
                    # Find smallest regular longitude value east of Element
                    tmp = nonzero(lon_reg > amax(elm.lon))[0]
                    if len(tmp) == 0:
                        # Element crosses the eastern boundary
                        iE = num_lon
                    else:
                        iE = tmp[0]
                    # Find largest regular latitude value south of Element
                    tmp = nonzero(lat_reg > amin(elm.lat))[0]
                    if len(tmp) == 0:
                        # Element crosses the southern boundary
                        jS = 0
                    else:
                        jS = tmp[0] - 1
                    # Find smallest regular latitude value north of Element
                    tmp = nonzero(lat_reg > amax(elm.lat))[0]
                    if len(tmp) == 0:
                        # Element crosses the northern boundary
                        jN = num_lat
                    else:
                        jN = tmp[0]
                    for i in range(iW+1,iE):
                        for j in range(jS+1,jN):
                            # There is a chance that the regular gridpoint at (i,j)
                            # lies within this element
                            lon0 = lon_reg[i]
                            lat0 = lat_reg[j]
                            if in_triangle(elm, lon0, lat0):
                                # Yes it does
                                # Get area of entire triangle
                                area = triangle_area(elm.lon, elm.lat)
                                # Get area of each sub-triangle formed by
                                # (lon0, lat0)
                                area0 = triangle_area([lon0, elm.lon[1], elm.lon[2]], [lat0, elm.lat[1], elm.lat[2]])
                                area1 = triangle_area([lon0, elm.lon[0], elm.lon[2]], [lat0, elm.lat[0], elm.lat[2]])
                                area2 = triangle_area([lon0, elm.lon[0], elm.lon[1]], [lat0, elm.lat[0], elm.lat[1]])
                                # Find fractional area of each
                                cff = [area0/area, area1/area, area2/area]
                                # Loop over seasons
                                for season in range(4):
                                    # Find value of uhice and vhice at each Node
                                    uhice_vals = []
                                    vhice_vals = []
                                    for n in range(3):
                                        uhice_vals.append(uhice_nodes_lr[season,elm.nodes[n].id])
                                        vhice_vals.append(vhice_nodes_lr[season,elm.nodes[n].id])
                                    # Barycentric interpolation to lon0, lat0
                                    uhice_reg_lr[season,j,i] = sum(array(cff)*array(uhice_vals))
                                    vhice_reg_lr[season,j,i] = sum(array(cff)*array(vhice_vals))
        # Save land mask: wherever identically zero
        mask_lr = ones([num_lat, num_lon])
        index = uhice_reg_lr[0,:,:] == 0.0
        mask_lr[index] = 0.0
        # Calculate divergence
        div_lr = ma.empty(shape(uhice_reg_lr))
        # Loop over seasons
        for season in range(4):
            # First calculate the two derivatives
            duhice_dx_lr = zeros([num_lat, num_lon])
            # Forward difference approximation
            duhice_dx_lr[:,:-1] = (uhice_reg_lr[season,:,1:] - uhice_reg_lr[season,:,:-1])/dx[:,:-1]
            # Backward difference for the last row
            duhice_dx_lr[:,-1] = (uhice_reg_lr[season,:,-1] - uhice_reg_lr[season,:,-2])/dx[:,-1]
            dvhice_dy_lr = zeros([num_lat, num_lon])
            dvhice_dy_lr[:-1,:] = (vhice_reg_lr[season,1:,:] - vhice_reg_lr[season,:-1,:])/dy[:-1,:]
            dvhice_dy_lr[-1,:] = (vhice_reg_lr[season,-1,:] - vhice_reg_lr[season,-2,:])/dy[-1,:]
            # Sum for the divergence
            div_lr_tmp = duhice_dx_lr + dvhice_dy_lr
            # Apply land mask
            div_lr[season,:,:] = ma.masked_where(mask_lr==0, div_lr_tmp)
        # Multiply by 10^7 so colourbar is more readable
        div_lr *= 1e7


    print 'Processing high-res FESOM'
    elements_hr, patches_hr = make_patches(mesh_path_hr, circumpolar, mask_cavities)
    if var_name in ['div', 'vel']:
        file = open(mesh_path_hr + 'nod2d.out', 'r')
        file.readline()
        rlon_hr = []
        rlat_hr = []
        for line in file:
            tmp = line.split()
            lon_tmp = float(tmp[1])
            lat_tmp = float(tmp[2])
            if lon_tmp < -180:
                lon_tmp += 360
            elif lon_tmp > 180:
                lon_tmp -= 360
            rlon_hr.append(lon_tmp)
            rlat_hr.append(lat_tmp)
        file.close()
        rlon_hr = array(rlon_hr)
        rlat_hr = array(rlat_hr)
        if var_name == 'vel':
            lon_hr, lat_hr = unrotate_grid(rlon_hr, rlat_hr)
            x_hr = -(lat_hr+90)*cos(lon_hr*deg2rad+pi/2)
            y_hr = (lat_hr+90)*sin(lon_hr*deg2rad+pi/2)
            cavity_hr = []
            f = open(mesh_path_hr + 'cavity_flag_nod2d.out', 'r')
            for line in f:
                tmp = int(line)
                if tmp == 1:
                    cavity_hr.append(True)
                elif tmp == 0:
                    cavity_hr.append(False)
                else:
                    print 'Problem'
            f.close()
    if var_name != 'f/h':
        id = Dataset(directory_hr + seasonal_file, 'r')
        if var_name == 'hi':
            data_nodes_hr = id.variables['hice'][:,:]
        elif var_name == 'thdgr':
            data_nodes_hr = id.variables['thdgr'][:,:]*1e7
        elif var_name == 'sst':
            data_nodes_hr = id.variables['temp'][:,:]
        elif var_name == 'div':
            urhice_nodes_hr = id.variables['uhice'][:,:]
            vrhice_nodes_hr = id.variables['vhice'][:,:]
            uhice_nodes_hr, vhice_nodes_hr = unrotate_vector(rlon_hr, rlat_hr, urhice_nodes_hr, vrhice_nodes_hr)
        elif var_name == 'vel':
            # Only read the surface nodes
            n2d_hr = size(rlon_hr)
            ur_nodes_hr = id.variables['u'][:,:n2d_hr]
            vr_nodes_hr = id.variables['v'][:,:n2d_hr]
            # Unrotate vectors
            u_nodes_hr, v_nodes_hr = unrotate_vector(rlon_hr, rlat_hr, ur_nodes_hr, vr_nodes_hr)
            data_nodes_hr = sqrt(u_nodes_hr**2 + v_nodes_hr**2)
        id.close()
    if var_name != 'div':
        num_elm_hr = 0
        for elm in elements_hr:
            if not elm.cavity:
                num_elm_hr += 1
        if var_name == 'f/h':
            data_hr = zeros(num_elm_hr)
        else:
            data_hr = zeros([4, num_elm_hr])
        i = 0
        for elm in elements_hr:
            if not elm.cavity:
                if var_name == 'f/h':
                    data_hr[i] = (abs(2*omega*sin(elm.nodes[0].lat*deg2rad)/elm.nodes[0].find_bottom().depth) + abs(2*omega*sin(elm.nodes[1].lat*deg2rad)/elm.nodes[1].find_bottom().depth) + abs(2*omega*sin(elm.nodes[2].lat*deg2rad)/elm.nodes[2].find_bottom().depth))/3
                else:
                    data_hr[:,i] = (data_nodes_hr[:,elm.nodes[0].id] + data_nodes_hr[:,elm.nodes[1].id] + data_nodes_hr[:,elm.nodes[2].id])/3
                i += 1
    if var_name == 'vel':
        ubin_hr = zeros([4, size(y_centres), size(x_centres)])
        vbin_hr = zeros([4, size(y_centres), size(x_centres)])
        num_pts_hr = zeros([4, size(y_centres), size(x_centres)])
        theta_hr = arctan2(v_nodes_hr, u_nodes_hr)
        theta_circ_hr = theta_hr - lon_hr*deg2rad
        u_circ_hr = data_nodes_hr*cos(theta_circ_hr)
        v_circ_hr = data_nodes_hr*sin(theta_circ_hr)
        for n in range(n2d_hr):
            if not cavity_hr[n]:
                if x_hr[n] > x_min and x_hr[n] < x_max and y_hr[n] > y_min and y_hr[n] < y_max:
                    x_index = nonzero(x_bins > x_hr[n])[0][0]-1
                    y_index = nonzero(y_bins > y_hr[n])[0][0]-1
                    ubin_hr[:,y_index, x_index] += u_circ_hr[:,n]
                    vbin_hr[:,y_index, x_index] += v_circ_hr[:,n]
                    num_pts_hr[:, y_index, x_index] += 1
        ubin_hr = ma.masked_where(num_pts_hr==0, ubin_hr)
        vbin_hr = ma.masked_where(num_pts_hr==0, vbin_hr)
        flag = num_pts_hr > 0
        ubin_hr[flag] = ubin_hr[flag]/num_pts_hr[flag]
        vbin_hr[flag] = vbin_hr[flag]/num_pts_hr[flag]
    if var_name == 'div':
        for elm in elements_hr:
            if not elm.cavity:
                if amax(elm.lon) > lon_min and amin(elm.lon) < lon_max and amax(elm.lat) > lat_min and amin(elm.lat) < lat_max:
                    tmp = nonzero(lon_reg > amin(elm.lon))[0]
                    if len(tmp) == 0:
                        iW = 0
                    else:
                        iW = tmp[0] - 1
                    tmp = nonzero(lon_reg > amax(elm.lon))[0]
                    if len(tmp) == 0:
                        iE = num_lon
                    else:
                        iE = tmp[0]
                    tmp = nonzero(lat_reg > amin(elm.lat))[0]
                    if len(tmp) == 0:
                        jS = 0
                    else:
                        jS = tmp[0] - 1
                    tmp = nonzero(lat_reg > amax(elm.lat))[0]
                    if len(tmp) == 0:
                        jN = num_lat
                    else:
                        jN = tmp[0]
                    for i in range(iW+1,iE):
                        for j in range(jS+1,jN):
                            lon0 = lon_reg[i]
                            lat0 = lat_reg[j]
                            if in_triangle(elm, lon0, lat0):
                                area = triangle_area(elm.lon, elm.lat)
                                area0 = triangle_area([lon0, elm.lon[1], elm.lon[2]], [lat0, elm.lat[1], elm.lat[2]])
                                area1 = triangle_area([lon0, elm.lon[0], elm.lon[2]], [lat0, elm.lat[0], elm.lat[2]])
                                area2 = triangle_area([lon0, elm.lon[0], elm.lon[1]], [lat0, elm.lat[0], elm.lat[1]])
                                cff = [area0/area, area1/area, area2/area]
                                for season in range(4):
                                    uhice_vals = []
                                    vhice_vals = []
                                    for n in range(3):
                                        uhice_vals.append(uhice_nodes_hr[season,elm.nodes[n].id])
                                        vhice_vals.append(vhice_nodes_hr[season,elm.nodes[n].id])
                                    uhice_reg_hr[season,j,i] = sum(array(cff)*array(uhice_vals))
                                    vhice_reg_hr[season,j,i] = sum(array(cff)*array(vhice_vals))
        mask_hr = ones([num_lat, num_lon])
        index = uhice_reg_hr[0,:,:] == 0.0
        mask_hr[index] = 0.0
        div_hr = ma.empty(shape(uhice_reg_hr))
        for season in range(4):
            duhice_dx_hr = zeros([num_lat, num_lon])
            duhice_dx_hr[:,:-1] = (uhice_reg_hr[season,:,1:] - uhice_reg_hr[season,:,:-1])/dx[:,:-1]
            duhice_dx_hr[:,-1] = (uhice_reg_hr[season,:,-1] - uhice_reg_hr[season,:,-2])/dx[:,-1]
            dvhice_dy_hr = zeros([num_lat, num_lon])
            dvhice_dy_hr[:-1,:] = (vhice_reg_hr[season,1:,:] - vhice_reg_hr[season,:-1,:])/dy[:-1,:]
            dvhice_dy_hr[-1,:] = (vhice_reg_hr[season,-1,:] - vhice_reg_hr[season,-2,:])/dy[-1,:]
            div_hr_tmp = duhice_dx_hr + dvhice_dy_hr
            div_hr[season,:,:] = ma.masked_where(mask_hr==0, div_hr_tmp)
        div_hr *= 1e7

    print 'Plotting'
    if var_name == 'f/h':
        fig = figure(figsize=(6,9))
        num_t = 1
    else:
        fig = figure(figsize=(19,9))
        num_t = 4
    for season in range(num_t):
        # Low-res
        ax = fig.add_subplot(2, num_t, season+1, aspect='equal')
        if var_name == 'div':
            img = contourf(x_reg, y_reg, div_lr[season,:,:], lev, cmap=colour_map, extend='both')
        else:
            img = PatchCollection(patches_lr, cmap=colour_map)
            if var_name == 'f/h':
                img.set_array(data_lr)
            else:
                img.set_array(data_lr[season,:])
            img.set_clim(vmin=bounds[0], vmax=bounds[1])
            img.set_edgecolor('face')
            ax.add_collection(img)
            if var_name == 'vel':
                # Overlay vectors
                quiver(x_centres, y_centres, ubin_lr[season,:,:], vbin_lr[season,:,:], scale=0.9, headwidth=8, headlength=9, color='black')
        if var_name != 'f/h':
            title(season_names[season], fontsize=24)
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        axis('off')
        if season == 0:
            if var_name == 'f/h':
                text(-24, 14, 'low-res', fontsize=24, ha='right', rotation=90)
            else:
                text(-24, 14, 'low-res', fontsize=24, ha='right')
        # High-res
        ax = fig.add_subplot(2, num_t, season+num_t+1, aspect='equal')
        if var_name == 'div':
            img = contourf(x_reg, y_reg, div_hr[season,:,:], lev, cmap=colour_map, extend='both')
        else:
            img = PatchCollection(patches_hr, cmap=colour_map)
            if var_name == 'f/h':
                img.set_array(data_hr)
            else:
                img.set_array(data_hr[season,:])
            img.set_clim(vmin=bounds[0], vmax=bounds[1])
            img.set_edgecolor('face')
            ax.add_collection(img)
            if var_name == 'vel':
                quiver(x_centres, y_centres, ubin_hr[season,:,:], vbin_hr[season,:,:], scale=0.9, headwidth=8, headlength=9, color='black')
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        axis('off')
        if season == 0:
            if var_name == 'f/h':
                text(-24, 14, 'high-res', fontsize=24, ha='right', rotation=90)
            else:
                text(-24, 14, 'high-res', fontsize=24, ha='right')
    cbaxes = fig.add_axes([0.35, 0.04, 0.3, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, ticks=cbar_ticks, extend='both')
    cbar.ax.tick_params(labelsize=20)
    if var_name == 'hi':
        suptitle('FESOM sea ice effective thickness (m), 1992-2016 average', fontsize=30)
    elif var_name == 'thdgr':
        suptitle(r'FESOM sea ice thermodynamic growth rate (10$^{-7}$ m/s), 1992-2016 average', fontsize=30)
    elif var_name == 'sst':
        suptitle(r'FESOM sea surface temperature ($^{\circ}$C), 1992-2016 average', fontsize=30)
    elif var_name == 'div':
        suptitle(r'FESOM sea ice flux divergence (10$^{-7}$ m/s), 1992-2016 average', fontsize=30)
    elif var_name == 'vel':
        suptitle('FESOM surface ocean velocity (m/s), 1992-2016 average', fontsize=30)
    elif var_name == 'f/h':
        suptitle('abs(f/h)', fontsize=30)
    subplots_adjust(wspace=0.025,hspace=0.025)
    fig.show()
    if var_name == 'f/h':
        fig.savefig('f_h_peninsula_res.png')
    else:
        fig.savefig(var_name + '_peninsula_res.png')


# Command-line interface
if __name__ == "__main__":

    var_name = raw_input("Enter variable to plot (hi, thdgr, sst, f/h, vel, or div): ")
    peninsula_res(var_name)
