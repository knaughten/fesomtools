from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from fesom_grid import *
from unrotate_vector import *
from in_triangle import *

def wind_stress_curl ():

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/meshB/'
    directory_beg =  '/short/y99/kaa561/FESOM/highres_spinup/'
    directories = ['/short/y99/kaa561/FESOM/rcp45_M/', '/short/y99/kaa561/FESOM/rcp45_A/', '/short/y99/kaa561/FESOM/rcp85_M/', '/short/y99/kaa561/FESOM/rcp85_A/']
    file_beg = 'annual_avg.forcing.diag.1996.2005.nc'
    file_end = 'annual_avg.forcing.diag.2091.2100.nc'
    # Titles for plotting
    expt_names = ['RCP 4.5 MMM', 'RCP 4.5 ACCESS', 'RCP 8.5 MMM', 'RCP 8.5 ACCESS']
    expt_filenames = ['rcp45_m', 'rcp45_a', 'rcp85_m', 'rcp85_a']
    num_expts = len(directories)
    colours = ['blue', 'cyan', 'green', 'magenta']
    # Bounds on regular grid
    lon_min = -180
    lon_max = 180
    lat_min = -75
    lat_max = -50
    # Number of points on regular grid
    num_lon = 1000
    num_lat = 200
    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians coversion factor
    deg2rad = pi/180.0
    # Don't consider values above this threshold (small, negative)
    threshold = -5e-8

    print 'Building mesh'
    elements = fesom_grid(mesh_path, circumpolar=True, cross_180=True)
    # Read (rotated) lon and lat at each 2D node
    f = open(mesh_path + 'nod2d.out', 'r')
    n2d = int(f.readline())
    rlon = []
    rlat = []
    for line in f:
        tmp = line.split()
        lon_tmp = float(tmp[1])
        lat_tmp = float(tmp[2])
        if lon_tmp < -180:
            lon_tmp += 360
        elif lon_tmp > 180:
            lon_tmp -= 360
        rlon.append(lon_tmp)
        rlat.append(lat_tmp)
    f.close()
    rlon = array(rlon)
    rlat = array(rlat)

    print 'Reading data'
    print '...1996-2005'
    # Read rotated wind stress components
    id = Dataset(directory_beg + file_beg, 'r')
    stress_xr = id.variables['stress_x'][0,:]
    stress_yr = id.variables['stress_y'][0,:]
    id.close()
    # Unrotate
    stress_x_beg, stress_y_beg = unrotate_vector(rlon, rlat, stress_xr, stress_yr)
    # Set up array for wind stress in each RCP experiment
    stress_x_end = zeros([num_expts, n2d])
    stress_y_end = zeros([num_expts, n2d])
    for expt in range(num_expts):
        print '...' + expt_names[expt]
        id = Dataset(directories[expt] + file_end, 'r')
        stress_xr = id.variables['stress_x'][0,:]
        stress_yr = id.variables['stress_y'][0,:]
        id.close()
        stress_x_tmp, stress_y_tmp = unrotate_vector(rlon, rlat, stress_xr, stress_yr)
        stress_x_end[expt,:] = stress_x_tmp
        stress_y_end[expt,:] = stress_y_tmp

    print 'Interpolating to regular grid'
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
    # Calculate differentials in Cartesian space
    dx = r*cos(lat_reg_2d*deg2rad)*dlon_2d*deg2rad
    dy = r*dlat_2d*deg2rad
    # Set up arrays for result
    stress_x_reg_beg = zeros([num_lat, num_lon])
    stress_y_reg_beg = zeros([num_lat, num_lon])
    stress_x_reg_end = zeros([num_expts, num_lat, num_lon])
    stress_y_reg_end = zeros([num_expts, num_lat, num_lon])
    # For each element, check if a point on the regular lat-lon grid lies
    # within. If so, do barycentric interpolation to that point.
    for elm in elements:
        # Check if we are within domain of regular grid
        if amin(elm.lat) > lat_max:
            continue
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
                    # Get area of entire triangle
                    area = triangle_area(elm.lon, elm.lat)
                    # Get area of each sub-triangle formed by
                    # (lon0, lat0)
                    area0 = triangle_area([lon0, elm.lon[1], elm.lon[2]], [lat0, elm.lat[1], elm.lat[2]])
                    area1 = triangle_area([lon0, elm.lon[0], elm.lon[2]], [lat0, elm.lat[0], elm.lat[2]])
                    area2 = triangle_area([lon0, elm.lon[0], elm.lon[1]], [lat0, elm.lat[0], elm.lat[1]])
                    # Find fractional area of each
                    cff = [area0/area, area1/area, area2/area]
                    # 1996-2005
                    # Find value of stress_x and stress_y at each Node
                    vals_x = []
                    vals_y = []
                    for n in range(3):
                        vals_x.append(stress_x_beg[elm.nodes[n].id])
                        vals_y.append(stress_y_beg[elm.nodes[n].id])
                    # Barycentric interpolation to lon0, lat0
                    stress_x_reg_beg[j,i] = sum(array(cff)*array(vals_x))
                    stress_y_reg_beg[j,i] = sum(array(cff)*array(vals_y))
                    # RCPs
                    for expt in range(num_expts):
                        vals_x = []
                        vals_y = []
                        for n in range(3):
                            vals_x.append(stress_x_end[expt,elm.nodes[n].id])
                            vals_y.append(stress_y_end[expt,elm.nodes[n].id])
                        stress_x_reg_end[expt,j,i] = sum(array(cff)*array(vals_x))
                        stress_y_reg_end[expt,j,i] = sum(array(cff)*array(vals_y))

    print 'Calculating curl'
    # 1996-2005
    # First calculate the two derivatives
    dv_dx = zeros(shape(stress_x_reg_beg))
    du_dy = zeros(shape(stress_x_reg_beg))
    # Forward difference approximation
    dv_dx[:,:-1] = (stress_y_reg_beg[:,1:] - stress_y_reg_beg[:,:-1])/dx[:,:-1]
    du_dy[:-1,:] = (stress_x_reg_beg[1:,:] - stress_x_reg_beg[:-1,:])/dy[:-1,:]
    # Backward difference for the last row
    dv_dx[:,-1] = (stress_y_reg_beg[:,-1] - stress_y_reg_beg[:,-2])/dx[:,-1]
    du_dy[-1,:] = (stress_x_reg_beg[-1,:] - stress_x_reg_beg[-2,:])/dy[-1,:]
    curl_beg = dv_dx - du_dy
    # RCPs
    curl_end_tmp = zeros(shape(stress_x_reg_end))
    for expt in range(num_expts):
        dv_dx = zeros(shape(stress_x_reg_beg))
        du_dy = zeros(shape(stress_x_reg_beg))
        dv_dx[:,:-1] = (stress_y_reg_end[expt,:,1:] - stress_y_reg_end[expt,:,:-1])/dx[:,:-1]
        du_dy[:-1,:] = (stress_x_reg_end[expt,1:,:] - stress_x_reg_end[expt,:-1,:])/dy[:-1,:]
        dv_dx[:,-1] = (stress_y_reg_end[expt,:,-1] - stress_y_reg_end[expt,:,-2])/dx[:,-1]
        du_dy[-1,:] = (stress_x_reg_end[expt,-1,:] - stress_x_reg_end[expt,-2,:])/dy[-1,:]
        curl_end_tmp[expt,:,:] = dv_dx - du_dy

    print 'Plotting zonal averages'
    # Calculate zonal averages
    curl_beg_avg = mean(curl_beg, axis=1)
    curl_end_avg = mean(curl_end_tmp, axis=2)
    # Plot zonal averages
    fig, ax = subplots(figsize=(10,6))
    ax.plot(curl_beg_avg, lat_reg, label='1996-2005', color='black', linewidth=2)
    for expt in range(num_expts):
        ax.plot(curl_end_avg[expt,:], lat_reg, label=expt_names[expt], color=colours[expt], linewidth=2)
    title('Curl of wind stress (2091-2100)', fontsize=18)
    xlabel(r'N/m$^3$', fontsize=14)
    ylabel('latitude', fontsize=14)
    ylim([lat_min, lat_max])
    grid(True)
    # Move plot over to make room for legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    # Make legend
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.show()
    fig.savefig('windstress_curl_rcp.png')
    # Plot anomalies in zonal averages
    fig, ax = subplots(figsize=(10,6))
    for expt in range(num_expts):
        ax.plot(curl_end_avg[expt,:]-curl_beg_avg, lat_reg, label=expt_names[expt], color=colours[expt], linewidth=2)
    title('Anomalies in curl of wind stress (2091-2100 minus 1996-2005)', fontsize=18)
    xlabel(r'N/m$^3$', fontsize=14)
    ylabel('latitude', fontsize=14)
    ylim([lat_min, lat_max])
    grid(True)
    # Move plot over to make room for legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    # Make legend
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.show()
    fig.savefig('windstress_curl_diff_rcp.png')
    # Plot percent change in zonal averages
    fig, ax = subplots(figsize=(10,6))
    for expt in range(num_expts):
        ax.plot((curl_end_avg[expt,:]-curl_beg_avg)/curl_beg_avg*100, lat_reg, label=expt_names[expt], color=colours[expt], linewidth=2)
    title('Percent change in curl of wind stress (2091-2100 minus 1996-2005)', fontsize=18)
    xlabel('%', fontsize=14)
    ylabel('latitude', fontsize=14)
    xlim([-20, 20])
    ylim([-65.5, -58])
    grid(True)
    # Move plot over to make room for legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    # Make legend
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.show()
    fig.savefig('windstress_curl_percent_rcp.png')

    print 'Plotting 2D fields'
    # First mask out regions above threshold at beginning
    curl_beg = ma.masked_where(curl_beg > threshold, curl_beg)
    curl_end = ma.empty(shape(curl_end_tmp))
    for expt in range(num_expts):
        curl_end[expt,:,:] = ma.masked_where(curl_beg > threshold, curl_end_tmp[expt,:,:])
    # Calculate percent change for each RCP
    percent_change = ma.empty(shape(curl_end))
    for expt in range(num_expts):
        percent_change[expt,:,:] = (curl_end[expt,:,:] - curl_beg)/curl_beg*100
    # 1996-2005
    fig, ax = subplots(figsize=(10,6))
    bound = 1e-6
    lev = linspace(-bound, bound, num=50)
    img = ax.contourf(lon_reg, lat_reg, curl_beg, lev, cmap='RdBu_r')
    xlabel('Longitude')
    ylabel('Latitude')
    xlim([lon_min, lon_max])
    ylim([lat_min, lat_max])
    title('Wind stress curl, 1996-2005 (N/m^3)', fontsize=18)
    colorbar(img)
    fig.show()
    fig.savefig('windstress_curl_2D_beg.png')
    for expt in range(num_expts):
        bound = 50
        lev = linspace(-bound, bound, num=50)
        fig, ax = subplots(figsize=(10,6))
        img = ax.contourf(lon_reg, lat_reg, percent_change[expt,:,:], lev, cmap='RdBu_r')
        xlabel('Longitude')
        ylabel('Latitude')
        xlim([lon_min, lon_max])
        ylim([lat_min, lat_max])
        title('Percent change in wind stress curl, 2091-2100 versus 1996-2005 (' + expt_names[expt] + ')', fontsize=18)
        colorbar(img)
        fig.show()
        fig.savefig('windstress_curl_2D_percent_'+expt_filenames[expt]+'.png')


# Command-line interface
if __name__ == "__main__":

    wind_stress_curl()
