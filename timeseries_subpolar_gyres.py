from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from fesom_grid import *
from unrotate_vector import *
from in_triangle import *

def timeseries_subpolar_gyres (mesh_path, output_path, start_year, end_year, log_file, fig_dir=''):

    # Lat-lon bounds on regions to search for each gyre
    # Weddell Sea gyre
    ws_wbdry = -60
    ws_ebdry = 30
    ws_sbdry = -90
    ws_nbdry = -50
    # Ross Sea gyre (crosses 180E)
    rs_wbdry = 150
    rs_ebdry = -140
    rs_sbdry = -90
    rs_nbdry = -60
    # Resolution of regular grid (degrees)
    res = 0.1
    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians coversion factor
    deg2rad = pi/180.0
    # Naming conventions for FESOM output files
    file_head = output_path + 'MK44005.'
    file_tail = '.oce.mean.nc'
    num_years = end_year - start_year + 1

    print 'Building mesh'
    elements = fesom_grid(mesh_path, circumpolar=True, cross_180=True)
    # Read number of nodes, 2D and 3D
    f = open(mesh_path + 'nod2d.out', 'r')
    n2d = int(f.readline())
    f.close()
    f = open(mesh_path + 'nod3d.out', 'r')
    n3d = int(f.readline())
    f.close()
    # Read (rotated) lon, lat, and depth, at each 3D node
    f = open(mesh_path + 'nod3d.out', 'r')
    f.readline()
    rlon = []
    rlat = []
    node_depth = []
    for line in f:
        tmp = line.split()
        lon_tmp = float(tmp[1])
        lat_tmp = float(tmp[2])
        node_depth_tmp = -1*float(tmp[3])
        if lon_tmp < -180:
            lon_tmp += 360
        elif lon_tmp > 180:
            lon_tmp -= 360
        rlon.append(lon_tmp)
        rlat.append(lat_tmp)
        node_depth.append(node_depth_tmp)
    f.close()
    rlon = array(rlon)
    rlat = array(rlat)
    node_depth = array(node_depth)
    # Read lists of which nodes are directly below which
    f = open(mesh_path + 'aux3d.out', 'r')
    max_num_layers = int(f.readline())
    node_columns = zeros([n2d, max_num_layers])
    for n in range(n2d):
        for k in range(max_num_layers):
            node_columns[n,k] = int(f.readline())
    node_columns = node_columns.astype(int)
    f.close()
    # Set up regular grids
    # Weddell Sea gyre
    # Start with boundaries
    ws_lon_reg_edges = arange(ws_wbdry, ws_ebdry+res, res)
    ws_lat_reg_edges = arange(ws_sbdry, ws_nbdry+res, res)
    # Now get centres
    ws_lon_reg = 0.5*(ws_lon_reg_edges[:-1] + ws_lon_reg_edges[1:])
    ws_lat_reg = 0.5*(ws_lat_reg_edges[:-1] + ws_lat_reg_edges[1:])
    # Also get differentials in lon-lat space
    ws_dlon = ws_lon_reg_edges[1:] - ws_lon_reg_edges[:-1]
    ws_dlat = ws_lat_reg_edges[1:] - ws_lat_reg_edges[:-1]
    # Make 2D versions
    ws_lon_reg_2d, ws_lat_reg_2d = meshgrid(ws_lon_reg, ws_lat_reg)
    ws_dlon_2d, ws_dlat_2d = meshgrid(ws_dlon, ws_dlat)
    # Calculate differentials in Cartesian space
    ws_dx = r*cos(ws_lat_reg_2d*deg2rad)*ws_dlon_2d*deg2rad
    ws_dy = r*ws_dlat_2d*deg2rad
    # Ross Sea gyre (split into 2 across 180E)
    rs_lon1_reg_edges = arange(rs_wbdry, 180, res)
    rs_lon2_reg_edges = arange(-180, rs_ebdry+res, res)
    rs_lat_reg_edges = arange(rs_sbdry, rs_nbdry+res, res)
    # Now get centres
    rs_lon1_reg = 0.5*(rs_lon1_reg_edges[:-1] + rs_lon1_reg_edges[1:])
    rs_lon2_reg = 0.5*(rs_lon2_reg_edges[:-1] + rs_lon2_reg_edges[1:])
    rs_lat_reg = 0.5*(rs_lat_reg_edges[:-1] + rs_lat_reg_edges[1:])
    # Also get differentials in lon-lat space
    rs_dlon1 = rs_lon1_reg_edges[1:] - rs_lon1_reg_edges[:-1]
    rs_dlon2 = rs_lon2_reg_edges[1:] - rs_lon2_reg_edges[:-1]
    rs_dlat = rs_lat_reg_edges[1:] - rs_lat_reg_edges[:-1]
    # Make 2D versions
    rs_lon1_reg_2d, rs_lat1_reg_2d = meshgrid(rs_lon1_reg, rs_lat_reg)
    rs_lon2_reg_2d, rs_lat2_reg_2d = meshgrid(rs_lon2_reg, rs_lat_reg)
    rs_dlon1_2d, rs_dlat1_2d = meshgrid(rs_dlon1, rs_dlat)
    rs_dlon2_2d, rs_dlat2_2d = meshgrid(rs_dlon2, rs_dlat)
    # Calculate differentials in Cartesian space
    rs_dx1 = r*cos(rs_lat1_reg_2d*deg2rad)*rs_dlon1_2d*deg2rad
    rs_dx2 = r*cos(rs_lat2_reg_2d*deg2rad)*rs_dlon2_2d*deg2rad
    rs_dy1 = r*rs_dlat1_2d*deg2rad
    rs_dy2 = r*rs_dlat2_2d*deg2rad

    print 'Reading data'
    u = empty([num_years, n3d])
    for year in range(start_year, end_year+1):
        print '...' + str(year)
        # Read horizontal velocity components for this year, annually average
        id = Dataset(file_head + str(year) + file_tail, 'r')
        ur = mean(id.variables['u'][:,:], axis=0)
        vr = mean(id.variables['v'][:,:], axis=0)
        id.close()
        # Unrotate
        u_tmp, v_tmp = unrotate_vector(rlon, rlat, ur, vr)
        # Save in array
        u[year-start_year,:] = u_tmp

    print 'Vertically integrating u*dz'
    int_udz = zeros([num_years, n2d])
    # Loop over nodes
    for n in range(n2d):
        # Loop over depth
        for k in range(max_num_layers-1):
            if node_columns[n,k+1] == -999:
                # Reached the bottom
                break
            # Trapezoidal rule
            top_id = node_columns[n,k]
            bot_id = node_columns[n,k+1]
            dz = node_depth[bot_id-1] - node_depth[top_id-1]
            # Now loop over years
            for year in range(num_years):
                int_udz[year,n] += 0.5*(u[year,top_id-1] + u[year,bot_id-1])*dz

    print 'Interpolating to regular grid'
    int_udz_reg_ws = zeros([num_years, size(ws_dy,0), size(ws_dy,1)])
    int_udz_reg_rs1 = zeros([num_years, size(rs_dy1,0), size(rs_dy1,1)])
    int_udz_reg_rs2 = zeros([num_years, size(rs_dy2,0), size(rs_dy2,1)])
    # For each element, check if a point on the regular lat-lon grid lies
    # within. If so, do barycentric interpolation to that point.
    for elm in elements:
        # Weddell Sea
        if amax(elm.lon) > ws_wbdry and amin(elm.lon) < ws_ebdry and amin(elm.lat) < ws_nbdry:
            # Find largest regular longitude value west of Element
            tmp = nonzero(ws_lon_reg > amin(elm.lon))[0]
            if len(tmp) == 0:
                # Element crosses the western boundary
                iW = 0
            else:
                iW = tmp[0] - 1
            # Find smallest regular longitude value east of Element
            tmp = nonzero(ws_lon_reg > amax(elm.lon))[0]
            if len(tmp) == 0:
                # Element crosses the eastern boundary
                iE = size(ws_lon_reg)
            else:
                iE = tmp[0]
            # Find largest regular latitude value south of Element
            tmp = nonzero(ws_lat_reg > amin(elm.lat))[0]
            if len(tmp) == 0:
                # Element crosses the southern boundary
                jS = 0
            else:
                jS = tmp[0] - 1
            # Find smallest regular latitude value north of Element
            tmp = nonzero(ws_lat_reg > amax(elm.lat))[0]
            if len(tmp) == 0:
                # Element crosses the northern boundary
                jN = size(ws_lat_reg)
            else:
                jN = tmp[0]
            for i in range(iW+1,iE):
                for j in range(jS+1,jN):
                    # There is a chance that the regular gridpoint at (i,j)
                    # lies within this element
                    lon0 = ws_lon_reg[i]
                    lat0 = ws_lat_reg[j]
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
                        # Now loop over years
                        for year in range(num_years):
                            # Find value of int_udz at each Node
                            vals = []
                            for n in range(3):
                                vals.append(int_udz[year,elm.nodes[n].id])
                            # Barycentric interpolation to lon0, lat0
                            int_udz_reg_ws[year,j,i] = sum(array(cff)*array(vals))
        # Ross Sea, part 1
        if amax(elm.lon) > rs_wbdry and amin(elm.lat) < rs_nbdry:
            tmp = nonzero(rs_lon1_reg > amin(elm.lon))[0]
            if len(tmp) == 0:
                iW = 0
            else:
                iW = tmp[0] - 1
            tmp = nonzero(rs_lon1_reg > amax(elm.lon))[0]
            if len(tmp) == 0:
                iE = size(rs_lon1_reg)
            else:
                iE = tmp[0]
            tmp = nonzero(rs_lat_reg > amin(elm.lat))[0]
            if len(tmp) == 0:
                jS = 0
            else:
                jS = tmp[0] - 1
            tmp = nonzero(rs_lat_reg > amax(elm.lat))[0]
            if len(tmp) == 0:
                jN = size(rs_lat_reg)
            else:
                jN = tmp[0]
            for i in range(iW+1,iE):
                for j in range(jS+1,jN):
                    lon0 = rs_lon1_reg[i]
                    lat0 = rs_lat_reg[j]
                    if in_triangle(elm, lon0, lat0):
                        area = triangle_area(elm.lon, elm.lat)
                        area0 = triangle_area([lon0, elm.lon[1], elm.lon[2]], [lat0, elm.lat[1], elm.lat[2]])
                        area1 = triangle_area([lon0, elm.lon[0], elm.lon[2]], [lat0, elm.lat[0], elm.lat[2]])
                        area2 = triangle_area([lon0, elm.lon[0], elm.lon[1]], [lat0, elm.lat[0], elm.lat[1]])
                        cff = [area0/area, area1/area, area2/area]
                        for year in range(num_years):
                            vals = []
                            for n in range(3):
                                vals.append(int_udz[year,elm.nodes[n].id])
                            int_udz_reg_rs1[year,j,i] = sum(array(cff)*array(vals))
        # Ross Sea, part 2
        if amin(elm.lon) < rs_ebdry and amin(elm.lat) < rs_nbdry:
            tmp = nonzero(rs_lon2_reg > amin(elm.lon))[0]
            if len(tmp) == 0:
                iW = 0
            else:
                iW = tmp[0] - 1
            tmp = nonzero(rs_lon2_reg > amax(elm.lon))[0]
            if len(tmp) == 0:
                iE = size(rs_lon2_reg)
            else:
                iE = tmp[0]
            tmp = nonzero(rs_lat_reg > amin(elm.lat))[0]
            if len(tmp) == 0:
                jS = 0
            else:
                jS = tmp[0] - 1
            tmp = nonzero(rs_lat_reg > amax(elm.lat))[0]
            if len(tmp) == 0:
                jN = size(rs_lat_reg)
            else:
                jN = tmp[0]
            for i in range(iW+1,iE):
                for j in range(jS+1,jN):
                    lon0 = rs_lon2_reg[i]
                    lat0 = rs_lat_reg[j]
                    if in_triangle(elm, lon0, lat0):
                        area = triangle_area(elm.lon, elm.lat)
                        area0 = triangle_area([lon0, elm.lon[1], elm.lon[2]], [lat0, elm.lat[1], elm.lat[2]])
                        area1 = triangle_area([lon0, elm.lon[0], elm.lon[2]], [lat0, elm.lat[0], elm.lat[2]])
                        area2 = triangle_area([lon0, elm.lon[0], elm.lon[1]], [lat0, elm.lat[0], elm.lat[1]])
                        cff = [area0/area, area1/area, area2/area]
                        for year in range(num_years):
                            vals = []
                            for n in range(3):
                                vals.append(int_udz[year,elm.nodes[n].id])
                            int_udz_reg_rs2[year,j,i] = sum(array(cff)*array(vals))

    # Indefinite integral from south to north of udz*dy, convert to Sv
    strf_ws = cumsum(int_udz_reg_ws*ws_dy, axis=1)*1e-6
    strf_rs1 = cumsum(int_udz_reg_rs1*rs_dy1, axis=1)*1e-6
    strf_rs2 = cumsum(int_udz_reg_rs2*rs_dy2, axis=1)*1e-6
    # Build timeseries
    ws_trans = empty(num_years)
    rs_trans = empty(num_years)
    for year in range(num_years):
        # Find most negative value    
        ws_trans[year] = -1*amin(strf_ws[year,:])
        rs_trans[year] = -1*min(amin(strf_rs1[year,:]), amin(strf_rs2[year,:]))

    # Make time axis
    time = range(start_year, end_year+1)

    print 'Plotting'
    # Weddell Sea
    fig = figure()
    plot(time, ws_trans)
    xlabel('year')
    ylabel('Sv')
    xlim([start_year, end_year])
    title('Weddell Sea Gyre transport')
    grid(True)
    fig.savefig(fig_dir + 'weddell_gyre.png')
    # Ross Sea
    fig = figure()
    plot(time, rs_trans)
    xlabel('year')
    ylabel('Sv')
    xlim([start_year, end_year])
    title('Ross Sea Gyre transport')
    grid(True)
    fig.savefig(fig_dir + 'ross_gyre.png')

    print 'Saving results to log file'
    f = open(log_file, 'w')
    f.write('Weddell Sea Gyre transport (Sv)\n')
    for t in range(num_years):
        f.write(str(ws_trans[t]) + '\n')
    f.write('Ross Sea Gyre transport (Sv)\n')
    for t in range(num_years):
        f.write(str(rs_trans[t]) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    output_path = raw_input("Path to FESOM output directory: ")
    start_year = int(raw_input("First year to process: "))
    end_year = int(raw_input("Last year to process: "))
    log_file = raw_input("Desired path to logfile: ")
    fig_dir = raw_input("Path to directories to save figures: ")
    timeseries_subpolar_gyres(mesh_path, output_path, start_year, end_year, log_file, fig_dir)
    
            
            

    
    
    
