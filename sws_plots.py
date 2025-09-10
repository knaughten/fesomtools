from netCDF4 import Dataset
from numpy import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.pyplot import *
from fesom_grid import *
from unrotate_grid import *
from unrotate_vector import *

def sws_plots ():

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/meshB/'
    oce_file_beg = '/short/y99/kaa561/FESOM/highres_spinup/annual_avg.oce.mean.1996.2005.nc'
    oce_file_end = '/short/y99/kaa561/FESOM/rcp85_A/annual_avg.oce.mean.2091.2100.nc'
    ice_file_beg = '/short/y99/kaa561/FESOM/highres_spinup/seasonal_climatology_ice_diag_1996_2005.nc'
    ice_file_end = '/short/y99/kaa561/FESOM/rcp85_A/seasonal_climatology_ice_diag_2091_2100.nc'
    # Bounds on plot (in polar coordinate transformation)
    x_min = -22
    x_max = -4
    y_min = 1
    y_max = 15
    # Number of bins for velocity vectors
    num_bins_x = 20
    num_bins_y = 20
    # Plotting parameters
    circumpolar = True
    # Season names for plot titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    # Degrees to radians conversion factor
    deg2rad = pi/180.0

    print('Building mesh')
    elements = fesom_grid(mesh_path, circumpolar)
    # Build one set of plotting patches with all elements, and one with
    # ice shelf cavities masked
    patches_all = []
    patches_ocn = []
    for elm in elements:
        coord = transpose(vstack((elm.x, elm.y)))
        patches_all.append(Polygon(coord, True, linewidth=0.))
        if not elm.cavity:
            patches_ocn.append(Polygon(coord, True, linewidth=0.))
    num_elm = len(patches_all)
    num_elm_ocn = len(patches_ocn)
    # Build ice shelf front contours
    contour_lines = []
    for elm in elements:
        # Select elements where exactly 2 of the 3 nodes are in a cavity
        if count_nonzero(elm.cavity_nodes) == 2:
            # Save the coastal flags and x- and y- coordinates of these 2
            coast_tmp = []
            x_tmp = []
            y_tmp = []
            for i in range(3):
                if elm.cavity_nodes[i]:
                    coast_tmp.append(elm.coast_nodes[i])
                    x_tmp.append(elm.x[i])
                    y_tmp.append(elm.y[i])
            # Select elements where at most 1 of these 2 nodes are coastal
            if count_nonzero(coast_tmp) < 2:
                # Draw a line between the 2 nodes
                contour_lines.append([(x_tmp[0], y_tmp[0]), (x_tmp[1], y_tmp[1])])

    print('Processing bathymetry')
    # Calculate bathymetry (depth of bottom node) averaged over 3 nodes making
    # up each element
    bathy = []
    for elm in elements:
        bathy.append(mean([elm.nodes[0].find_bottom().depth, elm.nodes[1].find_bottom().depth, elm.nodes[2].find_bottom().depth]))
    # Plot
    fig = figure(figsize=(12,8))
    ax = fig.add_subplot(1,1,1, aspect='equal')
    img = PatchCollection(patches_all, cmap='jet')
    img.set_array(array(bathy))
    img.set_edgecolor('face')
    img.set_clim(vmin=0, vmax=1500)
    ax.add_collection(img)
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax.add_collection(contours)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    title('Bathymetry (m)', fontsize=24)
    cbar = colorbar(img, extend='max')
    fig.show()
    fig.savefig('sws_bathy.png')

    print('Processing bottom water temperature')
    # Read annually averaged data
    id = Dataset(oce_file_beg, 'r')
    temp_nodes_beg = id.variables['temp'][0,:]
    id.close()
    id = Dataset(oce_file_end, 'r')
    temp_nodes_end = id.variables['temp'][0,:]
    id.close()
    # Now average bottom node temperatures over each element
    bwtemp_beg = []
    bwtemp_end = []
    for elm in elements:
        bwtemp_beg.append(mean([temp_nodes_beg[elm.nodes[0].find_bottom().id], temp_nodes_beg[elm.nodes[1].find_bottom().id], temp_nodes_beg[elm.nodes[2].find_bottom().id]]))
        bwtemp_end.append(mean([temp_nodes_end[elm.nodes[0].find_bottom().id], temp_nodes_end[elm.nodes[1].find_bottom().id], temp_nodes_end[elm.nodes[2].find_bottom().id]]))
    # Plot
    fig = figure(figsize=(16,7))
    # 1996-2005
    ax = fig.add_subplot(1, 2, 1, aspect='equal')
    img = PatchCollection(patches_all, cmap='jet')
    img.set_array(array(bwtemp_beg))
    img.set_edgecolor('face')
    img.set_clim(vmin=-2, vmax=-0.5)
    ax.add_collection(img)
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax.add_collection(contours)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    title('1996-2005', fontsize=20)
    # 2091-2100
    ax = fig.add_subplot(1, 2, 2, aspect='equal')
    img = PatchCollection(patches_all, cmap='jet')
    img.set_array(array(bwtemp_end))
    img.set_edgecolor('face')
    img.set_clim(vmin=-2, vmax=-0.5)
    ax.add_collection(img)
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax.add_collection(contours)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    title('2091-2100', fontsize=20)
    cbaxes = fig.add_axes([0.35, 0.04, 0.3, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='both')
    suptitle(r'Bottom water temperature ($^{\circ}$C)', fontsize=24)
    subplots_adjust(wspace=0.025, hspace=0.025)
    fig.show()
    fig.savefig('sws_bwtemp.png')

    print('Processing bottom water salinity')
    # Read annually averaged data
    id = Dataset(oce_file_beg, 'r')
    salt_nodes_beg = id.variables['salt'][0,:]
    id.close()
    id = Dataset(oce_file_end, 'r')
    salt_nodes_end = id.variables['salt'][0,:]
    id.close()
    # Now average bottom node salinities over each element
    bwsalt_beg = []
    bwsalt_end = []
    for elm in elements:
        bwsalt_beg.append(mean([salt_nodes_beg[elm.nodes[0].find_bottom().id], salt_nodes_beg[elm.nodes[1].find_bottom().id], salt_nodes_beg[elm.nodes[2].find_bottom().id]]))
        bwsalt_end.append(mean([salt_nodes_end[elm.nodes[0].find_bottom().id], salt_nodes_end[elm.nodes[1].find_bottom().id], salt_nodes_end[elm.nodes[2].find_bottom().id]]))
    # Plot beginning 
    fig = figure(figsize=(16,7))
    # 1996-2005
    ax = fig.add_subplot(1, 2, 1, aspect='equal')
    img = PatchCollection(patches_all, cmap='jet')
    img.set_array(array(bwsalt_beg))
    img.set_edgecolor('face')
    img.set_clim(vmin=34.2, vmax=34.7)
    ax.add_collection(img)
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax.add_collection(contours)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    title('1996-2005', fontsize=20)
    # 2091-2100
    ax = fig.add_subplot(1, 2, 2, aspect='equal')
    img = PatchCollection(patches_all, cmap='jet')
    img.set_array(array(bwsalt_end))
    img.set_edgecolor('face')
    img.set_clim(vmin=34.2, vmax=34.7)
    ax.add_collection(img)
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax.add_collection(contours)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    title('2091-2100', fontsize=20)
    cbaxes = fig.add_axes([0.35, 0.04, 0.3, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='both')
    suptitle('Bottom water salinity (psu)', fontsize=24)
    subplots_adjust(wspace=0.025, hspace=0.025)
    fig.show()
    fig.savefig('sws_bwsalt.png')

    print('Processing sea ice formation')
    # Read seasonally averaged data
    id = Dataset(ice_file_beg, 'r')
    thdgr_nodes_beg = id.variables['thdgr'][:,:]*1e7
    id.close()
    id = Dataset(ice_file_end, 'r')
    thdgr_nodes_end = id.variables['thdgr'][:,:]*1e7
    id.close()
    thdgr_nodes_diff = thdgr_nodes_end - thdgr_nodes_beg
    # Now average over each element
    thdgr_beg = empty([4, num_elm_ocn])
    thdgr_end = empty([4, num_elm_ocn])
    thdgr_diff = empty([4, num_elm_ocn])
    i = 0
    for elm in elements:
        if not elm.cavity:
            thdgr_beg[:,i] = (thdgr_nodes_beg[:,elm.nodes[0].id] + thdgr_nodes_beg[:,elm.nodes[1].id] + thdgr_nodes_beg[:,elm.nodes[2].id])/3.0
            thdgr_end[:,i] = (thdgr_nodes_end[:,elm.nodes[0].id] + thdgr_nodes_end[:,elm.nodes[1].id] + thdgr_nodes_end[:,elm.nodes[2].id])/3.0
            thdgr_diff[:,i] = (thdgr_nodes_diff[:,elm.nodes[0].id] + thdgr_nodes_diff[:,elm.nodes[1].id] + thdgr_nodes_diff[:,elm.nodes[2].id])/3.0
            i += 1
    # Plot beginning and end
    fig = figure(figsize=(19,9))
    for season in range(4):
        # 1996-2005
        ax = fig.add_subplot(2, 4, season+1, aspect='equal')
        img = PatchCollection(patches_ocn, cmap='RdBu_r')
        img.set_array(thdgr_beg[season,:])
        img.set_edgecolor('face')
        img.set_clim(vmin=-2, vmax=2)
        ax.add_collection(img)
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        ax.set_xticks([])
        ax.set_yticks([])
        title(season_names[season], fontsize=24)
        if season == 0:
            text(-22, 10, '1996-2005', fontsize=20, ha='right', rotation=90)
        # 2091-2100
        ax = fig.add_subplot(2, 4, season+5, aspect='equal')
        img = PatchCollection(patches_ocn, cmap='RdBu_r')
        img.set_array(thdgr_end[season,:])
        img.set_edgecolor('face')
        img.set_clim(vmin=-2, vmax=2)
        ax.add_collection(img)
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        ax.set_xticks([])
        ax.set_yticks([])
        if season == 0:
            text(-22, 10, '2091-2100', fontsize=20, ha='right', rotation=90)
        if season == 3:
            cbaxes = fig.add_axes([0.35, 0.04, 0.3, 0.02])
            cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='both')
    suptitle(r'Sea ice thermodynamic growth rate (10$^{-7}$ m/s)', fontsize=24)
    subplots_adjust(wspace=0.025, hspace=0.025)
    fig.show()
    fig.savefig('sws_thdgr.png')
    # Plot difference
    fig = figure(figsize=(19,6))
    for season in range(4):
        ax = fig.add_subplot(1, 4, season+1, aspect='equal')
        img = PatchCollection(patches_ocn, cmap='RdBu_r')
        img.set_array(thdgr_diff[season,:])
        img.set_edgecolor('face')
        img.set_clim(vmin=-1, vmax=1)
        ax.add_collection(img)
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        ax.set_xticks([])
        ax.set_yticks([])
        title(season_names[season], fontsize=24)
        if season == 3:
            cbaxes = fig.add_axes([0.35, 0.04, 0.3, 0.02])
            cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='both')
    suptitle(r'Sea ice thermodynamic growth rate (10$^{-7}$ m/s), 2091-2100 minus 1996-2005', fontsize=24)
    subplots_adjust(wspace=0.025)
    fig.show()
    fig.savefig('sws_thdgr_diff.png')       

    print('Processing vertically averaged velocity')
    # Need to read some more of the grid    
    # Read number of 2D nodes
    f  = open(mesh_path + 'nod2d.out', 'r')
    n2d = int(f.readline())
    f.close()
    # Read rotated lat and lon for each 3D node, also depth
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
    # For lat and lon, only care about the 2D nodes (the first n2d indices)
    rlon = array(rlon[0:n2d])
    rlat = array(rlat[0:n2d])
    node_depth = array(node_depth)
    # Unrotate lat and lon
    lon, lat = unrotate_grid(rlon, rlat)
    # Calculate polar coordinates of each node
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)
    # Read lists of which nodes are directly below which
    f = open(mesh_path + 'aux3d.out', 'r')
    max_num_layers = int(f.readline())
    node_columns = zeros([n2d, max_num_layers])
    for n in range(n2d):
        for k in range(max_num_layers):
            node_columns[n,k] = int(f.readline())
    node_columns = node_columns.astype(int)
    f.close()
    # Now read the data (full 3D field for both u and v)
    id = Dataset(oce_file_beg, 'r')
    node_ur_3d_beg = id.variables['u'][0,:]
    node_vr_3d_beg = id.variables['v'][0,:]
    id.close()
    id = Dataset(oce_file_end, 'r')
    node_ur_3d_end = id.variables['u'][0,:]
    node_vr_3d_end = id.variables['v'][0,:]
    id.close()
    # Vertically average
    node_ur_beg = zeros(n2d)
    node_vr_beg = zeros(n2d)
    node_ur_end = zeros(n2d)
    node_vr_end = zeros(n2d)
    for n in range(n2d):
        # Integrate udz, vdz, and dz over this water column
        udz_col_beg = 0
        vdz_col_beg = 0
        udz_col_end = 0
        vdz_col_end = 0
        dz_col = 0
        for k in range(max_num_layers-1):
            if node_columns[n,k+1] == -999:
                # Reached the bottom
                break
            # Trapezoidal rule
            top_id = node_columns[n,k]
            bot_id = node_columns[n,k+1]
            dz_tmp = node_depth[bot_id-1] - node_depth[top_id-1]
            udz_col_beg += 0.5*(node_ur_3d_beg[top_id-1]+node_ur_3d_beg[bot_id-1])*dz_tmp
            vdz_col_beg += 0.5*(node_vr_3d_beg[top_id-1]+node_vr_3d_beg[bot_id-1])*dz_tmp
            udz_col_end += 0.5*(node_ur_3d_end[top_id-1]+node_ur_3d_end[bot_id-1])*dz_tmp
            vdz_col_end += 0.5*(node_vr_3d_end[top_id-1]+node_vr_3d_end[bot_id-1])*dz_tmp
            dz_col += dz_tmp
        # Convert from integrals to averages
        node_ur_beg[n] = udz_col_beg/dz_col
        node_vr_beg[n] = vdz_col_beg/dz_col
        node_ur_end[n] = udz_col_end/dz_col
        node_vr_end[n] = vdz_col_end/dz_col
    # Unrotate
    node_u_beg, node_v_beg = unrotate_vector(rlon, rlat, node_ur_beg, node_vr_beg)
    node_u_end, node_v_end = unrotate_vector(rlon, rlat, node_ur_end, node_vr_end)
    # Calculate speed
    node_speed_beg = sqrt(node_u_beg**2 + node_v_beg**2)
    node_speed_end = sqrt(node_u_end**2 + node_v_end**2)
    # Now calculate speed at each element, averaged over 3 corners
    speed_beg = []
    speed_end = []
    for elm in elements:
        speed_beg.append(mean([node_speed_beg[elm.nodes[0].id], node_speed_beg[elm.nodes[1].id], node_speed_beg[elm.nodes[2].id]]))
        speed_end.append(mean([node_speed_end[elm.nodes[0].id], node_speed_end[elm.nodes[1].id], node_speed_end[elm.nodes[2].id]]))
    # Set up bins for vectors
    x_bins = linspace(x_min, x_max, num=num_bins_x+1)
    y_bins = linspace(y_min, y_max, num=num_bins_y+1)
    # Calculate centres of bins (for plotting)
    x_centres = 0.5*(x_bins[:-1] + x_bins[1:])
    y_centres = 0.5*(y_bins[:-1] + y_bins[1:])
    # Set up arrays to integrate velocity in each bin
    # Simple averaging of all the points inside each bin
    ubin_beg = zeros([size(y_centres), size(x_centres)])
    vbin_beg = zeros([size(y_centres), size(x_centres)])
    ubin_end = zeros([size(y_centres), size(x_centres)])
    vbin_end = zeros([size(y_centres), size(x_centres)])
    num_pts = zeros([size(y_centres), size(x_centres)])
    # First convert to polar coordinates, rotate to account for
    # longitude in circumpolar projection, and convert back to vector
    # components
    theta_beg = arctan2(node_v_beg, node_u_beg)
    theta_circ_beg = theta_beg - lon*deg2rad
    u_circ_beg = node_speed_beg*cos(theta_circ_beg)
    v_circ_beg = node_speed_beg*sin(theta_circ_beg)
    theta_end = arctan2(node_v_end, node_u_end)
    theta_circ_end = theta_end - lon*deg2rad
    u_circ_end = node_speed_end*cos(theta_circ_end)
    v_circ_end = node_speed_end*sin(theta_circ_end)
    # Loop over 2D nodes
    for n in range(n2d):
        if x[n] > x_min and x[n] < x_max and y[n] > y_min and y[n] < y_max:
            # Figure out which bins this falls into
            x_index = nonzero(x_bins > x[n])[0][0]-1
            y_index = nonzero(y_bins > y[n])[0][0]-1
            # Integrate
            ubin_beg[y_index, x_index] += u_circ_beg[n]
            vbin_beg[y_index, x_index] += v_circ_beg[n]
            ubin_end[y_index, x_index] += u_circ_end[n]
            vbin_end[y_index, x_index] += v_circ_end[n]
            num_pts[y_index, x_index] += 1
    # Convert from sums to averages
    # First mask out points with no data
    ubin_beg = ma.masked_where(num_pts==0, ubin_beg)
    vbin_beg = ma.masked_where(num_pts==0, vbin_beg)
    ubin_end = ma.masked_where(num_pts==0, ubin_end)
    vbin_end = ma.masked_where(num_pts==0, vbin_end)
    # Divide everything else by number of points
    flag = num_pts > 0
    ubin_beg[flag] = ubin_beg[flag]/num_pts[flag]
    vbin_beg[flag] = vbin_beg[flag]/num_pts[flag]
    ubin_end[flag] = ubin_end[flag]/num_pts[flag]
    vbin_end[flag] = vbin_end[flag]/num_pts[flag]
    # Plot
    fig = figure(figsize=(16,7))
    # 1996-2005
    ax = fig.add_subplot(1, 2, 1, aspect='equal')
    img = PatchCollection(patches_all, cmap='cool')
    img.set_array(array(speed_beg))
    img.set_edgecolor('face')
    img.set_clim(vmin=0, vmax=0.25)
    ax.add_collection(img)
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax.add_collection(contours)
    quiver(x_centres, y_centres, ubin_beg, vbin_beg, scale=0.9, headwidth=8, headlength=9, color='black')
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    title('1996-2005', fontsize=20)
    # 2091-2100
    ax = fig.add_subplot(1, 2, 2, aspect='equal')
    img = PatchCollection(patches_all, cmap='cool')
    img.set_array(array(speed_end))
    img.set_edgecolor('face')
    img.set_clim(vmin=0, vmax=0.25)
    ax.add_collection(img)
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax.add_collection(contours)
    quiver(x_centres, y_centres, ubin_end, vbin_end, scale=0.9, headwidth=8, headlength=9, color='black')
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    title('2091-2100', fontsize=20)
    cbaxes = fig.add_axes([0.35, 0.04, 0.3, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='max')
    suptitle('Vertically averaged velocity (m/s)', fontsize=24)
    subplots_adjust(wspace=0.025, hspace=0.025)
    fig.show()
    fig.savefig('sws_velavg.png')


# Command-line interface
if __name__ == "__main__":

    sws_plots()
    

    
    

