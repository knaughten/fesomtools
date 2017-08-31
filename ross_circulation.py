from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from patches import *
from unrotate_vector import *
from unrotate_grid import *

def ross_circulation ():

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/high_res/'
    file_beg = '/short/y99/kaa561/FESOM/highres_spinup/annual_avg.oce.mean.1996.2005.nc'
    file_end = '/short/y99/kaa561/FESOM/rcp85_M_highres/output/annual_avg.oce.mean.2091.2100.nc'
    deg2rad = pi/180.0
    # Limits on x and y (polar coordinate transformation)
    x_min = -6
    x_max = 4
    y_min = -13
    y_max = -4.75
    num_bins_x = 30
    num_bins_y = 30

    print 'Building mesh'
    # Mask open ocean
    elements, mask_patches = make_patches(mesh_path, circumpolar=True, mask_cavities=True)
    # Unmask ice shelves
    patches = iceshelf_mask(elements)
    # The overlaid vectors are based on nodes not elements, so many of the
    # fesom_grid data structures fail to apply and we need to read some of the
    # mesh files again.
    # Read the cavity flag for each 2D surface node
    node_cavity = []
    f = open(mesh_path + 'cavity_flag_nod2d.out', 'r')
    for line in f:
        tmp = int(line)
        if tmp == 1:
            node_cavity.append(True)
        elif tmp == 0:
            node_cavity.append(False)
        else:
            print 'Problem'
    f.close()
    # Save the number of 2D nodes
    n2d = len(node_cavity)
    # Read rotated lat and lon for each node, also depth
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
    # Unrotate longitude
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
    # Count the number of elements in ice shelf cavities
    num_cavity_elm = 0
    for elm in elements:
        if elm.cavity:
            num_cavity_elm += 1
    num_elm = len(elements)
    num_ice_elm = num_elm - num_cavity_elm

    print 'Calculating vertically averaged velocity'
    velavg_beg = zeros([num_cavity_elm])
    velavg_end = zeros([num_cavity_elm])
    # Read full 3D fields for both u and v
    id = Dataset(file_beg, 'r')
    node_ur_3d_beg = id.variables['u'][0,:]
    node_vr_3d_beg = id.variables['v'][0,:]
    id.close()
    id = Dataset(file_end, 'r')
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
    # Calculate speed at each element, averaged over 3 corners
    i = 0
    for elm in elements:
        if elm.cavity:
            velavg_beg[i] = mean([node_speed_beg[elm.nodes[0].id], node_speed_beg[elm.nodes[1].id], node_speed_beg[elm.nodes[2].id]])
            velavg_end[i] = mean([node_speed_end[elm.nodes[0].id], node_speed_end[elm.nodes[1].id], node_speed_end[elm.nodes[2].id]])
            i += 1
    # Find bounds on speed in this region
    # Initialise with something impossible
    var_min = amax(array(velavg_beg))
    var_max = amin(array(velavg_beg))
    # Modify as needed
    i = 0
    for elm in elements:
        if elm.cavity:
            if any(elm.x >= x_min) and any(elm.x <= x_max) and any(elm.y >= y_min) and any(elm.y <= y_max):
                if velavg_beg[i] < var_min:
                    var_min = velavg_beg[i]
                if velavg_beg[i] > var_max:
                    var_max = velavg_beg[i]
                if velavg_end[i] < var_min:
                    var_min = velavg_end[i]
                if velavg_end[i] > var_max:
                    var_max = velavg_end[i]
            i += 1

    print 'Making vectors for overlay'
    # Set up bins (edges)
    x_bins = linspace(x_min, x_max, num=num_bins_x+1)
    y_bins = linspace(y_min, y_max, num=num_bins_y+1)
    # Calculate centres of bins (for plotting)
    x_centres = 0.5*(x_bins[:-1] + x_bins[1:])
    y_centres = 0.5*(y_bins[:-1] + y_bins[1:])
    # First set up arrays to integrate velocity in each bin
    # Simple averaging of all the points inside each bin
    ubin_beg = zeros([size(y_centres), size(x_centres)])
    vbin_beg = zeros([size(y_centres), size(x_centres)])
    ubin_end = zeros([size(y_centres), size(x_centres)])
    vbin_end = zeros([size(y_centres), size(x_centres)])
    num_pts = zeros([size(y_centres), size(x_centres)])
    # Convert to polar coordinates, rotate to account for longitude in
    # circumpolar projection, and convert back to vector components
    theta_beg = arctan2(node_v_beg, node_u_beg)
    theta_circ_beg = theta_beg - lon*deg2rad
    u_circ_beg = node_speed_beg*cos(theta_circ_beg)
    v_circ_beg = node_speed_beg*sin(theta_circ_beg)
    theta_end = arctan2(node_v_end, node_u_end)
    theta_circ_end = theta_end - lon*deg2rad
    u_circ_end = node_speed_end*cos(theta_circ_end)
    v_circ_end = node_speed_end*sin(theta_circ_end)
    # Loop over 2D nodes to fill in the velocity bins
    for n in range(n2d):
        if node_cavity[n]:
            if x[n] > x_min and x[n] < x_max and y[n] > y_min and y[n] < y_max:
                x_index = nonzero(x_bins > x[n])[0][0]-1
                y_index = nonzero(y_bins > y[n])[0][0]-1
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
    # Divide everything else by the number of points
    flag = num_pts > 0
    ubin_beg[flag] = ubin_beg[flag]/num_pts[flag]
    vbin_beg[flag] = vbin_beg[flag]/num_pts[flag]
    ubin_end[flag] = ubin_end[flag]/num_pts[flag]
    vbin_end[flag] = vbin_end[flag]/num_pts[flag]

    print 'Plotting'
    fig = figure(figsize=(18,9))
    fig.patch.set_facecolor('white')
    # Set up a grey square to fill the background with land
    x_reg, y_reg = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
    land_square = zeros(shape(x_reg))
    # 1996-2005
    ax = fig.add_subplot(1, 2, 1, aspect='equal')
    # Start with land background
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    # Add ice shelf elements
    img = PatchCollection(patches, cmap='cool')
    img.set_array(array(velavg_beg))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    # Mask out the open ocean in white
    overlay = PatchCollection(mask_patches, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    # Overlay vectors
    quiver(x_centres, y_centres, ubin_beg, vbin_beg, scale=0.9, headwidth=8, headlength=9, color='black')
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    title('1996-2005', fontsize=20)
    # RCP
    ax = fig.add_subplot(1, 2, 2, aspect='equal')
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches, cmap='cool')
    img.set_array(array(velavg_end))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    overlay = PatchCollection(mask_patches, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    quiver(x_centres, y_centres, ubin_end, vbin_end, scale=0.9, headwidth=8, headlength=9, color='black')
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    title('RCP 8.5 (2091-2100)', fontsize=20)
    # Colourbar on the right
    cbaxes = fig.add_axes([0.93, 0.3, 0.02, 0.4])
    cbar = colorbar(img, cax=cbaxes)
    # Main title
    suptitle('Ross Sea vertically averaged velocity (m/s)', fontsize=30)
    subplots_adjust(wspace=0.01)
    fig.show()
    fig.savefig('ross_circulation_change.png')


# Command-line interface
if __name__ == "__main__":

    ross_circulation()
    
    

    
    
