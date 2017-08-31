from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.colors import LinearSegmentedColormap
from fesom_grid import *
from patches import *
from unrotate_vector import *
from unrotate_grid import *

# For each major ice shelf, make a 2x1 plot of the given field for the low-res
# spinup (left) and high-res spinup (right), both averaged over the third
# repetition of the forcing, zoomed into the region of interest. Current
# options are ice shelf draft, ice shelf melt rate, bottom water temperature,
# bottom water salinity, surface velocity (with vectors overlaid), or
# vertically averaged velocity (with vectors overlaid).
# Input: var_name = 'draft', 'melt', 'temp', 'salt', 'vsfc', 'vavg'
def cavity_fields_res (var_name):

    # Paths to mesh directories
    mesh_path_low = '../FESOM/mesh/low_res/'
    mesh_path_high = '../FESOM/mesh/high_res/'
    # Paths to output files
    output_path_low = '../FESOM/lowres_spinup/rep3/'
    output_path_high = '../FESOM/highres_spinup/rep3/'
    if var_name == 'melt':
        file_name = 'annual_avg.forcing.diag.nc'
    elif var_name in ['temp', 'salt', 'vsfc', 'vavg']:
        file_name = 'annual_avg.oce.mean.nc'

    # Name of each ice shelf
    shelf_names = ['Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    # Beginnings of filenames for figures
    fig_heads = ['larsen_d', 'larsen_c', 'wilkins_georgevi_stange', 'ronne_filchner', 'abbot', 'pig', 'thwaites', 'dotson', 'getz', 'nickerson', 'sulzberger', 'mertz', 'totten_moscowuni', 'shackleton', 'west', 'amery', 'prince_harald', 'baudouin_borchgrevink', 'lazarev', 'nivl', 'fimbul_jelbart_ekstrom', 'brunt_riiser_larsen', 'ross']
    # Limits on longitude and latitude for each ice shelf
    # Note Ross crosses 180W=180E
    lon_min = [-62.67, -65.5, -79.17, -85, -104.17, -102.5, -108.33, -114.5, -135.67, -149.17, -155, 144, 115, 94.17, 80.83, 65, 33.83, 19, 12.9, 9.33, -10.05, -28.33, 158.33]
    lon_max = [-59.33, -60, -66.67, -28.33, -88.83, -99.17, -103.33, -111.5, -114.33, -140, -145, 146.62, 123.33, 102.5, 89.17, 75, 37.67, 33.33, 16.17, 12.88, 7.6, -10.33, -146.67]
    lat_min = [-73.03, -69.35, -74.17, -83.5, -73.28, -75.5, -75.5, -75.33, -74.9, -76.42, -78, -67.83, -67.17, -66.67, -67.83, -73.67, -69.83, -71.67, -70.5, -70.75, -71.83, -76.33, -85]
    lat_max = [-69.37, -66.13, -69.5, -74.67, -71.67, -74.17, -74.67, -73.67, -73, -75.17, -76.41, -66.67, -66.5, -64.83, -66.17, -68.33, -68.67, -68.33, -69.33, -69.83, -69.33, -71.5, -77]
    num_shelves = len(shelf_names)

    # Constants
    sec_per_year = 365.25*24*3600
    deg2rad = pi/180.0
    # Number of bins in each direction for vector overlay
    num_bins = 50

    print 'Building FESOM mesh'
    # Mask open ocean
    elements_low, mask_patches_low = make_patches(mesh_path_low, circumpolar=True, mask_cavities=True)
    elements_high, mask_patches_high = make_patches(mesh_path_high, circumpolar=True, mask_cavities=True)
    # Unmask ice shelves
    patches_low = iceshelf_mask(elements_low)
    patches_high = iceshelf_mask(elements_high)
    if var_name == 'draft':
        # Nothing more to read
        pass
    else:
        print 'Reading data'
        id_low = Dataset(output_path_low + file_name, 'r')
        id_high = Dataset(output_path_high + file_name, 'r')
        if var_name == 'melt':
            # Convert from m/s to m/y
            node_data_low = id_low.variables['wnet'][0,:]*sec_per_year
            node_data_high = id_high.variables['wnet'][0,:]*sec_per_year
        elif var_name == 'temp':
            # Read full 3D field for now
            node_data_low = id_low.variables['temp'][0,:]
            node_data_high = id_high.variables['temp'][0,:]
        elif var_name == 'salt':
            # Read full 3D field for now
            node_data_low = id_low.variables['salt'][0,:]
            node_data_high = id_high.variables['salt'][0,:]
        elif var_name in ['vsfc', 'vavg']:
            # The overlaid vectors are based on nodes not elements, so many
            # of the fesom_grid data structures fail to apply and we need to
            # read some of the FESOM grid files again.
            # Read the cavity flag for each 2D surface node
            cavity_low = []
            f = open(mesh_path_low + 'cavity_flag_nod2d.out', 'r')
            for line in f:
                tmp = int(line)
                if tmp == 1:
                    cavity_low.append(True)
                elif tmp == 0:
                    cavity_low.append(False)
                else:
                    print 'Problem'
                    #return
            f.close()
            cavity_high = []
            f = open(mesh_path_high + 'cavity_flag_nod2d.out', 'r')
            for line in f:
                tmp = int(line)
                if tmp == 1:
                    cavity_high.append(True)
                elif tmp == 0:
                    cavity_high.append(False)
                else:
                    print 'Problem'
                    #return
            f.close()
            # Save the number of 2D nodes
            n2d_low = len(cavity_low)
            n2d_high = len(cavity_high)
            # Read rotated lat and lon for each node; also read depth which is
            # needed for vertically averaged velocity
            f = open(mesh_path_low + 'nod3d.out', 'r')
            f.readline()
            rlon_low = []
            rlat_low = []
            node_depth_low = []
            for line in f:
                tmp = line.split()
                lon_tmp = float(tmp[1])
                lat_tmp = float(tmp[2])
                node_depth_tmp = -1*float(tmp[3])
                if lon_tmp < -180:
                    lon_tmp += 360
                elif lon_tmp > 180:
                    lon_tmp -= 360
                rlon_low.append(lon_tmp)
                rlat_low.append(lat_tmp)
                node_depth_low.append(node_depth_tmp)
            f.close()
            # For lat and lon, only care about the 2D nodes (the first n2d
            # indices)
            rlon_low = array(rlon_low[0:n2d_low])
            rlat_low = array(rlat_low[0:n2d_low])
            node_depth_low = array(node_depth_low)
            # Repeat for high resolution
            f = open(mesh_path_high + 'nod3d.out', 'r')
            f.readline()
            rlon_high = []
            rlat_high = []
            node_depth_high = []
            for line in f:
                tmp = line.split()
                lon_tmp = float(tmp[1])
                lat_tmp = float(tmp[2])
                node_depth_tmp = -1*float(tmp[3])
                if lon_tmp < -180:
                    lon_tmp += 360
                elif lon_tmp > 180:
                    lon_tmp -= 360
                rlon_high.append(lon_tmp)
                rlat_high.append(lat_tmp)
                node_depth_high.append(node_depth_tmp)
            f.close()
            rlon_high = array(rlon_high[0:n2d_high])
            rlat_high = array(rlat_high[0:n2d_high])
            node_depth_high = array(node_depth_high)
            # Unrotate longitude
            lon_low, lat_low = unrotate_grid(rlon_low, rlat_low)
            lon_high, lat_high = unrotate_grid(rlon_high, rlat_high)
            # Calculate polar coordinates of each node
            x_low = -(lat_low+90)*cos(lon_low*deg2rad+pi/2)
            y_low = (lat_low+90)*sin(lon_low*deg2rad+pi/2)
            x_high = -(lat_high+90)*cos(lon_high*deg2rad+pi/2)
            y_high = (lat_high+90)*sin(lon_high*deg2rad+pi/2)
            if var_name == 'vavg':
                # Read lists of which nodes are directly below which
                f = open(mesh_path_low + 'aux3d.out', 'r')
                max_num_layers_low = int(f.readline())
                node_columns_low = zeros([n2d_low, max_num_layers_low])
                for n in range(n2d_low):
                    for k in range(max_num_layers_low):
                        node_columns_low[n,k] = int(f.readline())
                node_columns_low = node_columns_low.astype(int)
                f.close()
                # Repeat for high resolution
                f = open(mesh_path_high + 'aux3d.out', 'r')
                max_num_layers_high = int(f.readline())
                node_columns_high = zeros([n2d_high, max_num_layers_high])
                for n in range(n2d_high):
                    for k in range(max_num_layers_high):
                        node_columns_high[n,k] = int(f.readline())
                node_columns_high = node_columns_high.astype(int)
                f.close()
            # Now we can actually read the data
            # Read full 3D field for both u and v
            node_ur_3d_low = id_low.variables['u'][0,:]
            node_vr_3d_low = id_low.variables['v'][0,:]
            node_ur_3d_high = id_high.variables['u'][0,:]
            node_vr_3d_high = id_high.variables['v'][0,:]
            if var_name == 'vsfc':
                # Only care about the first n2d nodes
                node_ur_low = node_ur_3d_low[0:n2d_low]
                node_vr_low = node_vr_3d_low[0:n2d_low]
                node_ur_high = node_ur_3d_high[0:n2d_high]
                node_vr_high = node_vr_3d_high[0:n2d_high]
            elif var_name == 'vavg':
                # Vertically average
                node_ur_low = zeros(n2d_low)
                node_vr_low = zeros(n2d_low)
                for n in range(n2d_low):
                    # Integrate udz, vdz, and dz over this water column
                    udz_col = 0
                    vdz_col = 0
                    dz_col = 0
                    for k in range(max_num_layers_low-1):
                        if node_columns_low[n,k+1] == -999:
                            # Reached the bottom
                            break
                        # Trapezoidal rule
                        top_id = node_columns_low[n,k]
                        bot_id = node_columns_low[n,k+1]
                        dz_tmp = node_depth_low[bot_id-1] - node_depth_low[top_id-1]
                        udz_col += 0.5*(node_ur_3d_low[top_id-1]+node_ur_3d_low[bot_id-1])*dz_tmp
                        vdz_col += 0.5*(node_vr_3d_low[top_id-1]+node_vr_3d_low[bot_id-1])*dz_tmp
                        dz_col += dz_tmp
                    # Convert from integrals to averages
                    node_ur_low[n] = udz_col/dz_col
                    node_vr_low[n] = vdz_col/dz_col
                # Repeat for high resolution
                node_ur_high = zeros(n2d_high)
                node_vr_high = zeros(n2d_high)
                for n in range(n2d_high):
                    udz_col = 0
                    vdz_col = 0
                    dz_col = 0
                    for k in range(max_num_layers_high-1):
                        if node_columns_high[n,k+1] == -999:
                            break
                        top_id = node_columns_high[n,k]
                        bot_id = node_columns_high[n,k+1]
                        dz_tmp = node_depth_high[bot_id-1] - node_depth_high[top_id-1]
                        udz_col += 0.5*(node_ur_3d_high[top_id-1]+node_ur_3d_high[bot_id-1])*dz_tmp
                        vdz_col += 0.5*(node_vr_3d_high[top_id-1]+node_vr_3d_high[bot_id-1])*dz_tmp
                        dz_col += dz_tmp
                    node_ur_high[n] = udz_col/dz_col
                    node_vr_high[n] = vdz_col/dz_col
            # Unrotate
            node_u_low, node_v_low = unrotate_vector(rlon_low, rlat_low, node_ur_low, node_vr_low)
            node_u_high, node_v_high = unrotate_vector(rlon_high, rlat_high, node_ur_high, node_vr_high)
            # Calculate speed
            node_data_low = sqrt(node_u_low**2 + node_v_low**2)
            node_data_high = sqrt(node_u_high**2 + node_v_high**2)
        id_low.close()
        id_high.close()
    # Calculate given field at each element
    data_low = []
    for elm in elements_low:
        # For each element in an ice shelf cavity, append the mean value
        # for the 3 component Nodes
        if elm.cavity:
            if var_name == 'draft':
                # Ice shelf draft is depth of surface layer
                data_low.append(mean([elm.nodes[0].depth, elm.nodes[1].depth, elm.nodes[2].depth]))
            elif var_name in ['melt', 'vsfc', 'vavg']:
                # Surface nodes (or 2D in the case of vavg)
                data_low.append(mean([node_data_low[elm.nodes[0].id], node_data_low[elm.nodes[1].id], node_data_low[elm.nodes[2].id]]))
            elif var_name in ['temp', 'salt']:
                # Bottom nodes
                data_low.append(mean([node_data_low[elm.nodes[0].find_bottom().id], node_data_low[elm.nodes[1].find_bottom().id], node_data_low[elm.nodes[2].find_bottom().id]]))
    # Repeat for high resolution
    data_high = []
    for elm in elements_high:
        if elm.cavity:
            if var_name == 'draft':
                data_high.append(mean([elm.nodes[0].depth, elm.nodes[1].depth, elm.nodes[2].depth]))
            elif var_name in ['melt', 'vsfc', 'vavg']:
                data_high.append(mean([node_data_high[elm.nodes[0].id], node_data_high[elm.nodes[1].id], node_data_high[elm.nodes[2].id]]))
            elif var_name in ['temp', 'salt']:
                data_high.append(mean([node_data_high[elm.nodes[0].find_bottom().id], node_data_high[elm.nodes[1].find_bottom().id], node_data_high[elm.nodes[2].find_bottom().id]]))

    # Loop over ice shelves
    for index in range(num_shelves):
        print 'Processing ' + shelf_names[index]
        # Convert lat/lon bounds to polar coordinates for plotting
        x1 = -(lat_min[index]+90)*cos(lon_min[index]*deg2rad+pi/2)
        y1 = (lat_min[index]+90)*sin(lon_min[index]*deg2rad+pi/2)
        x2 = -(lat_min[index]+90)*cos(lon_max[index]*deg2rad+pi/2)
        y2 = (lat_min[index]+90)*sin(lon_max[index]*deg2rad+pi/2)
        x3 = -(lat_max[index]+90)*cos(lon_min[index]*deg2rad+pi/2)
        y3 = (lat_max[index]+90)*sin(lon_min[index]*deg2rad+pi/2)
        x4 = -(lat_max[index]+90)*cos(lon_max[index]*deg2rad+pi/2)
        y4 = (lat_max[index]+90)*sin(lon_max[index]*deg2rad+pi/2)
        # Find the new bounds on x and y
        x_min = amin(array([x1, x2, x3, x4]))
        x_max = amax(array([x1, x2, x3, x4]))
        y_min = amin(array([y1, y2, y3, y4]))
        y_max = amax(array([y1, y2, y3, y4]))
        # Now make the plot square: enlarge the smaller of delta_x and delta_y
        # so they are equal
        delta_x = x_max - x_min
        delta_y = y_max - y_min
        if delta_x > delta_y:
            diff = 0.5*(delta_x - delta_y)
            y_min -= diff
            y_max += diff
        elif delta_y > delta_x:
            diff = 0.5*(delta_y - delta_x)
            x_min -= diff
            x_max += diff
        # Set up a grey square to fill the background with land
        x_reg, y_reg = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
        land_square = zeros(shape(x_reg))
        # Find bounds on variables in this region
        # Start with something impossible
        var_min = amax(data_low)
        var_max = amin(data_low)
        # Modify with low-res
        i = 0
        for elm in elements_low:
            if elm.cavity:
                if any(elm.x >= x_min) and any(elm.x <= x_max) and any(elm.y >= y_min) and any(elm.y <= y_max):
                    if data_low[i] < var_min:
                        var_min = data_low[i]
                    if data_low[i] > var_max:
                        var_max = data_low[i]
                i += 1
        # Modify with high-res
        i = 0
        for elm in elements_high:
            if elm.cavity:
                if any(elm.x >= x_min) and any(elm.x <= x_max) and any(elm.y >= y_min) and any(elm.y <= y_max):
                    if data_high[i] < var_min:
                        var_min = data_high[i]
                    if data_high[i] > var_max:
                        var_max = data_high[i]
                i += 1
        if var_name == 'melt':
            # Special colour map
            if var_min < 0:
                # There is refreezing here; include blue for elements below 0
                cmap_vals = array([var_min, 0, 0.25*var_max, 0.5*var_max, 0.75*var_max, var_max])
                cmap_colors = [(0.26, 0.45, 0.86), (1, 1, 1), (1, 0.9, 0.4), (0.99, 0.59, 0.18), (0.5, 0.0, 0.08), (0.96, 0.17, 0.89)]
                cmap_vals_norm = (cmap_vals - var_min)/(var_max - var_min)
                cmap_list = []
                for i in range(size(cmap_vals)):
                    cmap_list.append((cmap_vals_norm[i], cmap_colors[i]))
                mf_cmap = LinearSegmentedColormap.from_list('melt_freeze', cmap_list)
            else:
                # No refreezing
                cmap_vals = array([0, 0.25*var_max, 0.5*var_max, 0.75*var_max, var_max])
                cmap_colors = [(1, 1, 1), (1, 0.9, 0.4), (0.99, 0.59, 0.18), (0.5, 0.0, 0.08), (0.96, 0.17, 0.89)]
                cmap_vals_norm = cmap_vals/var_max
                cmap_list = []
                for i in range(size(cmap_vals)):
                    cmap_list.append((cmap_vals_norm[i], cmap_colors[i]))
                mf_cmap = LinearSegmentedColormap.from_list('melt_freeze', cmap_list)
            colour_map = mf_cmap
        else:
            colour_map = 'jet'
        if var_name in ['vsfc', 'vavg']:
            # Make vectors for overlay
            # Set up bins (edges)
            x_bins = linspace(x_min, x_max, num=num_bins+1)
            y_bins = linspace(y_min, y_max, num=num_bins+1)
            # Calculate centres of bins (for plotting)
            x_centres = 0.5*(x_bins[:-1] + x_bins[1:])
            y_centres = 0.5*(y_bins[:-1] + y_bins[1:])
            # Low resolution
            # First set up arrays to integrate velocity in each bin
            # Simple averaging of all the points inside each bin
            u_low = zeros([size(y_centres), size(x_centres)])
            v_low = zeros([size(y_centres), size(x_centres)])
            num_pts_low = zeros([size(y_centres), size(x_centres)])
            # Convert to polar coordinates, rotate to account for longitude in
            # circumpolar projection, and convert back to vector components
            theta_low = arctan2(node_v_low, node_u_low)
            theta_circ_low = theta_low - lon_low*deg2rad
            u_circ_low = node_data_low*cos(theta_circ_low)  # node_data_low is speed
            v_circ_low = node_data_low*sin(theta_circ_low)
            # Loop over 2D nodes to fill in the velocity bins
            for n in range(n2d_low):
                # Make sure we're in an ice shelf cavity
                if cavity_low[n]:
                    # Check if we're in the region of interest
                    if x_low[n] > x_min and x_low[n] < x_max and y_low[n] > y_min and y_low[n] < y_max:
                        # Figure out which bins this falls into
                        x_index = nonzero(x_bins > x_low[n])[0][0] - 1
                        y_index = nonzero(y_bins > y_low[n])[0][0] - 1
                        # Integrate
                        u_low[y_index, x_index] += u_circ_low[n]
                        v_low[y_index, x_index] += v_circ_low[n]
                        num_pts_low[y_index, x_index] += 1
            # Convert from sums to averages
            # First mask out points with no data
            u_low = ma.masked_where(num_pts_low==0, u_low)
            v_low = ma.masked_where(num_pts_low==0, v_low)
            # Divide everything else by the number of points
            flag = num_pts_low > 0
            u_low[flag] = u_low[flag]/num_pts_low[flag]
            v_low[flag] = v_low[flag]/num_pts_low[flag]
            # Repeat for high resolution
            u_high = zeros([size(y_centres), size(x_centres)])
            v_high = zeros([size(y_centres), size(x_centres)])
            num_pts_high = zeros([size(y_centres), size(x_centres)])
            theta_high = arctan2(node_v_high, node_u_high)
            theta_circ_high = theta_high - lon_high*deg2rad
            u_circ_high = node_data_high*cos(theta_circ_high)
            v_circ_high = node_data_high*sin(theta_circ_high)
            for n in range(n2d_high):
                if cavity_high[n]:
                    if x_high[n] > x_min and x_high[n] < x_max and y_high[n] > y_min and y_high[n] < y_max:
                        x_index = nonzero(x_bins > x_high[n])[0][0] - 1
                        y_index = nonzero(y_bins > y_high[n])[0][0] - 1
                        u_high[y_index, x_index] += u_circ_high[n]
                        v_high[y_index, x_index] += v_circ_high[n]
                        num_pts_high[y_index, x_index] += 1
            u_high = ma.masked_where(num_pts_high==0, u_high)
            v_high = ma.masked_where(num_pts_high==0, v_high)
            flag = num_pts_high > 0
            u_high[flag] = u_high[flag]/num_pts_high[flag]
            v_high[flag] = v_high[flag]/num_pts_high[flag]
        # Plot
        fig = figure(figsize=(30,12))
        # Low resolution
        ax1 = fig.add_subplot(1,2,1, aspect='equal')
        # Start with land background
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        # Add ice shelf elements
        img1 = PatchCollection(patches_low, cmap=colour_map)
        img1.set_array(array(data_low))
        img1.set_edgecolor('face')
        img1.set_clim(vmin=var_min, vmax=var_max)
        ax1.add_collection(img1)
        # Mask out the open ocean in white
        overlay1 = PatchCollection(mask_patches_low, facecolor=(1,1,1))
        overlay1.set_edgecolor('face')
        ax1.add_collection(overlay1)
        if var_name in ['vsfc', 'vavg']:
            # Overlay vectors
            quiver(x_centres, y_centres, u_low, v_low, scale=1.5, color='black')
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        axis('off')
        title('Low-res', fontsize=24)
        # Repeat for high resolution
        ax2 = fig.add_subplot(1,2,2, aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img2 = PatchCollection(patches_high, cmap=colour_map)
        img2.set_array(array(data_high))
        img2.set_edgecolor('face')
        img2.set_clim(vmin=var_min, vmax=var_max)
        ax2.add_collection(img2)
        overlay2 = PatchCollection(mask_patches_high, facecolor=(1,1,1))
        overlay2.set_edgecolor('face')
        ax2.add_collection(overlay2)
        if var_name in ['vsfc', 'vavg']:
            quiver(x_centres, y_centres, u_high, v_high, scale=1.5, color='black')
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        axis('off')
        title('High-res', fontsize=24)
        # Colourbar on the right
        cbaxes = fig.add_axes([0.92, 0.2, 0.01, 0.6])
        cbar = colorbar(img2, cax=cbaxes)
        cbar.ax.tick_params(labelsize=20)
        # Main title
        if var_name == 'draft':
            title_string = ' draft (m)'
        elif var_name == 'melt':
            title_string = ' melt rate (m/y)'
        elif var_name == 'temp':
            title_string = r' bottom water temperature ($^{\circ}$C)'
        elif var_name == 'salt':
            title_string = ' bottom water salinity (psu)'
        elif var_name == 'vsfc':
            title_string = ' surface velocity (m/s)'
        elif var_name == 'vavg':
            title_string = ' vertically averaged velocity (m/s)'
        suptitle(shelf_names[index] + title_string, fontsize=30)
        subplots_adjust(wspace=0.05)
        #fig.show()
        fig.savefig(fig_heads[index] + '_' + var_name + '.png')


# Command-line interface
if __name__ == "__main__":

    var_key = int(raw_input("Ice shelf draft (1), melt rate (2), bottom water temperature (3), bottom water salinity (4), surface velocity (5), or vertically averaged velocity (6)? "))
    if var_key == 1:
        var_name = 'draft'
    elif var_key == 2:
        var_name = 'melt'
    elif var_key == 3:
        var_name = 'temp'
    elif var_key == 4:
        var_name = 'salt'
    elif var_key == 5:
        var_name = 'vsfc'
    elif var_key == 6:
        var_name = 'vavg'
    cavity_fields_res(var_name)
            
                    
            

    
