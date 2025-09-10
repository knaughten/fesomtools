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

# For each major ice shelf, make a 2x1 plot of the given field averaged over
# the first 10 years of the RCP (left) and the last 10 years (right), zoomed
# into the region of interest. Current options are ice shelf melt rate, bottom
# water temperature, bottom water salinity, surface velocity (with vectors
# overlaid), or vertically averaged velocity (with vectors overlaid).
# Input:
# var_name = 'melt', 'temp', 'salt', 'vsfc', 'vavg'
# mesh_path = path to FESOM mesh directory
# output_path = path to output directory for that RCP, containing the filenames
#               specified in file_name_beg and file_name_end (either
#               forcing.diag.nc or oce.mean.nc averaged over the first and last
#               10 years of the RCP)
# fig_dir = path to directory for saving figures (must end in '/' unless it is
#           the current directory, i.e. empty string)
def cavity_fields_rcp (var_name, mesh_path, output_path, fig_dir=''):

    if var_name == 'melt':
        file_name_beg = 'annual_avg.forcing.diag.2006.2015.nc'
        file_name_end = 'annual_avg.forcing.diag.2091.2100.nc'
    elif var_name in ['temp', 'salt', 'vsfc', 'vavg']:
        file_name_beg = 'annual_avg.oce.mean.2006.2015.nc'
        file_name_end = 'annual_avg.oce.mean.2091.2100.nc'

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

    print('Building FESOM mesh')
    # Mask open ocean
    elements, mask_patches = make_patches(mesh_path, circumpolar=True, mask_cavities=True)
    # Unmask ice shelves
    patches = iceshelf_mask(elements)
    print('Reading data')
    id_beg = Dataset(output_path + file_name_beg, 'r')
    id_end = Dataset(output_path + file_name_end, 'r')
    if var_name == 'melt':
        # Convert from m/s to m/y
        node_data_beg = id_beg.variables['wnet'][0,:]*sec_per_year
        node_data_end = id_end.variables['wnet'][0,:]*sec_per_year
    elif var_name == 'temp':
        # Read full 3D field for now
        node_data_beg = id_beg.variables['temp'][0,:]
        node_data_end = id_end.variables['temp'][0,:]
    elif var_name == 'salt':
        # Read full 3D field for now
        node_data_beg = id_beg.variables['salt'][0,:]
        node_data_end = id_end.variables['salt'][0,:]
    elif var_name in ['vsfc', 'vavg']:
        # The overlaid vectors are based on nodes not elements, so many
        # of the fesom_grid data structures fail to apply and we need to
        # read some of the FESOM grid files again.
        # Read the cavity flag for each 2D surface node
        cavity = []
        f = open(mesh_path + 'cavity_flag_nod2d.out', 'r')
        for line in f:
            tmp = int(line)
            if tmp == 1:
                cavity.append(True)
            elif tmp == 0:
                cavity.append(False)
            else:
                print('Problem')
                #return
        f.close()
        # Save the number of 2D nodes
        n2d = len(cavity)
        # Read rotated lat and lon for each node; also read depth which is
        # needed for vertically averaged velocity
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
        # For lat and lon, only care about the 2D nodes (the first n2d
        # indices)
        rlon = array(rlon[0:n2d])
        rlat = array(rlat[0:n2d])
        node_depth = array(node_depth)
        # Unrotate longitude
        lon, lat = unrotate_grid(rlon, rlat)
        # Calculate polar coordinates of each node
        x = -(lat+90)*cos(lon*deg2rad+pi/2)
        y = (lat+90)*sin(lon*deg2rad+pi/2)
        if var_name == 'vavg':
            # Read lists of which nodes are directly below which
            f = open(mesh_path + 'aux3d.out', 'r')
            max_num_layers = int(f.readline())
            node_columns = zeros([n2d, max_num_layers])
            for n in range(n2d):
                for k in range(max_num_layers):
                    node_columns[n,k] = int(f.readline())
            node_columns = node_columns.astype(int)
            f.close()
        # Now we can actually read the data
        # Read full 3D field for both u and v
        node_ur_3d_beg = id_beg.variables['u'][0,:]
        node_vr_3d_beg = id_beg.variables['v'][0,:]
        node_ur_3d_end = id_end.variables['u'][0,:]
        node_vr_3d_end = id_end.variables['v'][0,:]
        if var_name == 'vsfc':
            # Only care about the first n2d nodes
            node_ur_beg = node_ur_3d_beg[0:n2d]
            node_vr_beg = node_vr_3d_beg[0:n2d]
            node_ur_end = node_ur_3d_end[0:n2d]
            node_vr_end = node_vr_3d_end[0:n2d]
        elif var_name == 'vavg':
            # Vertically average
            node_ur_beg = zeros(n2d)
            node_vr_beg = zeros(n2d)
            node_ur_end = zeros(n2d)
            node_vr_end = zeros(n2d)
            for n in range(n2d):
                # Integrate udz, vdz, and dz over this water column
                udz_beg_col = 0
                vdz_beg_col = 0
                udz_end_col = 0
                vdz_end_col = 0
                dz_col = 0
                for k in range(max_num_layers-1):
                    if node_columns[n,k+1] == -999:
                        # Reached the bottom
                        break
                    # Trapezoidal rule
                    top_id = node_columns[n,k]
                    bot_id = node_columns[n,k+1]
                    dz_tmp = node_depth[bot_id-1] - node_depth[top_id-1]
                    udz_beg_col += 0.5*(node_ur_3d_beg[top_id-1]+node_ur_3d_beg[bot_id-1])*dz_tmp
                    vdz_beg_col += 0.5*(node_vr_3d_beg[top_id-1]+node_vr_3d_beg[bot_id-1])*dz_tmp
                    udz_end_col += 0.5*(node_ur_3d_end[top_id-1]+node_ur_3d_end[bot_id-1])*dz_tmp
                    vdz_end_col += 0.5*(node_vr_3d_end[top_id-1]+node_vr_3d_end[bot_id-1])*dz_tmp
                    dz_col += dz_tmp
                # Convert from integrals to averages
                node_ur_beg[n] = udz_beg_col/dz_col
                node_vr_beg[n] = vdz_beg_col/dz_col
                node_ur_end[n] = udz_end_col/dz_col
                node_vr_end[n] = vdz_end_col/dz_col
        # Unrotate
        node_u_beg, node_v_beg = unrotate_vector(rlon, rlat, node_ur_beg, node_vr_beg)
        node_u_end, node_v_end = unrotate_vector(rlon, rlat, node_ur_end, node_vr_end)
        # Calculate speed
        node_data_beg = sqrt(node_u_beg**2 + node_v_beg**2)
        node_data_end = sqrt(node_u_end**2 + node_v_end**2)
    id_beg.close()
    id_end.close()
    # Calculate given field at each element
    data_beg = []
    data_end = []
    for elm in elements:
        # For each element in an ice shelf cavity, append the mean value
        # for the 3 component Nodes
        if elm.cavity:
            if var_name in ['melt', 'vsfc', 'vavg']:
                # Surface nodes (or 2D in the case of vavg)
                data_beg.append(mean([node_data_beg[elm.nodes[0].id], node_data_beg[elm.nodes[1].id], node_data_beg[elm.nodes[2].id]]))
                data_end.append(mean([node_data_end[elm.nodes[0].id], node_data_end[elm.nodes[1].id], node_data_end[elm.nodes[2].id]]))
            elif var_name in ['temp', 'salt']:
                # Bottom nodes
                data_beg.append(mean([node_data_beg[elm.nodes[0].find_bottom().id], node_data_beg[elm.nodes[1].find_bottom().id], node_data_beg[elm.nodes[2].find_bottom().id]]))
                data_end.append(mean([node_data_end[elm.nodes[0].find_bottom().id], node_data_end[elm.nodes[1].find_bottom().id], node_data_end[elm.nodes[2].find_bottom().id]]))

    # Loop over ice shelves
    for index in range(num_shelves):
        print('Processing ' + shelf_names[index])
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
        var_min = amax(data_beg)
        var_max = amin(data_beg)
        # Modify
        i = 0
        for elm in elements:
            if elm.cavity:
                if any(elm.x >= x_min) and any(elm.x <= x_max) and any(elm.y >= y_min) and any(elm.y <= y_max):
                    if data_beg[i] < var_min:
                        var_min = data_beg[i]
                    if data_beg[i] > var_max:
                        var_max = data_beg[i]
                    if data_end[i] < var_min:
                        var_min = data_end[i]
                    if data_end[i] > var_max:
                        var_max = data_end[i]
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
            # Set up arrays to integrate velocity in each bin
            # Simple averaging of all the points inside each bin
            u_beg = zeros([size(y_centres), size(x_centres)])
            v_beg = zeros([size(y_centres), size(x_centres)])
            u_end = zeros([size(y_centres), size(x_centres)])
            v_end = zeros([size(y_centres), size(x_centres)])
            num_pts = zeros([size(y_centres), size(x_centres)])
            # Convert to polar coordinates, rotate to account for longitude in
            # circumpolar projection, and convert back to vector components
            theta_beg = arctan2(node_v_beg, node_u_beg)
            theta_circ_beg = theta_beg - lon*deg2rad
            u_circ_beg = node_data_beg*cos(theta_circ_beg)  # node_data_beg is speed
            v_circ_beg = node_data_beg*sin(theta_circ_beg)
            theta_end = arctan2(node_v_end, node_u_end)
            theta_circ_end = theta_end - lon*deg2rad
            u_circ_end = node_data_end*cos(theta_circ_end)
            v_circ_end = node_data_end*sin(theta_circ_end)
            # Loop over 2D nodes to fill in the velocity bins
            for n in range(n2d):
                # Make sure we're in an ice shelf cavity
                if cavity[n]:
                    # Check if we're in the region of interest
                    if x[n] > x_min and x[n] < x_max and y[n] > y_min and y[n] < y_max:
                        # Figure out which bins this falls into
                        x_index = nonzero(x_bins > x[n])[0][0] - 1
                        y_index = nonzero(y_bins > y[n])[0][0] - 1
                        # Integrate
                        u_beg[y_index, x_index] += u_circ_beg[n]
                        v_beg[y_index, x_index] += v_circ_beg[n]
                        u_end[y_index, x_index] += u_circ_end[n]
                        v_end[y_index, x_index] += v_circ_end[n]
                        num_pts[y_index, x_index] += 1
            # Convert from sums to averages
            # First mask out points with no data
            u_beg = ma.masked_where(num_pts==0, u_beg)
            v_beg = ma.masked_where(num_pts==0, v_beg)
            u_end = ma.masked_where(num_pts==0, u_end)
            v_end = ma.masked_where(num_pts==0, v_end)
            # Divide everything else by the number of points
            flag = num_pts > 0
            u_beg[flag] = u_beg[flag]/num_pts[flag]
            v_beg[flag] = v_beg[flag]/num_pts[flag]
            u_end[flag] = u_end[flag]/num_pts[flag]
            v_end[flag] = v_end[flag]/num_pts[flag]
        # Plot
        fig = figure(figsize=(30,12))
        # First 10 years
        ax1 = fig.add_subplot(1,2,1, aspect='equal')
        # Start with land background
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        # Add ice shelf elements
        img1 = PatchCollection(patches, cmap=colour_map)
        img1.set_array(array(data_beg))
        img1.set_edgecolor('face')
        img1.set_clim(vmin=var_min, vmax=var_max)
        ax1.add_collection(img1)
        # Mask out the open ocean in white
        overlay1 = PatchCollection(mask_patches, facecolor=(1,1,1))
        overlay1.set_edgecolor('face')
        ax1.add_collection(overlay1)
        if var_name in ['vsfc', 'vavg']:
            # Overlay vectors
            quiver(x_centres, y_centres, u_beg, v_beg, scale=1.5, color='black')
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        axis('off')
        title('2006-2015', fontsize=24)
        # Repeat for last 10 years
        ax2 = fig.add_subplot(1,2,2, aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img2 = PatchCollection(patches, cmap=colour_map)
        img2.set_array(array(data_end))
        img2.set_edgecolor('face')
        img2.set_clim(vmin=var_min, vmax=var_max)
        ax2.add_collection(img2)
        overlay2 = PatchCollection(mask_patches, facecolor=(1,1,1))
        overlay2.set_edgecolor('face')
        ax2.add_collection(overlay2)
        if var_name in ['vsfc', 'vavg']:
            quiver(x_centres, y_centres, u_end, v_end, scale=1.5, color='black')
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        axis('off')
        title('2091-2100', fontsize=24)
        # Colourbar on the right
        cbaxes = fig.add_axes([0.92, 0.2, 0.01, 0.6])
        cbar = colorbar(img2, cax=cbaxes)
        cbar.ax.tick_params(labelsize=20)
        # Main title
        if var_name == 'melt':
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
        fig.savefig(fig_dir + fig_heads[index] + '_' + var_name + '.png')


# Command-line interface
if __name__ == "__main__":

    var_key = int(input("Ice shelf melt rate (1), bottom water temperature (2), bottom water salinity (3), surface velocity (4), or vertically averaged velocity (5)? "))
    if var_key == 1:
        var_name = 'melt'
    elif var_key == 2:
        var_name = 'temp'
    elif var_key == 3:
        var_name = 'salt'
    elif var_key == 4:
        var_name = 'vsfc'
    elif var_key == 5:
        var_name = 'vavg'
    mesh_path = input("Path to FESOM mesh directory: ")
    file_path = input("Path to output directory for RCP: ")
    fig_dir = input("Path to directory to save figures into: ")
    cavity_fields_rcp(var_name, mesh_path, output_path, fig_dir)
            
                    
            

    
