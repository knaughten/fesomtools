from numpy import *
from netCDF4 import Dataset
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.pyplot import *
from patches import *

def cavity_warming_depth (rcp, model):

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/meshB/'
    directory_beg = '/short/y99/kaa561/FESOM/highres_spinup/'
    directory_end = '/short/y99/kaa561/FESOM/rcp' + rcp + '_' + model + '/'
    annual_file_beg = 'annual_avg.oce.mean.1996.2005.nc'
    annual_file_end = 'annual_avg.oce.mean.2091.2100.nc'
    monthly_file_beg = 'monthly_climatology_temp_1996_2005.nc'
    monthly_file_end = 'monthly_climatology_temp_2091_2100.nc'
    # Bound on plot
    lat_max = -64+90

    print 'Building mesh'
    # Mask open ocean
    elements, mask_patches = make_patches(mesh_path, circumpolar=True, mask_cavities=True)
    # Unmask ice shelf patches
    patches = iceshelf_mask(elements)
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
    # Set up a grey square covering the domain, anything that isn't covered
    # up later is land
    x_reg, y_reg = meshgrid(linspace(-lat_max, lat_max, num=100), linspace(-lat_max, lat_max, num=100))
    land_square = zeros(shape(x_reg))
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
    # Read the depth of each 3D node
    f = open(mesh_path + 'nod3d.out', 'r')
    f.readline()
    node_depth = []
    for line in f:
        tmp = line.split()
        node_depth.append(-1*float(tmp[3]))
    f.close()
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

    print 'Reading data'
    # Annual average
    id = Dataset(directory_beg + annual_file_beg, 'r')
    temp_beg = id.variables['temp'][0,:]
    id.close()
    id = Dataset(directory_end + annual_file_end, 'r')
    temp_end = id.variables['temp'][0,:]
    id.close()
    # Monthly climatology
    id = Dataset(directory_beg + monthly_file_beg, 'r')
    monthly_temp_beg = id.variables['temp'][:,:]
    id.close()
    id = Dataset(directory_end + monthly_file_end, 'r')
    monthly_temp_end = id.variables['temp'][:,:]
    id.close()

    print 'Processing nodes'
    max_warming_nodes = zeros(n2d)
    fractional_depth_nodes = zeros(n2d)
    seasonality_nodes = zeros(n2d)
    for n in range(n2d):
        if node_cavity[n]:
            # Calculate the warming at each depth down the water column
            # Also save the depths
            warming_column = []
            depth_column = []
            for k in range(max_num_layers):
                if node_columns[n,k] == -999:
                    # Reached the bottom
                    break
                node_id = node_columns[n,k] - 1
                warming_column.append(temp_end[node_id] - temp_beg[node_id])
                depth_column.append(node_depth[node_id])
            warming_column = array(warming_column)
            depth_column = array(depth_column)
            # Save surface depth and water column thickness
            sfc_depth = depth_column[0]
            wct = depth_column[-1] - depth_column[0]
            # Select index of maximum warming
            k0 = argmax(warming_column)
            node_id = node_columns[n,k0] - 1
            # Save maximum warming
            max_warming_nodes[n] = warming_column[k0]
            # Save fractional depth in cavity
            # 0 is ice shelf base, 1 is bottom
            fractional_depth_nodes[n] = (depth_column[k0] - sfc_depth)/wct
            # For this node, get warming for every month in monthly climatology
            monthly_warming = monthly_temp_end[:,node_id] - monthly_temp_beg[:,node_id]
            # Calculate seasonality metric
            # Difference between max and min warming over months, divided
            # by annual warming. 0 means no seasonality.
            seasonality_nodes[n] = (amax(monthly_warming) - amin(monthly_warming))/max_warming_nodes[n]

    print 'Calculating element averages'
    max_warming = []
    fractional_depth = []
    seasonality = []
    cooling_patches = []
    i = 0
    for elm in elements:
        if elm.cavity:
            if any(array([max_warming_nodes[elm.nodes[0].id], max_warming_nodes[elm.nodes[1].id], max_warming_nodes[elm.nodes[2].id]]) < 0):
                max_warming.append(NaN)
                fractional_depth.append(NaN)
                seasonality.append(NaN)
                cooling_patches.append(patches[i])
            else:
                max_warming.append(mean(array([max_warming_nodes[elm.nodes[0].id], max_warming_nodes[elm.nodes[1].id], max_warming_nodes[elm.nodes[2].id]])))
                fractional_depth.append(mean(array([fractional_depth_nodes[elm.nodes[0].id], fractional_depth_nodes[elm.nodes[1].id], fractional_depth_nodes[elm.nodes[2].id]])))
                seasonality.append(mean(array([seasonality_nodes[elm.nodes[0].id], seasonality_nodes[elm.nodes[1].id], seasonality_nodes[elm.nodes[2].id]])))
            i += 1
    max_warming = array(max_warming)
    fractional_depth = array(fractional_depth)
    seasonality = array(seasonality)
    # Mask out regions which cool
    max_warming = ma.masked_where(isnan(max_warming), max_warming)
    fractional_depth = ma.masked_where(isnan(fractional_depth), fractional_depth)
    seasonality = ma.masked_where(isnan(seasonality), seasonality)

    print 'Plotting'
    fig = figure(figsize=(8,14))
    fig.patch.set_facecolor('white')
    gs = GridSpec(3,1)
    gs.update(left=0, right=1, bottom=0, top=0.9, hspace=0)    
    # Maximum warming
    ax = subplot(gs[0,0], aspect='equal')
    # Start with grey square background for land
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    # Mask out cooling patches in white
    cooling = PatchCollection(cooling_patches, facecolor=(1,1,1))
    cooling.set_edgecolor('face')
    ax.add_collection(cooling)
    img = PatchCollection(patches, cmap='jet')
    img.set_array(max_warming)
    img.set_clim(vmin=0, vmax=1.5)
    img.set_edgecolor('face')
    ax.add_collection(img)
    # Mask out the open ocean in white
    overlay = PatchCollection(mask_patches, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    # Contour ice shelf fronts
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax.add_collection(contours)
    # Configure plot
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    axis('off')
    title(r'a) Maximum warming over depth ($^{\circ}$C)', fontsize=22)
    cbaxes_a = fig.add_axes([0.8, 0.67, 0.03, 0.16])
    cbar = colorbar(img, cax=cbaxes_a, extend='max', ticks=arange(0, 1.5+0.5, 0.5))
    cbar.ax.tick_params(labelsize=14)
    # Fractional depth
    ax = subplot(gs[1,0], aspect='equal')
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    cooling = PatchCollection(cooling_patches, facecolor=(1,1,1))
    cooling.set_edgecolor('face')
    ax.add_collection(cooling)
    img = PatchCollection(patches, cmap='jet')
    img.set_array(fractional_depth)
    img.set_edgecolor('face')
    ax.add_collection(img)
    overlay = PatchCollection(mask_patches, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax.add_collection(contours)
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    axis('off')
    title('b) Fractional depth of maximum warming', fontsize=22)
    cbaxes_b = fig.add_axes([0.8, 0.37, 0.03, 0.16])
    cbar = colorbar(img, cax=cbaxes_b, ticks=arange(0, 1+0.25, 0.25))
    cbar.ax.tick_params(labelsize=14)
    # Seasonality
    ax = subplot(gs[2,0], aspect='equal')
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches, cmap='jet')
    cooling = PatchCollection(cooling_patches, facecolor=(1,1,1))
    cooling.set_edgecolor('face')
    ax.add_collection(cooling)
    img.set_array(seasonality)
    img.set_edgecolor('face')
    img.set_clim(vmin=0, vmax=3)
    ax.add_collection(img)
    overlay = PatchCollection(mask_patches, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax.add_collection(contours)
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    axis('off')
    title('c) Seasonality metric', fontsize=22)
    cbaxes_c = fig.add_axes([0.8, 0.07, 0.03, 0.16])
    cbar = colorbar(img, cax=cbaxes_c, extend='max', ticks=arange(0, 3+1, 1))
    cbar.ax.tick_params(labelsize=14)
    suptitle('RCP ' + rcp[0] + '.' + rcp[1] + ' ' + model + ', 2091-2100 minus 1996-2005', fontsize=26)
    fig.show()
    fig.savefig('warming_depth_rcp'+rcp+'_'+model+'.png')


# Command-line interface
if __name__ == "__main__":

    key = int(raw_input('RCP 4.5 (4) or 8.5 (8)? '))
    if key == 4:
        rcp = '45'
    elif key == 8:
        rcp = '85'
    key = int(raw_input('Multi-model mean (1) or ACCESS 1.0 (2)? '))
    if key == 1:
        model = 'M'
    elif key == 2:
        model = 'A'
    cavity_warming_depth(rcp, model)
