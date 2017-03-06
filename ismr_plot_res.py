from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.colors import LinearSegmentedColormap
from patches import *

def ismr_plot_res (elements_lr, mask_patches_lr, patches_lr, elements_hr, mask_patches_hr, patches_hr, file_lr, file_hr, tstep, bounds, save=False, fig_name=None):

    sec_per_year = 365*24*3600
    deg2rad = pi/180.0

    # Choose bounds on circumpolar plot
    lon_min = bounds[0]
    lon_max = bounds[1]
    lat_min = bounds[2]
    lat_max = bounds[3]
    x1 = -(lat_min+90)*cos(lon_min*deg2rad+pi/2)
    y1 = (lat_min+90)*sin(lon_min*deg2rad+pi/2)
    x2 = -(lat_min+90)*cos(lon_max*deg2rad+pi/2)
    y2 = (lat_min+90)*sin(lon_max*deg2rad+pi/2)
    x3 = -(lat_max+90)*cos(lon_min*deg2rad+pi/2)
    y3 = (lat_max+90)*sin(lon_min*deg2rad+pi/2)
    x4 = -(lat_max+90)*cos(lon_max*deg2rad+pi/2)
    y4 = (lat_max+90)*sin(lon_max*deg2rad+pi/2)
    x_min = amin(array([x1, x2, x3, x4]))
    x_max = amax(array([x1, x2, x3, x4]))
    y_min = amin(array([y1, y2, y3, y4]))
    y_max = amax(array([y1, y2, y3, y4]))
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

    # Read freshwater flux
    file = Dataset(file_lr, 'r')
    # Convert from m/s to m/y
    if tstep == -1:
        data_lr = mean(file.variables['wnet'][:,:], axis=0)*sec_per_year
    else:
        data_lr = file.variables['wnet'][tstep-1,:]*sec_per_year
    file.close()
    values_lr = []
    var_min = 0
    var_max = 0
    # Loop over elements
    for elm in elements_lr:
        # For each element in an ice shelf cavity, append the mean value
        # for the 3 component Nodes
        if elm.cavity:
            values_lr.append(mean([data_lr[elm.nodes[0].id], data_lr[elm.nodes[1].id], data_lr[elm.nodes[2].id]]))
            if any(elm.x >= x_min) and any(elm.x <= x_max) and any(elm.y >= y_min) and any(elm.y <= y_max):
                if values_lr[-1] < var_min:
                    var_min = values_lr[-1]
                if values_lr[-1] > var_max:
                    var_max = values_lr[-1]
    # Repeat for high res
    file = Dataset(file_hr, 'r')
    if tstep == -1:
        data_hr = mean(file.variables['wnet'][:,:], axis=0)*sec_per_year
    else:
        data_hr = file.variables['wnet'][tstep-1,:]*sec_per_year
    file.close()
    values_hr = []
    for elm in elements_hr:
        if elm.cavity:
            values_hr.append(mean([data_hr[elm.nodes[0].id], data_hr[elm.nodes[1].id], data_hr[elm.nodes[2].id]]))
            if any(elm.x >= x_min) and any(elm.x <= x_max) and any(elm.y >= y_min) and any(elm.y <= y_max):
                if values_hr[-1] < var_min:
                    var_min = values_hr[-1]
                if values_hr[-1] > var_max:
                    var_max = values_hr[-1]

    # Set colour map
    if var_min < 0:
        cmap_vals = array([var_min, 0, 0.25*var_max, 0.5*var_max, 0.75*var_max, var_max])
        cmap_colors = [(0.26, 0.45, 0.86), (1, 1, 1), (1, 0.9, 0.4), (0.99, 0.59, 0.18), (0.5, 0.0, 0.08), (0.96, 0.17, 0.89)]
        cmap_vals_norm = (cmap_vals - var_min)/(var_max - var_min)
        cmap_list = []
        for i in range(size(cmap_vals)):
            cmap_list.append((cmap_vals_norm[i], cmap_colors[i]))
        mf_cmap = LinearSegmentedColormap.from_list('melt_freeze', cmap_list)
    else:
        cmap_vals = array([0, 0.25*var_max, 0.5*var_max, 0.75*var_max, var_max])
        cmap_colors = [(1, 1, 1), (1, 0.9, 0.4), (0.99, 0.59, 0.18), (0.5, 0.0, 0.08), (0.96, 0.17, 0.89)]
        cmap_vals_norm = cmap_vals/var_max
        cmap_list = []
        for i in range(size(cmap_vals)):
            cmap_list.append((cmap_vals_norm[i], cmap_colors[i]))
        mf_cmap = LinearSegmentedColormap.from_list('melt_freeze', cmap_list)

    x_reg, y_reg = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
    land_square = zeros(shape(x_reg))

    fig = figure(figsize=(20,10))
    ax1 = fig.add_subplot(1,2,1,aspect='equal')
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img1 = PatchCollection(patches_lr, cmap=mf_cmap)
    img1.set_array(array(values_lr))
    img1.set_edgecolor('face')
    img1.set_clim(vmin=var_min, vmax=var_max)
    ax1.add_collection(img1)
    # Mask out the open ocean in white
    overlay = PatchCollection(mask_patches_lr, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax1.add_collection(overlay)
    # Contour ice shelf fronts
    contour_lines = []
    for elm in elements_lr:
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
    # Add all the lines to the plot
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax1.add_collection(contours)
    # Configure plot
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    axis('off')
    title('Low res', fontsize=20)

    ax2 = fig.add_subplot(1,2,2,aspect='equal')
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img2 = PatchCollection(patches_hr, cmap=mf_cmap)
    img2.set_array(array(values_hr))
    img2.set_edgecolor('face')
    img2.set_clim(vmin=var_min, vmax=var_max)
    ax2.add_collection(img2)
    # Mask out the open ocean in white
    overlay = PatchCollection(mask_patches_hr, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax2.add_collection(overlay)
    # Contour ice shelf fronts
    contour_lines = []
    for elm in elements_hr:
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
    # Add all the lines to the plot
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax2.add_collection(contours)
    # Configure plot
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    axis('off')
    title('High res', fontsize=20)
    cbaxes = fig.add_axes([0.93, 0.2, 0.03, 0.6])
    cbar = colorbar(img2, cax=cbaxes)
    cbar.ax.tick_params(labelsize=20)
    suptitle('Ice shelf melt rate (m/y)', fontsize=30)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


if __name__ == "__main__":

    file_lr = raw_input("Path to FESOM low res forcing.diag file: ")
    file_hr = raw_input("Path to FESOM high res forcing.diag file: ")
    tstep = int(raw_input("Time index to plot (starting at 1) or -1 for annual average: "))
    lon_min = float(raw_input("Minimum longitude (-180 to 180): "))
    lon_max = float(raw_input("Maximum longitude (-180 to 180): "))
    lat_min = float(raw_input("Minimum latitude (-90 to -60): "))
    lat_max = float(raw_input("Maximum latitude (-90 to -60): "))
    bounds = [lon_min, lon_max, lat_min, lat_max]
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save=True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None

    circumpolar = True
    mask_cavities = True
    mesh_path_lr = '../FESOM/mesh/low_res/'
    mesh_path_hr = '../FESOM/mesh/high_res/'
    elements_lr, mask_patches_lr = make_patches(mesh_path_lr, circumpolar, mask_cavities)
    patches_lr = iceshelf_mask(elements_lr)
    elements_hr, mask_patches_hr = make_patches(mesh_path_hr, circumpolar, mask_cavities)
    patches_hr = iceshelf_mask(elements_hr)
    ismr_plot_res(elements_lr, mask_patches_lr, patches_lr, elements_hr, mask_patches_hr, patches_hr, file_lr, file_hr, tstep, bounds, save, fig_name)

    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                changes = raw_input("Enter a parameter to change: (1) file paths, (2) timestep, (3) lat/lon bounds, (4) save/display; or enter to continue: ")
                if len(changes) == 0:
                    break
                else:
                    if int(changes) == 1:
                        file_lr = raw_input("Path to FESOM low res forcing.diag file: ")
                        file_hr = raw_input("Path to FESOM high res forcing.diag file: ")
                    elif int(changes) == 2:
                        tstep = int(raw_input("Time index to plot (starting at 1) or -1 for annual average: "))
                    elif int(changes) == 3:
                        lon_min = float(raw_input("Minimum longitude (-180 to 180): "))
                        lon_max = float(raw_input("Maximum longitude (-180 to 180): "))
                        lat_min = float(raw_input("Minimum latitude (-90 to -60): "))
                        lat_max = float(raw_input("Maximum latitude (-90 to -60): "))
                        bounds = [lon_min, lon_max, lat_min, lat_max]
                    elif int(changes) == 4:
                        save = not save
            if save:
                fig_name = raw_input("File name for figure: ")
            ismr_plot_res(elements_lr, mask_patches_lr, patches_lr, elements_hr, mask_patches_hr, patches_hr, file_lr, file_hr, tstep, bounds, save, fig_name)
        else:
            break
                      
                        

    

    

    
    

    
