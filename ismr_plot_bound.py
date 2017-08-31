from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.colors import LinearSegmentedColormap
from patches import *

# Plot the ice shelf melt rate field in the given region.
# Input:
# elements = FESOM grid elements (created by fesom_grid.py)
# mask_patches = plotting patches for all elements (created by make_patches in
#                patches.py)
# patches = plotting patches for ice shelf elements (created by iceshelf_mask in
#           patches.py)
# file_path = path to FESOM forcing.diag.nc file
# tstep = time index in file_path to plot
# bounds = list of size 4 containing minimum longitude, maximum longitude,
#          minimum latitude, maximum latitude to plot. The actual bounds will
#          be larger due to the circumpolar projection and the need for a square
#          plot.
# save = optional boolean indicating to save the figure, rather than display
# fig_name = if save=True, filename for figure
def ismr_plot_bound (elements, mask_patches, patches, file_path, tstep, bounds, save=False, fig_name=None):

    # Constants
    sec_per_year = 365.25*24*3600
    deg2rad = pi/180.0

    # Choose bounds on circumpolar plot
    # User-defined lon-lat bounds
    lon_min = bounds[0]
    lon_max = bounds[1]
    lat_min = bounds[2]
    lat_max = bounds[3]
    # Convert to polar coordinates for plotting
    x1 = -(lat_min+90)*cos(lon_min*deg2rad+pi/2)
    y1 = (lat_min+90)*sin(lon_min*deg2rad+pi/2)
    x2 = -(lat_min+90)*cos(lon_max*deg2rad+pi/2)
    y2 = (lat_min+90)*sin(lon_max*deg2rad+pi/2)
    x3 = -(lat_max+90)*cos(lon_min*deg2rad+pi/2)
    y3 = (lat_max+90)*sin(lon_min*deg2rad+pi/2)
    x4 = -(lat_max+90)*cos(lon_max*deg2rad+pi/2)
    y4 = (lat_max+90)*sin(lon_max*deg2rad+pi/2)
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

    # Read freshwater flux
    file = Dataset(file_path, 'r')
    # Convert from m/s to m/y
    data = file.variables['wnet'][tstep-1,:]*sec_per_year
    file.close()
    values = []
    var_min = 0
    var_max = 0
    # Loop over elements
    for elm in elements:
        # For each element in an ice shelf cavity, append the mean value
        # for the 3 component Nodes
        if elm.cavity:
            values.append(mean([data[elm.nodes[0].id], data[elm.nodes[1].id], data[elm.nodes[2].id]]))
            if any(elm.x >= x_min) and any(elm.x <= x_max) and any(elm.y >= y_min) and any(elm.y <= y_max):
                if values[-1] < var_min:
                    var_min = values[-1]
                if values[-1] > var_max:
                    var_max = values[-1]

    # Set colour map
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

    # Make grey square to fill in the background as land
    x_reg, y_reg = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
    land_square = zeros(shape(x_reg))

    # Plot
    fig = figure(figsize=(16,12))
    ax = fig.add_subplot(1,1,1,aspect='equal')
    # Start with land background
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    # Add ice shelf elements
    img = PatchCollection(patches, cmap=mf_cmap)
    img.set_array(array(values))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    # Mask out the open ocean in white
    overlay = PatchCollection(mask_patches, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)

    # Contour ice shelf fronts
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
    # Add all the lines to the plot
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax.add_collection(contours)

    # Configure plot
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    axis('off')
    title('Ice shelf melt rate (m/y)', fontsize=30)
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=20)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    circumpolar = True
    mask_cavities = True

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    file_path = raw_input("Path to FESOM forcing.diag file: ")
    tstep = int(raw_input("Time index to plot (starting at 1): "))
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

    # Make patches ahead of time
    elements, mask_patches = make_patches(mesh_path, circumpolar, mask_cavities)
    # Mask out open ocean
    patches = iceshelf_mask(elements)
    ismr_plot_bound(elements, mask_patches, patches, file_path, tstep, bounds, save, fig_name)

    # Repeat until the user is finished
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            update_mesh = False
            while True:
                changes = raw_input("Enter a parameter to change: (1) mesh path, (2) file path, (3) timestep, (4) lat/lon bounds, (5) save/display; or enter to continue: ")
                if len(changes) == 0:
                    break
                else:
                    if int(changes) == 1:
                        update_mesh = True
                        mesh_path = raw_input("Path to FESOM mesh directory: ")
                        file_path = raw_input("Path to FESOM forcing.diag file: ")
                    elif int(changes) == 2:
                        file_path = raw_input("Path to FESOM forcing.diag file: ")
                    elif int(changes) == 3:
                        tstep = int(raw_input("Time index to plot (starting at 1): "))
                    elif int(changes) == 4:
                        lon_min = float(raw_input("Minimum longitude (-180 to 180): "))
                        lon_max = float(raw_input("Maximum longitude (-180 to 180): "))
                        lat_min = float(raw_input("Minimum latitude (-90 to -60): "))
                        lat_max = float(raw_input("Maximum latitude (-90 to -60): "))
                        bounds = [lon_min, lon_max, lat_min, lat_max]
                    elif int(changes) == 5:
                        save = not save
            if save:
                fig_name = raw_input("File name for figure: ")
            if update_mesh:
                elements, mask_patches = make_patches(mesh_path, circumpolar, mask_cavities)
                patches = iceshelf_mask(elements)
            ismr_plot_bound(elements, mask_patches, patches, file_path, tstep, bounds, save, fig_name)
        else:
            break
                      
                        

    

    

    
    

    
