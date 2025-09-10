from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.colors import LinearSegmentedColormap
from patches import *

# Make a circumpolar plot of ice shelf melt rates annually averaged over the
# given year of FESOM simulation.
# Input:
# mesh_path = path to FESOM mesh directory
# file_path = path to output file containing one year of forcing averages
# save = optional boolean indicating to save the figure, rather than display
#        it on the screen
# fig_name = if save=True, filename for figure
def ismr_plot (mesh_path, file_path, save=False, fig_name=None):

    # Plotting parameters
    lat_max = -63 + 90
    circumpolar = True
    mask_cavities = True
    # Seconds per year
    sec_per_year = 365.25*24*3600

    # Set colour map
    # Values for change points
    cmap_vals = array([-0.1, 0, 1, 2, 5, 8])
    # Colours for change points
    # (blue, white, yellow-orange, red-orange, dark red, purple)
    cmap_colors = [(0.26, 0.45, 0.86), (1, 1, 1), (1, 0.9, 0.4), (0.99, 0.59, 0.18), (0.5, 0.0, 0.08), (0.96, 0.17, 0.89)]
    # Map to 0-1
    cmap_vals_norm = (cmap_vals + 0.1)/(8 + 0.1)
    # Combine into a list
    cmap_list = []
    for i in range(size(cmap_vals)):
        cmap_list.append((cmap_vals_norm[i], cmap_colors[i]))
    # Make colour map    
    mf_cmap = LinearSegmentedColormap.from_list('melt_freeze', cmap_list)

    # Build FESOM mesh
    # Get separate patches for the open ocean elements so we can mask them out
    elements, mask_patches = make_patches(mesh_path, circumpolar, mask_cavities)
    patches = iceshelf_mask(elements)

    # Read freshwater flux
    file = Dataset(file_path, 'r')
    # Annually average and convert from m/s to m/y
    data = mean(file.variables['wnet'][:,:]*sec_per_year, axis=0)
    file.close()
    values = []
    # Loop over elements
    for elm in elements:
        # For each element in an ice shelf cavity, append the mean value
        # for the 3 component Nodes
        if elm.cavity:
            values.append(mean([data[elm.nodes[0].id], data[elm.nodes[1].id], data[elm.nodes[2].id]]))

    # Set up a grey square covering the domain, anything that isn't covered
    # up later is land
    x_reg, y_reg = meshgrid(linspace(-lat_max, lat_max, num=100), linspace(-lat_max, lat_max, num=100))
    land_square = zeros(shape(x_reg))

    # Plot
    fig = figure(figsize=(16,12))
    ax = fig.add_subplot(1,1,1,aspect='equal')
    # Start with grey square background for land
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches, cmap=mf_cmap)
    img.set_array(array(values))
    img.set_edgecolor('face')
    img.set_clim(vmin=-0.1, vmax=8)
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
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    axis('off')
    title('Ice shelf melt rate (m/y), annual average', fontsize=30)
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=20)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    mesh_path = input("Path to FESOM mesh directory: ")
    file_path = input("Path to FESOM forcing.diag.nc containing one year of output: ")
    action = input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save=True
        fig_name = input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    ismr_plot(mesh_path, file_path, save, fig_name)

    
