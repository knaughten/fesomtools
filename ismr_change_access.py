from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.colors import LinearSegmentedColormap
from patches import *

def ismr_change_access():

    directory_head = '/short/y99/kaa561/FESOM/'
    mesh_path = directory_head + 'mesh/low_res/'
    control_file = directory_head + 'lowres_spinup/rep3/avg.forcing.diag.nc'
    access_file = directory_head + 'rcp45_A/output/avg.forcing.diag.2006.2015.nc'

    # Plotting parameters
    lat_max = -63 + 90
    circumpolar = True
    mask_cavities = True
    # Seconds per year
    sec_per_year = 365*24*3600

    # Build FESOM mesh
    # Get separate patches for the open ocean elements so we can mask them out
    elements, mask_patches = make_patches(mesh_path, circumpolar, mask_cavities)
    patches = iceshelf_mask(elements)

    # Read freshwater flux in control simulation (1992-2005 climatology)
    file = Dataset(control_file, 'r')
    # Annually average and convert from m/s to m/y
    data1 = mean(file.variables['wnet'][:,:]*sec_per_year, axis=0)
    file.close()
    # Same for RCP 4.5 A simulation (2006-2015)
    file = Dataset(access_file, 'r')
    data2 = mean(file.variables['wnet'][:,:]*sec_per_year, axis=0)
    file.close()
    # Get difference
    data = data2 - data1
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

    max_val = 3.0

    # Plot
    fig = figure(figsize=(16,12))
    ax = fig.add_subplot(1,1,1,aspect='equal')
    # Start with grey square background for land
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches, cmap='YlOrRd')
    img.set_array(array(values))
    img.set_edgecolor('face')
    img.set_clim(0, max_val)
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
    title('Change in ice shelf melt rate (m/y)\nRCP 4.5 A (2006-2015) vs control (1992-2005)', fontsize=30)
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=20)

    fig.show()
    #fig.savefig('ismr_change_access.png')


# Command-line interface
if __name__ == "__main__":

    ismr_change_access()

    
