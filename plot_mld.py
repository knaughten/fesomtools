from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from patches import *
from unesco import *

# Plot the mixed layer depth (defined as in Sallee et al 2013: depth where
# potential density exceeds surface density by 0.03 kg/m^3).
# Input:
# elements = array of Elements for the FESOM grid (created using fesom_grid)
# patches = array of Polygon patches corresponding to elements (created using
#           make_patches)
# file_path = string containing path to FESOM output file
# tstep = int specifying index of time axis in file_path (1-based)
# circumpolar = boolean flag indicating whether to use the circumpolar domain
#               or the global domain
# save = optional boolean flag indicating whether to save plot to file
#        (otherwise will display on screen)
# fig_name = optional string containing name of figure file, if save = True
# limit = optional float containing upper bound for colour scale
def plot_mld (elements, patches, file_path, tstep, circumpolar, save=False, fig_name=None, limit=None):

    # Definition of mixed layer depth: where potential density exceeds
    # surface density by this amount (kg/m^3) as in Sallee et al 2013
    density_anom = 0.03
    mask_cavities=True

    # Set bounds for domain
    if circumpolar:
        # Northern boundary 30S
        lat_max = -30+90
    else:
        lon_min = -180
        lon_max = 180
        lat_min = -90
        lat_max = 90
        # Configure position of latitude and longitude labels
        lon_ticks = arange(-120, 120+1, 60)
        lat_ticks = arange(-60, 60+1, 30)
    # Set font sizes
    font_sizes = [30, 24, 20]

    # Read temperature and salinity at each node
    print('Reading data')
    id = Dataset(file_path, 'r')
    temp = id.variables['temp'][tstep-1,:]
    salt = id.variables['salt'][tstep-1,:]
    id.close()
    # Calculate potential density (depth 0)
    print('Calculating density')
    density = unesco(temp, salt, zeros(shape(temp)))

    # Build an array of mixed layer depth corresponding to each 2D Element
    print('Calculating mixed layer depth')
    values = []
    for elm in elements:
        if (mask_cavities and not elm.cavity) or (not mask_cavities):
            # Get mixed layer depth at each node
            mld_nodes = []
            # Make sure we exclude ice shelf cavity nodes from element mean
            # (an Element can be  a non-cavity element and still have up to 2
            # cavity nodes)
            for i in range(3):
                if (mask_cavities and not elm.cavity_nodes[i]) or (not mask_cavities):
                    node = elm.nodes[i]
                    density_sfc = density[node.id]
                    temp_depth = node.depth
                    curr_node = node.below
                    while True:
                        if curr_node is None:
                            mld_nodes.append(temp_depth)
                            break
                        if density[curr_node.id] >= density_sfc + density_anom:
                            mld_nodes.append(curr_node.depth)
                            break
                        temp_depth = curr_node.depth
                        curr_node = curr_node.below
            # For this element, save the mean mixed layer depth across
            # non-cavity nodes (up to 3)
            values.append(mean(array(mld_nodes)))

    if mask_cavities:
        # Get mask array of patches for ice shelf cavity elements
        mask_patches = iceshelf_mask(elements)

    # Choose colour bounds
    var_min = 0
    if limit is not None:
        var_max = limit
    else:
        var_max = amax(array(values))

    print('Plotting')
    # Set up plot
    if circumpolar:
        fig = figure(figsize=(16,12))
        ax = fig.add_subplot(1,1,1, aspect='equal')
    else:
        fig = figure(figsize=(16,8))
        ax = fig.add_subplot(1,1,1)
    # Set colourmap for patches, and refer it to the values array
    img = PatchCollection(patches, cmap='jet')
    img.set_array(array(values))
    img.set_edgecolor('face')
    # Add patches to plot
    ax.add_collection(img)
    if mask_cavities:
        # Set colour to light grey for patches in mask
        overlay = PatchCollection(mask_patches, facecolor=(0.6, 0.6, 0.6))
        overlay.set_edgecolor('face')
        # Add mask to plot
        ax.add_collection(overlay)

    # Configure plot
    if circumpolar:
        xlim([-lat_max, lat_max])
        ylim([-lat_max, lat_max])
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        axis('off')
    else:
        xlim([lon_min, lon_max])
        ylim([lat_min, lat_max])
        xticks(lon_ticks)
        yticks(lat_ticks)
        xlabel('Longitude', fontsize=font_sizes[1])
        ylabel('Latitude', fontsize=font_sizes[1])
        setp(ax.get_xticklabels(), fontsize=font_sizes[2])
        setp(ax.get_yticklabels(), fontsize=font_sizes[2])
    title('Mixed layer depth (m)', fontsize=font_sizes[0])
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=font_sizes[2])
    img.set_clim(vmin=var_min, vmax=var_max)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    mask_cavities=True

    mesh_path = input("Path to mesh directory: ")
    file_path = input("Path to FESOM output oce.mean.nc file: ")
    tstep = int(input("Time index to plot (starting at 1): "))
    domain = input("Global (g) or circumpolar (c)? ")
    if domain == 'c':
        circumpolar = True
    elif domain == 'g':
        circumpolar = False
    get_bound = input("Set upper bound on colour scale (y/n)? ")
    if get_bound == 'y':
        limit = float(input("Enter upper bound (positive, in metres): "))
    elif get_bound == 'n':
        limit = None
    action = input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    print("Building grid")
    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)
    plot_mld(elements, patches, file_path, tstep, circumpolar, save, fig_name, limit)

    # Repeat until user wants to exit
    while True:
        repeat = input("Make another plot (y/n)? ")
        if repeat == 'y':
            new_grid = False
            while True:
                changes = input("Enter another parameter to change: (1) mesh path, (2) file path, (3) time index, (4) global/circumpolar, (5) colour bound, (6) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        new_grid = True
                        mesh_path = input("Path to mesh directory: ")
                        # New mesh implies new data file
                        file_path = input("Path to FESOM output oce.mean.nc file: ")
                    elif int(changes) == 2:
                        file_path = input("Path to FESOM output oce.mean.nc file: ")
                    elif int(changes) == 3:
                        tstep = int(input("Time index to plot (starting at 1): "))
                    elif int(changes) == 4:
                        new_grid = True
                        circumpolar = not circumpolar
                    elif int(changes) == 5:
                        get_bound = input("Set upper bound on colour scale (y/n)? ")
                        if get_bound == 'y':
                            limit = float(input("Enter upper bound (positive, in metres): "))
                        else:
                            limit = None
                    elif int(changes) == 6:
                        save = not save
            if save:
                fig_name = input("File name for figure: ")
            if new_grid:
                print("Building grid")
                elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)
            plot_mld (elements, patches, file_path, tstep, circumpolar, save, fig_name, limit)
        else:
            break
            
                    
    
        
