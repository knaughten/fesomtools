from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from patches import *
from seasonal_avg import *
from unesco import *

# Make a circumpolar Antarctic plot of the change in the winter (JJA average)
# mixed layer depth (defined as in Sallee et al 2013: depth at which potential
# density is 0.03 kg/m^3 higher than at the surface) in the last 10 years of the
# RCP compared to the first 10 years.
# Input:
# mesh_path = path to FESOM mesh directory
# file_path_beg, file_path_end = paths to FESOM oce.mean.nc files containing
#                                5-day climatologies over the first 10 years
#                                and last 10 years of the RCP respectively.
#                                These can be created using average_years.py.
# save = optional boolean indicating to save the figure, rather than display it
#        on the screen
# fig_name = if save=True, filename for figure
# limit = optional upper bound for colour scale; the scale will go from -limit
#         to +limit, centered on 0.
def mld_jja_diff (mesh_path, file_path_beg, file_path_end, save=False, fig_name=None, limit=None):

    # Definition of mixed layer depth: where potential density exceeds
    # surface density by this amount (kg/m^3) as in Sallee et al 2013
    density_anom = 0.03
    # Plotting parameters
    circumpolar=True
    mask_cavities=True
    lat_max = -30+90
    font_sizes = [30, 24, 20]

    print 'Building grid'
    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)

    print 'Reading data'
    # First 10 years
    # Read temperature and salinity at each node, seasonally averaged over JJA
    tmp = seasonal_avg(file_path_beg, file_path_beg, 'temp')
    temp_beg = tmp[2,:]
    tmp = seasonal_avg(file_path_beg, file_path_beg, 'salt')
    salt_beg = tmp[2,:]
    # Last 10 years
    tmp = seasonal_avg(file_path_beg, file_path_end, 'temp')
    temp_end = tmp[2,:]
    tmp = seasonal_avg(file_path_beg, file_path_end, 'salt')
    salt_end = tmp[2,:]
    # Calculate potential density (depth 0)
    print 'Calculating density'
    density_beg = unesco(temp_beg, salt_beg, zeros(shape(temp_beg)))
    density_end = unesco(temp_end, salt_end, zeros(shape(temp_end)))

    # Calculate mixed layer depth at each element
    print 'Calculating mixed layer depth'
    # First 10 years
    mld_beg = []
    for elm in elements:
        if (mask_cavities and not elm.cavity) or (not mask_cavities):
            # Get mixed layer depth at each node
            mld_nodes = []
            # Make sure we exclude ice shelf cavity nodes from element mean
            # (an Element can be a non-cavity element and still have up to 2
            # cavity nodes)
            for i in range(3):
                if (mask_cavities and not elm.cavity_nodes[i]) or (not mask_cavities):
                    node = elm.nodes[i]
                    density_sfc = density_beg[node.id]
                    temp_depth = node.depth
                    curr_node = node.below
                    while True:
                        if curr_node is None:
                            # Reached bottom
                            mld_nodes.append(temp_depth)
                            break
                        if density_beg[curr_node.id] >= density_sfc + density_anom:
                            # Reached critical density anomaly
                            mld_nodes.append(curr_node.depth)
                            break
                        temp_depth = curr_node.depth
                        curr_node = curr_node.below
            # For this element, save the mean mixed layer depth across
            # non-cavity nodes (up to 3)
            mld_beg.append(mean(array(mld_nodes)))
    # Last 10 years
    mld_end = []
    for elm in elements:
        if (mask_cavities and not elm.cavity) or (not mask_cavities):
            # Get mixed layer depth at each node
            mld_nodes = []
            # Make sure we exclude ice shelf cavity nodes from element mean
            # (an Element can be a non-cavity element and still have up to 2
            # cavity nodes)
            for i in range(3):
                if (mask_cavities and not elm.cavity_nodes[i]) or (not mask_cavities):
                    node = elm.nodes[i]
                    density_sfc = density_end[node.id]
                    temp_depth = node.depth
                    curr_node = node.below
                    while True:
                        if curr_node is None:
                            mld_nodes.append(temp_depth)
                            break
                        if density_end[curr_node.id] >= density_sfc + density_anom:
                            mld_nodes.append(curr_node.depth)
                            break
                        temp_depth = curr_node.depth
                        curr_node = curr_node.below
            # For this element, save the mean mixed layer depth across
            # non-cavity nodes (up to 3)
            mld_end.append(mean(array(mld_nodes)))
    # Calculate change in mixed layer depth
    mld_change = array(mld_end) - array(mld_beg)

    if mask_cavities:
        # Get mask array of patches for ice shelf cavity elements
        mask_patches = iceshelf_mask(elements)

    # Choose colour bounds
    if limit is not None:
        bound = limit
    else:
        bound = amax(array(mld_change))

    print 'Plotting'
    # Set up plot
    fig = figure(figsize=(16,12))
    ax = fig.add_subplot(1,1,1, aspect='equal')
    # Set colourmap for patches, and refer it to the values array
    img = PatchCollection(patches, cmap='RdBu_r')
    img.set_array(array(mld_change))
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
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    axis('off')
    title('Change in JJA mixed layer depth (m)\n2091-2100 vs 2006-2015', fontsize=font_sizes[0])
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=font_sizes[2])
    img.set_clim(vmin=-bound, vmax=bound)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()
    

# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    file_path_beg = raw_input("Path to FESOM oce.mean.nc file with climatology for first 10 years of RCP: ")
    file_path_end = raw_input("Path to FESOM oce.mean.nc file with climatology for last 10 years of RCP: ")
    get_bound = raw_input("Set upper bound on colour scale (y/n)? ")
    if get_bound == 'y':
        limit = float(raw_input("Enter upper bound (positive, in metres): "))
    elif get_bound == 'n':
        limit = None
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    mld_jja_diff(mesh_path, file_path_beg, file_path_end, save, fig_name, limit)
