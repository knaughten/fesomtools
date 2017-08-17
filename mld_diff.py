from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from patches import *
from unesco import *

def mld_diff (mesh_path, seasonal_file_beg, seasonal_file_end, expt_name, fig_name):

    # Definition of mixed layer depth: where potential density exceeds
    # surface density by this amount (kg/m^3) as in Sallee et al 2013
    density_anom = 0.03
    # Northern boundary for plot: 63S
    nbdry = -64 + 90
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Mesh parameters
    circumpolar = True
    mask_cavities = False
    # Season names
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    # Maximum for colour scale in each season
    max_bound_summer = 150
    max_bound_winter = 600
    diff_bound_summer = 40
    diff_bound_winter = 500

    print 'Building mesh'
    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)
    print 'Reading data'
    id = Dataset(seasonal_file_beg, 'r')
    temp_nodes_beg = id.variables['temp'][:,:]
    salt_nodes_beg = id.variables['salt'][:,:]
    id.close()
    id = Dataset(seasonal_file_end, 'r')
    temp_nodes_end = id.variables['temp'][:,:]
    salt_nodes_end = id.variables['salt'][:,:]
    id.close()
    print 'Calculating density'
    density_nodes_beg = unesco(temp_nodes_beg, salt_nodes_beg, zeros(shape(temp_nodes_beg)))
    density_nodes_end = unesco(temp_nodes_end, salt_nodes_end, zeros(shape(temp_nodes_end)))
    print 'Calculating mixed layer depth'
    # Set up array for mixed layer depth at each element, at each season
    mld_beg = zeros([4, len(elements)])
    mld_end = zeros([4, len(elements)])
    # Loop over seasons and elements to fill these in
    for season in [0,2]:
        print '...' + season_names[season]
        # Beginning
        mld_season = []
        for elm in elements:
            # Get mixed layer depth at each node
            mld_nodes = []
            for i in range(3):
                node = elm.nodes[i]
                density_sfc = density_nodes_beg[season,node.id]
                # Save surface depth (only nonzero in ice shelf cavities)
                depth_sfc = node.depth
                temp_depth = node.depth
                curr_node = node.below
                while True:
                    if curr_node is None:
                        # Reached the bottom
                        mld_nodes.append(temp_depth-depth_sfc)
                        break
                    if density_nodes_beg[season,curr_node.id] >= density_sfc + density_anom:
                        # Reached the critical density anomaly
                        mld_nodes.append(curr_node.depth-depth_sfc)
                        break
                    temp_depth = curr_node.depth
                    curr_node = curr_node.below
            # For this element, save the mean mixed layer depth
            mld_season.append(mean(array(mld_nodes)))
        mld_beg[season,:] = array(mld_season)
        # End
        mld_season = []
        for elm in elements:
            mld_nodes = []
            for i in range(3):
                node = elm.nodes[i]
                density_sfc = density_nodes_end[season,node.id]
                depth_sfc = node.depth
                temp_depth = node.depth
                curr_node = node.below
                while True:
                    if curr_node is None:
                        mld_nodes.append(temp_depth-depth_sfc)
                        break
                    if density_nodes_end[season,curr_node.id] >= density_sfc + density_anom:
                        mld_nodes.append(curr_node.depth-depth_sfc)
                        break
                    temp_depth = curr_node.depth
                    curr_node = curr_node.below
            mld_season.append(mean(array(mld_nodes)))
        mld_end[season,:] = array(mld_season)
    # Calculate difference
    mld_diff = mld_end - mld_beg

    print 'Plotting'
    fig = figure(figsize=(16,10))
    # Summer, beginning
    ax = fig.add_subplot(2, 3, 1, aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(mld_beg[0,:])
    img.set_clim(vmin=0, vmax=max_bound_summer)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    ax.set_xticks([])
    ax.set_yticks([])
    title('DJF (2006-2015)', fontsize=20)
    # Add a colorbar on the left
    cbaxes = fig.add_axes([0.05, 0.57, 0.02, 0.3])
    cbar = colorbar(img, cax=cbaxes, extend='max', ticks=arange(0, max_bound_summer+50, 50))
    cbar.ax.tick_params(labelsize=20)
    # Summer, end
    ax = fig.add_subplot(2, 3, 2, aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(mld_end[0,:])
    img.set_clim(vmin=0, vmax=max_bound_summer)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    ax.set_xticks([])
    ax.set_yticks([])
    title('DJF (2091-2100)', fontsize=20)
    # Summer, difference
    ax = fig.add_subplot(2, 3, 3, aspect='equal')
    img = PatchCollection(patches, cmap='RdBu_r')
    img.set_array(mld_diff[0,:])
    img.set_clim(vmin=-diff_bound_summer, vmax=diff_bound_summer)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    ax.set_xticks([])
    ax.set_yticks([])
    title('DJF (change)', fontsize=20)
    # Add a colorbar on the right
    cbaxes = fig.add_axes([0.92, 0.57, 0.02, 0.3])
    cbar = colorbar(img, cax=cbaxes, extend='both', ticks=arange(-diff_bound_summer, diff_bound_summer+20, 20))
    cbar.ax.tick_params(labelsize=20)
    # Winter, beginning
    ax = fig.add_subplot(2, 3, 4, aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(mld_beg[2,:])
    img.set_clim(vmin=0, vmax=max_bound_winter)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    ax.set_xticks([])
    ax.set_yticks([])
    title('JJA (2006-2015)', fontsize=20)
    # Add a colorbar on the left
    cbaxes = fig.add_axes([0.05, 0.13, 0.02, 0.3])
    cbar = colorbar(img, cax=cbaxes, extend='max', ticks=arange(0, max_bound_winter+200, 200))
    cbar.ax.tick_params(labelsize=20)
    # Summer, end
    ax = fig.add_subplot(2, 3, 5, aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(mld_end[2,:])
    img.set_clim(vmin=0, vmax=max_bound_winter)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    ax.set_xticks([])
    ax.set_yticks([])
    title('JJA (2091-2100)', fontsize=20)
    # Winter, difference
    ax = fig.add_subplot(2, 3, 6, aspect='equal')
    img = PatchCollection(patches, cmap='RdBu_r')
    img.set_array(mld_diff[2,:])
    img.set_clim(vmin=-diff_bound_winter, vmax=diff_bound_winter)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    ax.set_xticks([])
    ax.set_yticks([])
    title('JJA (change)', fontsize=20)
    # Add a colorbar on the right
    cbaxes = fig.add_axes([0.92, 0.13, 0.02, 0.3])
    cbar = colorbar(img, cax=cbaxes, extend='both', ticks=arange(-diff_bound_winter, diff_bound_winter+250, 250))
    cbar.ax.tick_params(labelsize=20)
    suptitle('Mixed layer depth (m): '+expt_name, fontsize=30)
    fig.savefig(fig_name)


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    seasonal_file_beg = raw_input("Path to seasonal climatology file containing 3D temp and salt for first 10 years: ")
    seasonal_file_end = raw_input("Path to seasonal climatology file containing 3D temp and salt for last 10 years: ")
    expt_name = raw_input("Experiment name for title: ")
    fig_name = raw_input("Filename for figure: ")
    mld_diff(mesh_path, seasonal_file_beg, seasonal_file_end, expt_name, fig_name)
