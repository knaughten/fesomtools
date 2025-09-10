from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from patches import *
from unesco import *

def mld_diff ():

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/high_res/'
    directory_beg = '/short/y99/kaa561/FESOM/highres_spinup/'
    directories = ['/short/y99/kaa561/FESOM/rcp45_M_highres/output/', '/short/y99/kaa561/FESOM/rcp45_A_highres/output/', '/short/y99/kaa561/FESOM/rcp85_M_highres/output/', '/short/y99/kaa561/FESOM/rcp85_A_highres/output/', '/short/y99/kaa561/FESOM/highres_spinup/']
    seasonal_file_beg = 'seasonal_climatology_oce_1996_2005.nc'
    seasonal_file_end = 'seasonal_climatology_oce_2091_2100.nc'
    # Titles for plotting
    expt_names = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A', 'CONTROL']
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    num_expts = len(directories)
    middle_expt = (num_expts+1)//2 - 1
    # Start and end years for each period
    beg_years = [1996, 2005]
    end_years = [2091, 2100]
    # Definition of mixed layer depth: where potential density exceeds
    # surface density by this amount (kg/m^3) as in Sallee et al 2013
    density_anom = 0.03
    # Northern boundary for plot: 63S
    nbdry = -64 + 90
    # Mesh parameters
    circumpolar = True
    mask_cavities = False
    # Maximum for colour scale in each season
    max_bound_summer = 150
    max_bound_winter = 600
    diff_bound_summer = 80
    diff_bound_winter = 1000

    print('Building mesh')
    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)
    print('Processing 1996-2005')
    print('Reading data')
    id = Dataset(directory_beg + seasonal_file_beg, 'r')
    temp_nodes_beg = id.variables['temp'][:,:]
    salt_nodes_beg = id.variables['salt'][:,:]
    id.close()
    print('Calculating density')
    density_nodes_beg = unesco(temp_nodes_beg, salt_nodes_beg, zeros(shape(temp_nodes_beg)))
    print('Calculating mixed layer depth')
    # Set up arrays for mixed layer depth at each element, at each season
    mld_summer_beg = zeros(len(elements))
    mld_winter_beg = zeros(len(elements))
    # Loop over seasons and elements to fill these in
    for season in [0,2]:
        print('...' + season_names[season])
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
        if season == 0:
            mld_summer_beg[:] = mld_season
        elif season == 2:
            mld_winter_beg[:] = mld_season
    # Now calculate anomalies for each experiment
    mld_summer_diff = zeros([num_expts, len(elements)])
    mld_winter_diff = zeros([num_expts, len(elements)])
    for expt in range(num_expts):
        print('Processing ' + expt_names[expt])
        print('Reading data')
        id = Dataset(directories[expt] + seasonal_file_end, 'r')
        temp_nodes_end = id.variables['temp'][:,:]
        salt_nodes_end = id.variables['salt'][:,:]
        id.close()
        print('Calculating density')
        density_nodes_end = unesco(temp_nodes_end, salt_nodes_end, zeros(shape(temp_nodes_end)))
        print('Calculating mixed layer depth')
        for season in [0,2]:
            print('...' + season_names[season])
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
                        if density_nodes_beg[season,curr_node.id] >= density_sfc + density_anom:
                            mld_nodes.append(curr_node.depth-depth_sfc)
                            break
                        temp_depth = curr_node.depth
                        curr_node = curr_node.below
                mld_season.append(mean(array(mld_nodes)))                
            if season == 0:
                mld_summer_diff[expt,:] = array(mld_season) - mld_summer_beg
            elif season == 2:
                mld_winter_diff[expt,:] = array(mld_season) - mld_winter_beg

    print('Plotting')
    fig = figure(figsize=(24,8))
    # Summer, beginning
    ax = fig.add_subplot(2, num_expts+1, 1, aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(mld_summer_beg)
    img.set_clim(vmin=0, vmax=max_bound_summer)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    ax.set_xticks([])
    ax.set_yticks([])
    title(str(beg_years[0])+'-'+str(beg_years[1]), fontsize=20)
    text(-35, 0, season_names[0], fontsize=24)
    # Add a colorbar on the left
    cbaxes = fig.add_axes([0.05, 0.57, 0.02, 0.3])
    cbar = colorbar(img, cax=cbaxes, extend='max', ticks=arange(0, max_bound_summer+50, 50))
    cbar.ax.tick_params(labelsize=16)
    # Summer, anomalies for each experiment
    for expt in range(num_expts):
        ax = fig.add_subplot(2, num_expts+1, expt+2, aspect='equal')
        img = PatchCollection(patches, cmap='RdBu_r')
        img.set_array(mld_summer_diff[expt,:])
        img.set_clim(vmin=-diff_bound_summer, vmax=diff_bound_summer)
        img.set_edgecolor('face')
        ax.add_collection(img)
        xlim([-nbdry, nbdry])
        ylim([-nbdry, nbdry])
        ax.set_xticks([])
        ax.set_yticks([])
        title(expt_names[expt], fontsize=20)
        if expt == num_expts-1:
            # Add a colorbar on the right
            cbaxes = fig.add_axes([0.92, 0.57, 0.02, 0.3])
            cbar = colorbar(img, cax=cbaxes, extend='both', ticks=arange(-diff_bound_summer, diff_bound_summer+40, 40))
            cbar.ax.tick_params(labelsize=16)
    # Winter, beginning
    ax = fig.add_subplot(2, num_expts+1, num_expts+2, aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(mld_winter_beg)
    img.set_clim(vmin=0, vmax=max_bound_winter)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    ax.set_xticks([])
    ax.set_yticks([])
    text(-35, 0, season_names[2], fontsize=24)
    # Add a colorbar on the left
    cbaxes = fig.add_axes([0.05, 0.13, 0.02, 0.3])
    cbar = colorbar(img, cax=cbaxes, extend='max', ticks=arange(0, max_bound_winter+200, 200))
    cbar.ax.tick_params(labelsize=16)
    # Winter, anomalies for each experiment
    for expt in range(num_expts):
        ax =fig.add_subplot(2, num_expts+1, num_expts+expt+3, aspect='equal')
        img = PatchCollection(patches, cmap='RdBu_r')
        img.set_array(mld_winter_diff[expt,:])
        img.set_clim(vmin=-diff_bound_winter, vmax=diff_bound_winter)
        img.set_edgecolor('face')
        ax.add_collection(img)
        xlim([-nbdry, nbdry])
        ylim([-nbdry, nbdry])
        ax.set_xticks([])
        ax.set_yticks([])
        if expt == middle_expt:
            # Add subtitle for anomalies as xlabel
            xlabel(str(end_years[0])+'-'+str(end_years[1])+' anomalies', fontsize=20)
        if expt == num_expts-1:
            # Add a colorbar on the right
            cbaxes = fig.add_axes([0.92, 0.13, 0.02, 0.3])
            cbar = colorbar(img, cax=cbaxes, extend='both', ticks=arange(-diff_bound_winter, diff_bound_winter+500, 500))
            cbar.ax.tick_params(labelsize=16)
    suptitle('Mixed layer depth (m)', fontsize=30)
    subplots_adjust(wspace=0.025, hspace=0.025)
    fig.show()
    fig.savefig('mld_scenarios.png')


# Command-line interface
if __name__ == "__main__":

    mld_diff()
