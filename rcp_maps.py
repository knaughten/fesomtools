from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from patches import *
from unesco import *

def rcp_maps (var):

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/meshB/'
    directory_beg = '/short/y99/kaa561/FESOM/highres_spinup/'
    directories = ['/short/y99/kaa561/FESOM/rcp45_M/', '/short/y99/kaa561/FESOM/rcp45_A/', '/short/y99/kaa561/FESOM/rcp85_M/', '/short/y99/kaa561/FESOM/rcp85_A/', '/short/y99/kaa561/FESOM/highres_spinup/']
    if var in ['bwtemp', 'bwsalt', 'velsfc']:
        file_beg = 'annual_avg.oce.mean.1996.2005.nc'
        file_end = 'annual_avg.oce.mean.2091.2100.nc'
    elif var in ['sst', 'mld']:
        file_beg = 'seasonal_climatology_oce_1996_2005.nc'
        file_end = 'seasonal_climatology_oce_2091_2100.nc'
    elif var == 'aice':
        file_beg = 'seasonal_climatology_ice_1996_2005.nc'
        file_end = 'seasonal_climatology_ice_2091_2100.nc'
    elif var == 'fwflux':
        file_beg = 'annual_avg.forcing.diag.1996.2005.nc'
        file_end = 'annual_avg.forcing.diag.2091.2100.nc'
    elif var == 'thdgr':
        file_beg = 'annual_avg.ice.diag.1996.2005.nc'
        file_end = 'annual_avg.ice.diag.2091.2100.nc'
    num_expts = len(directories)
    expt_names = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A', 'CONTROL']
    # Northern boundary for plot: 64S
    nbdry = -64 + 90
    # Mesh parameters
    circumpolar = True
    if var in ['aice', 'sst', 'thdgr']:
        mask_cavities = True
    else:
        mask_cavities = False
    # Density anomaly for definition of mixed layer (from Sallee et al 2013)
    density_anom = 0.03
    # Seconds to years conversion factor
    sec_per_year = 365.25*24*60*60
    # Bounds on colour scale, and colourmap for initial plot
    if var == 'bwtemp':
        abs_min = -2
        abs_max = 0.6
        diff_max = 1.8
        abs_cmap = 'jet'
        abs_extend = 'both'
    elif var == 'bwsalt':
        abs_min = 34.3
        abs_max = 34.9
        diff_max = 0.5
        abs_cmap = 'jet'
        abs_extend = 'both'
    elif var == 'velsfc':
        abs_min = 0
        abs_max = 0.2
        diff_max = 0.12
        abs_cmap = 'cool'
        abs_extend = 'max'
    elif var == 'sst':
        abs_min = -2
        abs_max = 5
        diff_max = 4
        abs_cmap = 'jet'
        abs_extend = 'both'
    elif var == 'aice':
        abs_min = 0
        abs_max = 1
        diff_max = 1
        abs_cmap = 'jet'
        abs_extend = 'neither'
    elif var == 'mld':
        abs_min = 0
        abs_max = 600
        diff_max = 400
        abs_cmap = 'jet'
        abs_extend = 'max'
    elif var == 'thdgr':
        abs_min = -15
        abs_max = 15
        diff_max = 3
        abs_cmap = 'PiYG_r'
        abs_extend = 'both'

    print('Building mesh')
    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)
    num_patches = len(patches)
    if var == 'mld':
        # Read a couple of extra things for the mesh
        # Number of 2D nodes
        f = open(mesh_path + 'nod2d.out', 'r')
        n2d = int(f.readline())
        f.close()
        # Lists of which nodes are directly below which
        f = open(mesh_path + 'aux3d.out', 'r')
        max_num_layers = int(f.readline())
        node_columns = zeros([n2d, max_num_layers])
        for n in range(n2d):
            for k in range(max_num_layers):
                node_columns[n,k] = int(f.readline())
        node_columns = node_columns.astype(int)
        f.close()
        # Depth of each 3D node
        f = open(mesh_path + 'nod3d.out', 'r')
        f.readline()
        node_depth = []
        for line in f:
            tmp = line.split()
            node_depth.append(-1*float(tmp[3]))
        f.close()
        node_depth = array(node_depth)

    print('Processing 1996-2005')
    id = Dataset(directory_beg + file_beg, 'r')
    if var in ['bwtemp', 'sst']:
        # bwtemp is annually averaged, sst is DJF. Either way, index 0.
        var_nodes_beg = id.variables['temp'][0,:]
    elif var == 'bwsalt':
        var_nodes_beg = id.variables['salt'][0,:]
    elif var == 'velsfc':
        u_nodes_beg = id.variables['u'][0,:]
        v_nodes_beg = id.variables['v'][0,:]
        var_nodes_beg = sqrt(u_nodes_beg**2 + v_nodes_beg**2)
    elif var == 'thdgr':
        # Convert from m/s to m/y
        var_nodes_beg = id.variables['thdgr'][0,:]*sec_per_year        
    elif var == 'aice':
        # Read DJF
        var_nodes_beg = id.variables['area'][0,:]
    elif var == 'mld':
        # Read JJA temp and salt
        temp_nodes_beg = id.variables['temp'][2,:]
        salt_nodes_beg = id.variables['salt'][2,:]
        density_nodes_beg = unesco(temp_nodes_beg, salt_nodes_beg, zeros(shape(temp_nodes_beg)))
        var_nodes_beg = zeros(n2d)
        for n in range(n2d):
            node_id = node_columns[n,0] - 1
            # Save surface density and depth (might be nonzero because cavities)
            density_sfc = density_nodes_beg[node_id]
            depth_sfc = node_depth[node_id]
            # Now loop down
            for k in range(1, max_num_layers):
                if node_columns[n,k] == -999:
                    # Reached the bottom
                    # Save the last depth
                    var_nodes_beg[n] = node_depth[node_id] - depth_sfc
                    break
                node_id = node_columns[n,k] - 1
                if density_nodes_beg[node_id] >= density_sfc + density_anom:
                    # Reached the critical density anomaly
                    # Save this depth
                    var_nodes_beg[n] = node_depth[node_id] - depth_sfc
                    break
    id.close()

    # Now set up arrays for anomalies
    var_nodes_diff = zeros([num_expts, len(var_nodes_beg)])
    for expt in range(num_expts):
        print('Processing ' + expt_names[expt])
        id = Dataset(directories[expt] + file_end, 'r')
        if var in ['bwtemp', 'sst']:
            var_nodes_diff[expt,:] = id.variables['temp'][0,:] - var_nodes_beg
        elif var == 'bwsalt':
            var_nodes_diff[expt,:] = id.variables['salt'][0,:] - var_nodes_beg
        elif var == 'velsfc':
            u_nodes_end = id.variables['u'][0,:]
            v_nodes_end = id.variables['v'][0,:]
            var_nodes_diff[expt,:] = sqrt(u_nodes_end**2 + v_nodes_end**2) - var_nodes_beg
        elif var == 'thdgr':
            var_nodes_diff[expt,:] = id.variables['thdgr'][0,:]*sec_per_year - var_nodes_beg
        elif var == 'aice':
            var_nodes_diff[expt,:] = id.variables['area'][0,:] - var_nodes_beg
        elif var == 'mld':
            temp_nodes_end = id.variables['temp'][2,:]
            salt_nodes_end = id.variables['salt'][2,:]
            density_nodes_end = unesco(temp_nodes_end, salt_nodes_end, zeros(shape(temp_nodes_end)))
            var_nodes_end = zeros(n2d)
            for n in range(n2d):
                node_id = node_columns[n,0] - 1
                density_sfc = density_nodes_end[node_id]
                depth_sfc = node_depth[node_id]
                for k in range(1, max_num_layers):
                    if node_columns[n,k] == -999:
                        var_nodes_end[n] = node_depth[node_id] - depth_sfc
                        break
                    node_id = node_columns[n,k] - 1
                    if density_nodes_end[node_id] >= density_sfc + density_anom:
                        var_nodes_end[n] = node_depth[node_id] - depth_sfc
                        break
            var_nodes_diff[expt,:] = var_nodes_end[:] - var_nodes_beg
        id.close()

    print('Calculating element-averages')
    var_beg = zeros(num_patches)
    var_diff = zeros([num_expts, num_patches])
    i = 0
    for elm in elements:
        # Skip cavity elements for some variables
        if mask_cavities and elm.cavity:
            continue
        if var in ['bwtemp', 'bwsalt']:
            # Bottom nodes
            var_beg[i] = (var_nodes_beg[elm.nodes[0].find_bottom().id] + var_nodes_beg[elm.nodes[1].find_bottom().id] + var_nodes_beg[elm.nodes[2].find_bottom().id])/3.0
            var_diff[:,i] = (var_nodes_diff[:,elm.nodes[0].find_bottom().id] + var_nodes_diff[:,elm.nodes[1].find_bottom().id] + var_nodes_diff[:,elm.nodes[2].find_bottom().id])/3.0
        else:
            # Surface or 2D nodes
            var_beg[i] = (var_nodes_beg[elm.nodes[0].id] + var_nodes_beg[elm.nodes[1].id] + var_nodes_beg[elm.nodes[2].id])/3.0
            var_diff[:,i] = (var_nodes_diff[:,elm.nodes[0].id] + var_nodes_diff[:,elm.nodes[1].id] + var_nodes_diff[:,elm.nodes[2].id])/3.0
        i += 1

    print('Plotting')
    fig = figure(figsize=(24,5))
    gs = GridSpec(1, num_expts+1)
    gs.update(left=0.07, right=0.93, bottom=0.05, top=0.85, wspace=0.02)
    # Beginning
    ax = subplot(gs[0,0], aspect='equal')
    img = PatchCollection(patches, cmap=abs_cmap)
    img.set_array(var_beg)
    img.set_clim(vmin=abs_min, vmax=abs_max)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    ax.set_xticks([])
    ax.set_yticks([])
    title('1996-2005', fontsize=20)
    # Add a colourbar on the left
    cbaxes = fig.add_axes([0.02, 0.25, 0.015, 0.4])
    cbar = colorbar(img, cax=cbaxes, extend=abs_extend)
    cbar.ax.tick_params(labelsize=16)
    # Anomalies for each experiment
    for expt in range(num_expts):
        ax = subplot(gs[0,expt+1], aspect='equal')
        img = PatchCollection(patches, cmap='RdBu_r')
        img.set_array(var_diff[expt,:])
        img.set_clim(vmin=-diff_max, vmax=diff_max)
        img.set_edgecolor('face')
        ax.add_collection(img)
        xlim([-nbdry, nbdry])
        ylim([-nbdry, nbdry])
        ax.set_xticks([])
        ax.set_yticks([])
        title(expt_names[expt], fontsize=20)
        if expt == 0:
            xlabel('2091-2100 anomalies', fontsize=18)
        if expt == num_expts-1:
            # Add a colourbar on the right
            cbaxes = fig.add_axes([0.95, 0.25, 0.015, 0.4])
            cbar = colorbar(img, cax=cbaxes, extend='both')
            cbar.ax.tick_params(labelsize=16)
    if var == 'bwtemp':
        suptitle(r'Bottom water temperature ($^{\circ}$C)', fontsize=30)
    elif var == 'bwsalt':
        suptitle('Bottom water salinity (psu)', fontsize=30)
    elif var == 'velsfc':
        suptitle('Surface velocity (m/s)', fontsize=30)
    elif var == 'thdgr':
        suptitle('Net sea ice growth rate (m/y)', fontsize=30)
    elif var == 'sst':
        suptitle(r'DJF sea surface temperature ($^{\circ}$C)', fontsize=30)
    elif var == 'aice':
        suptitle('DJF sea ice concentration', fontsize=30)
    elif var == 'mld':
        suptitle('JJA mixed layer depth (m)', fontsize=30)
    fig.show()
    fig.savefig(var + '_maps.png')


# Command-line interface
if __name__ == "__main__":

    key = int(input("Bottom water temperature (1), bottom water salinity (2), surface velocity (3), DJF sea surface temperature (4), DJF sea ice concentration (5), JJA mixed layer depth (6), or sea ice growth rate (7)? "))
    if key == 1:
        var = 'bwtemp'
    elif key == 2:
        var = 'bwsalt'
    elif key == 3:
        var = 'velsfc'
    elif key == 4:
        var = 'sst'
    elif key == 5:
        var = 'aice'
    elif key == 6:
        var = 'mld'
    elif key == 7:
        var = 'thdgr'
    else:
        print('Invalid response')
        exit
    rcp_maps(var)

    
        
        
        
    
                    
        
    
    
