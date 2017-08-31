from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from patches import *

def bw_ts_rcp ():

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/high_res/'
    directory_beg = '/short/y99/kaa561/FESOM/highres_spinup/'
    directories = ['/short/y99/kaa561/FESOM/rcp45_M_highres/output/', '/short/y99/kaa561/FESOM/rcp45_A_highres/output/', '/short/y99/kaa561/FESOM/rcp85_M_highres/output/', '/short/y99/kaa561/FESOM/rcp85_A_highres/output/', '/short/y99/kaa561/FESOM/highres_spinup/']
    file_beg = 'annual_avg.oce.mean.1996.2005.nc'
    file_end = 'annual_avg.oce.mean.2091.2100.nc'
    # Titles for plotting
    expt_names = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A', 'CONTROL']
    num_expts = len(directories)
    middle_expt = (num_expts+1)/2 - 1
    # Start and end years for each period
    beg_years = [1996, 2005]
    end_years = [2091, 2100]
    # Northern boundary for plot: 64S
    nbdry = -64 + 90
    # Mesh parameters
    circumpolar = True
    mask_cavities = False
    # Bounds on variables
    temp_min = -2.5
    temp_max = 0.5
    salt_min = 34
    salt_max = 35
    temp_max_diff = 2
    salt_max_diff = 0.5

    print 'Building mesh'
    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)
    print 'Processing 1996-2005'
    id = Dataset(directory_beg + file_beg, 'r')
    temp_nodes_beg = id.variables['temp'][0,:]
    salt_nodes_beg = id.variables['salt'][0,:]
    id.close()
    # Get temperature and salinity averaged over three corners of each bottom
    # Element
    temp_beg = []
    salt_beg = []
    for elm in elements:
        temp_beg.append(mean([temp_nodes_beg[elm.nodes[0].find_bottom().id], temp_nodes_beg[elm.nodes[1].find_bottom().id], temp_nodes_beg[elm.nodes[2].find_bottom().id]]))
        salt_beg.append(mean([salt_nodes_beg[elm.nodes[0].find_bottom().id], salt_nodes_beg[elm.nodes[1].find_bottom().id], salt_nodes_beg[elm.nodes[2].find_bottom().id]]))
    temp_beg = array(temp_beg)
    salt_beg = array(salt_beg)
    # Now calculate anomalies for each experiment
    temp_diff = zeros([num_expts, len(elements)])
    salt_diff = zeros([num_expts, len(elements)])
    for expt in range(num_expts):
        print 'Processing ' + expt_names[expt]
        id = Dataset(directories[expt] + file_end, 'r')
        temp_nodes_diff = id.variables['temp'][0,:] - temp_nodes_beg
        salt_nodes_diff = id.variables['salt'][0,:] - salt_nodes_beg
        id.close()
        i = 0
        for elm in elements:
            temp_diff[expt,i] = mean([temp_nodes_diff[elm.nodes[0].find_bottom().id], temp_nodes_diff[elm.nodes[1].find_bottom().id], temp_nodes_diff[elm.nodes[2].find_bottom().id]])
            salt_diff[expt,i] = mean([salt_nodes_diff[elm.nodes[0].find_bottom().id], salt_nodes_diff[elm.nodes[1].find_bottom().id], salt_nodes_diff[elm.nodes[2].find_bottom().id]])
            i += 1

    print 'Plotting'
    fig = figure(figsize=(24,8))
    # Temperature, beginning
    ax = fig.add_subplot(2, num_expts+1, 1, aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(temp_beg)
    img.set_clim(vmin=temp_min, vmax=temp_max)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    ax.set_xticks([])
    ax.set_yticks([])
    title(str(beg_years[0])+'-'+str(beg_years[1]), fontsize=20)
    text(-35, 20, r'Temperature ($^{\circ}$C)', fontsize=20, rotation=90)
    # Add a colorbar on the left
    cbaxes = fig.add_axes([0.05, 0.57, 0.015, 0.3])
    cbar = colorbar(img, cax=cbaxes, extend='both')
    cbar.ax.tick_params(labelsize=16)
    # Temperature, anomalies for each experiment
    for expt in range(num_expts):
        ax = fig.add_subplot(2, num_expts+1, expt+2, aspect='equal')
        img = PatchCollection(patches, cmap='RdBu_r')
        img.set_array(temp_diff[expt,:])
        img.set_clim(vmin=-temp_max_diff, vmax=temp_max_diff)
        img.set_edgecolor('face')
        ax.add_collection(img)
        xlim([-nbdry, nbdry])
        ylim([-nbdry, nbdry])
        ax.set_xticks([])
        ax.set_yticks([])
        title(expt_names[expt], fontsize=20)
        if expt == num_expts-1:
            # Add a colorbar on the right
            cbaxes = fig.add_axes([0.92, 0.57, 0.015, 0.3])
            cbar = colorbar(img, cax=cbaxes, extend='both')
            cbar.ax.tick_params(labelsize=16)
    # Salinity, beginning
    ax = fig.add_subplot(2, num_expts+1, num_expts+2, aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(salt_beg)
    img.set_clim(vmin=salt_min, vmax=salt_max)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    ax.set_xticks([])
    ax.set_yticks([])
    text(-35, 10, 'Salinity (psu)', fontsize=20, rotation=90)
    # Add a colorbar on the left
    cbaxes = fig.add_axes([0.05, 0.13, 0.015, 0.3])
    cbar = colorbar(img, cax=cbaxes, extend='both')
    cbar.ax.tick_params(labelsize=16)
    # Salinity, anomalies for each experiment
    for expt in range(num_expts):
        ax = fig.add_subplot(2, num_expts+1, num_expts+expt+3, aspect='equal')
        img = PatchCollection(patches, cmap='RdBu_r')
        img.set_array(salt_diff[expt,:])
        img.set_clim(vmin=-salt_max_diff, vmax=salt_max_diff)
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
            cbaxes = fig.add_axes([0.92, 0.13, 0.015, 0.3])
            cbar = colorbar(img, cax=cbaxes, extend='both')
            cbar.ax.tick_params(labelsize=16)
    subplots_adjust(wspace=0.025, hspace=0.025)
    fig.show()
    fig.savefig('bw_ts_rcp.png')


# Command-line interface
if __name__ == "__main__":

    bw_ts_rcp()
    
