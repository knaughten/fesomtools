from netCDF4 import Dataset, num2date
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from patches import *
from monthly_avg import *

def rcp_aice_minmax ():

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/high_res/'
    directory_beg = '/short/y99/kaa561/FESOM/highres_spinup/'
    directories = ['/short/y99/kaa561/FESOM/rcp45_M_highres/output/', '/short/y99/kaa561/FESOM/rcp45_A_highres/output/', '/short/y99/kaa561/FESOM/rcp85_M_highres/output/', '/short/y99/kaa561/FESOM/rcp85_A_highres/output/', '/short/y99/kaa561/FESOM/highres_spinup/']
    file_beg = 'avg.ice.mean.1996.2005.nc'
    file_end = 'avg.ice.mean.2091.2100.nc'
    # Titles for plotting
    expt_names = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A', 'CONTROL']
    num_expts = len(directories)
    # Start and end years for each period
    beg_years = [1996, 2005]
    end_years = [2091, 2100]
    # FESOM plotting parameters
    circumpolar = True
    mask_cavities = True
    # Boundaries on plot (under polar coordinate transformation)
    x_min = -36.25
    x_max = 36.25
    y_min = -34.5
    y_max = 38

    print 'Building mesh'
    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)

    print 'Reading data'
    # Get averages for February and September
    print '...1996-2005'
    fesom_feb_nodes_beg = monthly_avg(directory_beg + file_beg, 'area', 1)
    fesom_sep_nodes_beg = monthly_avg(directory_beg + file_beg, 'area', 8)
    # Calculate anomalies for the rest of the experiments
    fesom_feb_nodes_diff = empty([num_expts, size(fesom_feb_nodes_beg)])
    fesom_sep_nodes_diff = empty([num_expts, size(fesom_sep_nodes_beg)])
    for expt in range(num_expts):
        print '...' + expt_names[expt]
        fesom_feb_nodes_diff[expt,:] = monthly_avg(directories[expt] + file_end, 'area', 1) - fesom_feb_nodes_beg
        fesom_sep_nodes_diff[expt,:] = monthly_avg(directories[expt] + file_end, 'area', 8) - fesom_sep_nodes_beg
    # Find element-averages
    fesom_feb_beg = empty(len(patches))
    fesom_sep_beg = empty(len(patches))
    fesom_feb_diff = empty([num_expts, len(patches)])
    fesom_sep_diff = empty([num_expts, len(patches)])
    i = 0
    for elm in elements:
        if not elm.cavity:
            fesom_feb_beg[i] = mean(array([fesom_feb_nodes_beg[elm.nodes[0].id], fesom_feb_nodes_beg[elm.nodes[1].id], fesom_feb_nodes_beg[elm.nodes[2].id]]))
            fesom_sep_beg[i] = mean(array([fesom_sep_nodes_beg[elm.nodes[0].id], fesom_sep_nodes_beg[elm.nodes[1].id], fesom_sep_nodes_beg[elm.nodes[2].id]]))
            for expt in range(num_expts):
                fesom_feb_diff[expt,i] = mean(array([fesom_feb_nodes_diff[expt,elm.nodes[0].id], fesom_feb_nodes_diff[expt,elm.nodes[1].id], fesom_feb_nodes_diff[expt,elm.nodes[2].id]]))
                fesom_sep_diff[expt,i] = mean(array([fesom_sep_nodes_diff[expt,elm.nodes[0].id], fesom_sep_nodes_diff[expt,elm.nodes[1].id], fesom_sep_nodes_diff[expt,elm.nodes[2].id]]))
            i += 1

    print 'Plotting'
    fig = figure(figsize=(24,8))
    # Feburary
    # 1996-2005
    ax = fig.add_subplot(2, num_expts+1, 1, aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(fesom_feb_beg)
    img.set_clim(vmin=0, vmax=1)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    title(str(beg_years[0])+'-'+str(beg_years[1]), fontsize=20)
    text(-39, 0, 'February', ha='center', va='center', rotation=90, fontsize=20)
    # Colourbar on the left
    cbaxes = fig.add_axes([0.06, 0.4, 0.015, 0.3])
    cbar = colorbar(img, cax=cbaxes, ticks=arange(0, 1+0.25, 0.25))
    cbar.ax.tick_params(labelsize=16)
    # Loop over the rest of the experiments    
    for expt in range(num_expts):
        ax = fig.add_subplot(2, num_expts+1, expt+2, aspect='equal')
        img = PatchCollection(patches, cmap='RdBu_r')
        img.set_array(fesom_feb_diff[expt,:])
        img.set_clim(vmin=-1, vmax=1)
        img.set_edgecolor('face')
        ax.add_collection(img)
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        ax.set_xticks([])
        ax.set_yticks([])
        title(expt_names[expt], fontsize=20)
        if expt == num_expts-1:
            # Colourbar on the right
            cbaxes = fig.add_axes([0.92, 0.4, 0.015, 0.3])
            cbar = colorbar(img, cax=cbaxes, ticks=arange(-1, 1+0.5, 0.5))
            cbar.ax.tick_params(labelsize=16)
    # September
    # 1996-2005
    ax = fig.add_subplot(2, num_expts+1, num_expts+2, aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(fesom_sep_beg)
    img.set_clim(vmin=0, vmax=1)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    text(-39, 0, 'September', ha='center', va='center', rotation=90, fontsize=20)
    # Loop over the rest of the experiments    
    for expt in range(num_expts):
        ax = fig.add_subplot(2, num_expts+1, num_expts+expt+3, aspect='equal')
        img = PatchCollection(patches, cmap='RdBu_r')
        img.set_array(fesom_sep_diff[expt,:])
        img.set_clim(vmin=-1, vmax=1)
        img.set_edgecolor('face')
        ax.add_collection(img)
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        ax.set_xticks([])
        ax.set_yticks([])
        if expt == num_expts-1:
            xlabel(str(end_years[0])+'-'+str(end_years[1])+' anomalies', fontsize=20)
    suptitle('Sea ice concentration', fontsize=30)
    subplots_adjust(wspace=0.025, hspace=0.025)
    fig.show()
    fig.savefig('rcp_aice_minmax.png')


# Command-line interface
if __name__ == "__main__":

    rcp_aice_minmax()
        
    
