from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.cm import *
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import *
from patches import *
from monthly_avg import *

def rcp_seaice ():

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/meshB/'
    directory_beg = '/short/y99/kaa561/FESOM/highres_spinup/'
    directories = ['/short/y99/kaa561/FESOM/rcp45_M/', '/short/y99/kaa561/FESOM/rcp45_A/', '/short/y99/kaa561/FESOM/rcp85_M/', '/short/y99/kaa561/FESOM/rcp85_A/', '/short/y99/kaa561/FESOM/highres_spinup/']
    file_beg = 'avg.ice.mean.1996.2005.nc'
    file_end = 'avg.ice.mean.2091.2100.nc'
    # Titles for plotting
    expt_names = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A', 'CONTROL']
    num_expts = len(directories)
    # FESOM plotting parameters
    circumpolar = True
    mask_cavities = True
    # Boundaries on plot (under polar coordinate transformation)
    x_min = -36.25
    x_max = 36.25
    y_min = -34.5
    y_max = 38
    # Locations to plot each experiment
    j_plot = [0,0,1,1,1]
    i_plot = [1,2,1,2,0]

    print 'Building mesh'
    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)
    num_ocn_elm = len(patches)

    print 'Reading data'
    print '...1996-2005'
    aice_nodes_beg = monthly_avg(directory_beg + file_beg, 'area', 8)
    n2d = size(aice_nodes_beg)
    # Anomalies for the rest of the experiments
    aice_nodes_diff = empty([num_expts, n2d])
    for expt in range(num_expts):
        print '...' + expt_names[expt]
        aice_nodes_diff[expt,:] = monthly_avg(directories[expt] + file_end, 'area', 8) - aice_nodes_beg

    print 'Calculating element-averages'
    aice_beg = empty(num_ocn_elm)
    aice_diff = empty([num_expts, num_ocn_elm])
    i = 0
    for elm in elements:
        if not elm.cavity:
            aice_beg[i] = mean(array([aice_nodes_beg[elm.nodes[0].id], aice_nodes_beg[elm.nodes[1].id], aice_nodes_beg[elm.nodes[2].id]]))
            for expt in range(num_expts):
                aice_diff[expt,i] = mean(array([aice_nodes_diff[expt,elm.nodes[0].id], aice_nodes_diff[expt,elm.nodes[1].id], aice_nodes_diff[expt,elm.nodes[2].id]]))
            i += 1

    # Truncate difference colourmap
    min_colour = 0
    max_colour = (amax(aice_diff)+1)/2.0
    diff_cmap = truncate_colormap(get_cmap('RdBu_r'), min_colour, max_colour)

    print 'Plotting'
    fig = figure(figsize=(10,8))
    gs = GridSpec(2,3)
    gs.update(left=0.1, right=0.9, bottom=0.1, top=0.85, wspace=0.01, hspace=0.15)
    # 1996-2005
    ax = subplot(gs[0,0], aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(aice_beg)
    img.set_clim(vmin=0, vmax=1)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    title('1996-2005', fontsize=20)
    # Colourbar on the left
    cbaxes = fig.add_axes([0.02, 0.58, 0.02, 0.2])
    cbar = colorbar(img, cax=cbaxes, ticks=arange(0, 1+0.25, 0.25))
    cbar.ax.tick_params(labelsize=14)
    # Loop over the rest of the experiments
    for expt in range(num_expts):
        ax = subplot(gs[j_plot[expt],i_plot[expt]], aspect='equal')
        img = PatchCollection(patches, cmap=diff_cmap)
        img.set_array(aice_diff[expt,:])
        img.set_clim(vmin=-1, vmax=amax(aice_diff))
        img.set_edgecolor('face')
        ax.add_collection(img)
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        ax.set_xticks([])
        ax.set_yticks([])
        title(expt_names[expt], fontsize=20)
        if expt == 1:
            # Colourbar on the right
            cbaxes = fig.add_axes([0.92, 0.58, 0.02, 0.2])
            cbar = colorbar(img, cax=cbaxes, ticks=arange(-1, 1+0.5, 0.5))
            cbar.ax.tick_params(labelsize=14)
        if expt == num_expts-1:
            # Text to indicate anomalies
            text(x_min, y_min-3, 'anomalies (2091-2100 minus 1996-2005)', ha='left', va='top', fontsize=18)
    # Main title
    suptitle('September sea ice concentration', fontsize=28)
    fig.show()
    fig.savefig('rcp_seaice.png')


# Truncate colourmap function from https://stackoverflow.com/questions/40929467/how-to-use-and-plot-only-a-part-of-a-colorbar-in-matplotlib
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=-1):
    if n== -1:
        n = cmap.N
    new_cmap = LinearSegmentedColormap.from_list('trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval), cmap(linspace(minval, maxval, n)))
    return new_cmap


# Command-line interface
if __name__ == "__main__":

    rcp_seaice()
            
    
    
