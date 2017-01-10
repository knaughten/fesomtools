from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from patches import *

# Creates a 4x2 plot of seasonally averaged sea ice concentration (top row) and
# thickness (bottom row) over the last year of simulation.
# Input:
# mesh_path = path to FESOM mesh directory
# file_path1 = path to a FESOM output file containing one year of 5-day
#              averages for sea ice variables (we will just use December)
# file_path2 = path to a FESOM output file containing the following year of
#              5-day averages for sea ice variables (we will use January
#              through November)
# save = optional boolean indicating to save the figure to a file, rather than
#        display it on the screen
# fig_name = if save=True, filename for figure
def aice_hi_seasonal (mesh_path, file_path1, file_path2, save=False, fig_name=None):

    # FESOM parameters
    circumpolar=True
    mask_cavities=True
    # Season names for plot titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']

    # Build FESOM mesh
    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)

    # Get seasonal averages of the FESOM output
    # This is hard-coded and ugly
    id = Dataset(file_path1, 'r')
    n2d = id.variables['area'].shape[1]
    aice = zeros([4, n2d])
    hi = zeros([4, n2d])
    # DJF: 1/5 of index 67 (1-based) and indices 68-73 in file1; indices 1-11
    # and 4/5 of index 12 in file2; 90 days in total
    aice[0,:] = id.variables['area'][66,:] + sum(id.variables['area'][67:73,:]*5, axis=0)
    hi[0,:] = id.variables['hice'][66,:] + sum(id.variables['hice'][67:73,:]*5, axis=0)
    id.close()
    id = Dataset(file_path2, 'r')
    aice[0,:] += sum(id.variables['area'][0:11,:]*5, axis=0) + id.variables['area'][11,:]*4
    aice[0,:] /= 90
    hi[0,:] += sum(id.variables['hice'][0:11,:]*5, axis=0) + id.variables['hice'][11,:]*4
    hi[0,:] /= 90
    # MAM: 1/5 of index 12, indices 13-30, and 1/5 of index 31 in file2;
    # 92 days in total
    aice[1,:] = id.variables['area'][11,:] + sum(id.variables['area'][12:30,:]*5, axis=0) + id.variables['area'][30,:]
    aice[1,:] /= 92
    hi[1,:] = id.variables['hice'][11,:] + sum(id.variables['hice'][12:30,:]*5, axis=0) + id.variables['hice'][30,:]
    hi[1,:] /= 92
    # JJA: 4/5 of index 31, indices 32-48, and 3/5 of index 49 in file2;
    # 92 days in total
    aice[2,:] = id.variables['area'][30,:]*4 + sum(id.variables['area'][31:48,:]*5, axis=0) + id.variables['area'][48,:]*3
    aice[2,:] /= 92
    hi[2,:] = id.variables['hice'][30,:]*4 + sum(id.variables['hice'][31:48,:]*5, axis=0) + id.variables['hice'][48,:]*3
    hi[2,:] /= 92
    # SON: 2/5 of index 49, indices 50-66, and 4/5 of index 67 in file2;
    # 91 days in total
    aice[3,:] = id.variables['area'][48,:]*2 + sum(id.variables['area'][49:66,:]*5, axis=0) + id.variables['area'][66,:]*4
    aice[3,:] /= 91
    hi[3,:] = id.variables['hice'][48,:]*2 + sum(id.variables['hice'][49:66,:]*5, axis=0) + id.variables['hice'][66,:]*4
    hi[3,:] /= 91
    id.close()

    # Plot
    fig = figure(figsize=(20,9))
    # Loop over seasons
    for season in range(4):
        # aice
        # Build an array of FESOM data values corresponding to each Element
        values1 = []
        for elm in elements:
            # For each element not in an ice shelf cavity, append the mean
            # value for the 3 component Nodes
            if not elm.cavity:
                values1.append(mean([aice[season,elm.nodes[0].id], aice[season,elm.nodes[1].id], aice[season,elm.nodes[2].id]]))
        ax = fig.add_subplot(2, 4, season+1, aspect='equal')
        img = PatchCollection(patches, cmap=jet)
        img.set_array(array(values1))
        img.set_clim(vmin=0, vmax=1)
        img.set_edgecolor('face')
        ax.add_collection(img)
        xlim([-35, 35])
        ylim([-33, 37])
        axis('off')
        if season == 0:
            text(-39, 0, 'aice (%)', fontsize=21, ha='right')
        title(season_names[season], fontsize=24)
        if season == 3:
            cbaxes1 = fig.add_axes([0.92, 0.55, 0.01, 0.3])
            cbar1 = colorbar(img, ticks=arange(0,1+0.25,0.25), cax=cbaxes1)
            cbar1.ax.tick_params(labelsize=16)
        # hi
        values2 = []
        for elm in elements:
            # For each element not in an ice shelf cavity, append the mean
            # value for the 3 component Nodes
            if not elm.cavity:
                values2.append(mean([hi[season,elm.nodes[0].id], hi[season,elm.nodes[1].id], hi[season,elm.nodes[2].id]]))
        ax = fig.add_subplot(2, 4, season+5, aspect='equal')
        img = PatchCollection(patches, cmap=jet)
        img.set_array(array(values2))
        img.set_clim(vmin=0, vmax=2)
        img.set_edgecolor('face')
        ax.add_collection(img)
        xlim([-35, 35])
        ylim([-33, 37])
        axis('off')
        if season == 0:
            text(-39, 0, 'hi (m)', fontsize=21, ha='right')
        if season == 3:
            cbaxes2 = fig.add_axes([0.92, 0.15, 0.01, 0.3])
            cbar2 = colorbar(img, ticks=arange(0,2+0.5,0.5), cax=cbaxes2)
            cbar2.ax.tick_params(labelsize=16)
    # Decrease space between plots
    subplots_adjust(wspace=0.025,hspace=0.025)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()
        

# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    file_path1 = raw_input("Path to output ice.mean.nc containing one year of 5-day averages (December will be used): ")
    file_path2 = raw_input("Path to the following ice.mean.nc containing 5-day averages for the next year (January through November will be used): ")
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None

    aice_hi_seasonal(mesh_path, file_path1, file_path2, save, fig_name)

    

    

    
