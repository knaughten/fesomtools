from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from patches import *
from seasonal_avg import *

# Creates a 4x2 plot of seasonally averaged sea surface temperature (top row)
# and salinity (bottom row) over the last year of simulation.
# Input:
# mesh_path = path to FESOM mesh directory
# file_path1 = path to a FESOM output file containing one year of 5-day
#              averages for ocean variables (we will just use December)
# file_path2 = path to a FESOM output file containing the following year of
#              5-day averages for ocean variables (we will use January
#              through November)
# save = optional boolean indicating to save the figure to a file, rather than
#        display it on the screen
# fig_name = if save=True, filename for figure
def sst_sss_seasonal (mesh_path, file_path1, file_path2, save=False, fig_name=None):

    # FESOM parameters
    circumpolar=True
    mask_cavities=True
    # Season names for plot titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']

    # Build FESOM mesh
    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)

    # Figure out how many 2D nodes there are
    file = open(mesh_path + 'nod2d.out', 'r')
    file.readline()
    n2d = 0
    for line in file:
        n2d += 1
    file.close()

    # Get seasonal averages of the 3D FESOM output
    temp = seasonal_avg(file_path1, file_path2, 'temp')
    salt = seasonal_avg(file_path1, file_path2, 'salt')
    # Select the surface layer
    sst = temp[:,:n2d]
    sss = salt[:,:n2d]

    # Plot
    fig = figure(figsize=(20,9))
    # Loop over seasons
    for season in range(4):
        # SST
        # Build an array of FESOM data values corresponding to each Element
        values1 = []
        for elm in elements:
            # For each element not in an ice shelf cavity, append the mean
            # value for the 3 component Nodes
            if not elm.cavity:
                values1.append(mean([sst[season,elm.nodes[0].id], sst[season,elm.nodes[1].id], sst[season,elm.nodes[2].id]]))
        ax = fig.add_subplot(2, 4, season+1, aspect='equal')
        img = PatchCollection(patches, cmap='RdBu_r' )#jet)
        img.set_array(array(values1))
        #img.set_clim(vmin=-2, vmax=10)
        img.set_clim(vmin=-4, vmax=4)
        img.set_edgecolor('face')
        ax.add_collection(img)
        xlim([-35, 35])
        ylim([-33, 37])
        axis('off')
        if season == 0:
            text(-39, 0, r'SST ($^{\circ}$C)', fontsize=21, ha='right')
        title(season_names[season], fontsize=24)
        if season == 3:
            cbaxes1 = fig.add_axes([0.92, 0.55, 0.01, 0.3])            
            #cbar1 = colorbar(img, ticks=arange(-2,10+4,4), cax=cbaxes1)
            cbar1 = colorbar(img, ticks=arange(-4,4+2,2), cax=cbaxes1)
            cbar1.ax.tick_params(labelsize=16)
        # SSS
        values2 = []
        for elm in elements:
            # For each element not in an ice shelf cavity, append the mean
            # value for the 3 component Nodes
            if not elm.cavity:
                values2.append(mean([sss[season,elm.nodes[0].id], sss[season,elm.nodes[1].id], sss[season,elm.nodes[2].id]]))
        ax = fig.add_subplot(2, 4, season+5, aspect='equal')
        img = PatchCollection(patches, cmap='RdBu_r') #jet)
        img.set_array(array(values2))
        #img.set_clim(vmin=33, vmax=35)
        img.set_clim(vmin=-0.5, vmax=0.5)
        img.set_edgecolor('face')
        ax.add_collection(img)
        xlim([-35, 35])
        ylim([-33, 37])
        axis('off')
        if season == 0:
            text(-39, 0, 'SSS (psu)', fontsize=21, ha='right')
        if season == 3:
            cbaxes2 = fig.add_axes([0.92, 0.15, 0.01, 0.3])
            #cbar2 = colorbar(img, ticks=arange(33,35+0.5,0.5), cax=cbaxes2)
            cbar2 = colorbar(img, ticks=arange(-0.5,0.5+0.25,0.25), cax=cbaxes2)
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

    mesh_path = input("Path to FESOM mesh directory: ")
    file_path1 = input("Path to output oce.mean.nc containing one year of 5-day averages (December will be used): ")
    file_path2 = input("Path to the following oce.mean.nc containing 5-day averages for the next year (January through November will be used): ")
    action = input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None

    sst_sss_seasonal(mesh_path, file_path1, file_path2, save, fig_name)

    

    

    
