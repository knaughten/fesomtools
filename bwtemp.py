from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from patches import *

# Make a circumpolar plot of bottom water temperature annually averaged over
# the given year of FESOM simulation.
# Input:
# mesh_path = path to FESOM mesh directory
# file_path = path to output file containing one year of ocean averages
# save = optional boolean indicating to save the figure, rather than display
#        it on the screen
# fig_name = if save=True, filename for figure
def bwtemp (mesh_path, file_path, save=False, fig_name=None):

    # Plotting parameters
    lat_max = -30 + 90
    circumpolar=True

    # Build FESOM mesh
    elements, patches = make_patches(mesh_path, circumpolar)

    # Calculate annual average of bottom water temperature
    file = Dataset(file_path, 'r')
    data = mean(file.variables['temp'][:,:], axis=0)
    file.close()
    values = []
    # Loop over elements
    for elm in elements:
        values_tmp = []
        # For each component node, find the bottom index; average over
        # all 3 such indices to get the value for this element
        for node in elm.nodes:
            id = node.find_bottom().id
            values_tmp.append(data[id])
        values.append(mean(values_tmp))

    # Plot
    fig = figure(figsize=(16,12))
    ax = fig.add_subplot(1,1,1,aspect='equal')
    img = PatchCollection(patches, cmap=jet)
    img.set_array(array(values))
    img.set_edgecolor('face')
    img.set_clim(vmin=-2.5, vmax=2.5)
    ax.add_collection(img)
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    axis('off')
    title(r'Bottom water temperature ($^{\circ}$C), annual average', fontsize=30)
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=20)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()
    

# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    file_path = raw_input("Path to FESOM oce.mean.nc containing one year of output: ")
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save=True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    bwtemp(mesh_path, file_path, save, fig_name)

    

    


    
