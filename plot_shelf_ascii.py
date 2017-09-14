from patches import *
from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.patches import Polygon

def plot_shelf_ascii (mesh_path, shelf_file, fig_name, circumpolar=True):

    # Plotting parameters
    if circumpolar:
        lat_max = -30 + 90
        font_sizes = [240, 192, 160]
    else:
        font_sizes = [30, 24, 20]

    # Build triangular patches for each element
    elements, patches = make_patches(mesh_path, circumpolar)

    # Read the ice shelf draft for each node
    node_shelf = []
    f = open(shelf_file, 'r')
    for line in f:
        node_shelf.append(-1*float(line))
    f.close()

    # Find the ice shelf draft for each element
    elm_shelf = []
    for elm in elements:
        #shelf1 = (elm.nodes[0].find_bottom()).shelf
        #shelf2 = (elm.nodes[1].find_bottom()).shelf
        #shelf3 = (elm.nodes[2].find_bottom()).shelf
        elm_shelf.append(mean(array([node_shelf[elm.nodes[0].id], node_shelf[elm.nodes[1].id], node_shelf[elm.nodes[2].id]])))

    # Set up figure
    if circumpolar:
        fig = figure(figsize=(128, 96))
        ax = fig.add_subplot(1,1,1, aspect='equal')
    else:
        fig = figure(figsize=(16, 8))
        ax = fig.add_subplot(1,1,1)    
    # Set colours for patches and add them to plot
    img = PatchCollection(patches, cmap=jet)
    img.set_array(array(elm_shelf))
    img.set_edgecolor('face')
    ax.add_collection(img)

    # Configure plot
    if circumpolar:
        xlim([-lat_max, lat_max])
        ylim([-lat_max, lat_max])
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        axis('off')
    else:
        xlim([-180, 180])
        ylim([-90, 90])
        ax.get_xaxis().set_ticks(arange(-120,120+1,60))
        ax.get_yaxis().set_ticks(arange(-60,60+1,30))
    title('Ice shelf draft (m)', fontsize=font_sizes[0])
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=font_sizes[2])
    img.set_clim(vmin=0, vmax=max(elm_shelf))

    savefig(fig_name)


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to mesh directory: ")
    shelf_file = raw_input("Path to shelf.out file: ")
    fig_name = raw_input("Filename for figure: ")
    domain = raw_input("Global (g) or circumpolar (c)? ")
    if domain == 'c':
        circumpolar = True
    elif domain == 'g':
        circumpolar = False
    plot_shelf_ascii(mesh_path, shelf_file, fig_name, circumpolar)
