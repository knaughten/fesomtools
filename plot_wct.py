from patches import *
from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.patches import Polygon

# Calculate and plot the water column thickness for the given FESOM grid.
# Input:
# mesh_path = path to FESOM mesh directory
# fig_name = filename for figure
def plot_wct (mesh_path, fig_name):

    # Plotting parameters
    circumpolar = True
    lat_max = -30 + 90
    font_sizes = [240, 192, 160]

    # Build triangular patches for each element
    elements, patches = make_patches(mesh_path, circumpolar)

    # Read the ice shelf draft for each node
    file = open(mesh_path + 'shelf.out', 'r')
    node_shelf = []
    for line in file:
        node_shelf.append(-float(line))
    file.close()

    # Read the depth for each node
    file = open(mesh_path + 'depth.out', 'r')
    node_depth = []
    for line in file:
        node_depth.append(-float(line))
    file.close()

    # Calculate water column thickness
    node_wct = abs(array(node_depth)) - abs(array(node_shelf))

    # For each element, calculate the minimum wct of the 3 component nodes
    elm_wct = []
    for elm in elements:
        wct1 = node_wct[(elm.nodes[0]).id-1]
        wct2 = node_wct[(elm.nodes[1]).id-1]
        wct3 = node_wct[(elm.nodes[2]).id-1]
        elm_wct.append(amin(array([wct1, wct2, wct3])))

    # Set up figure
    fig = figure(figsize=(128, 96))
    ax = fig.add_subplot(1,1,1, aspect='equal')
    # Set colours for patches and add them to plot
    img = PatchCollection(patches, cmap=jet)
    img.set_array(array(elm_wct))
    img.set_edgecolor('face')
    ax.add_collection(img)

    # Configure plot
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    title('Water column thickness (m)', fontsize=font_sizes[0])
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=font_sizes[2])
    axis('off')

    # Plot specified points
    #problem_ids = [291, 297]
    #problem_x = []
    #problem_y = []
    #for elm in elements:
        #for i in range(3):
            #if elm.nodes[i].id in problem_ids:
                #problem_x.append(elm.x[i])
                #problem_y.append(elm.y[i])
                #problem_ids.remove(elm.nodes[i].id)
    #ax.plot(problem_x, problem_y, 'or')    

    savefig(fig_name)


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to mesh directory: ")
    fig_name = raw_input("Filename for figure: ")
    plot_wct(mesh_path, fig_name)
