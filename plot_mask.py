from patches import *
from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.patches import Polygon

# Plot the mask created by restoring_mask.py and display on the screen.
# Input:
# mesh_path = path to FESOM mesh directory
# circumpolar = optional boolean indicating the plot should be circumpolar
#               Antarctic (True) as opposed to global lon-lat (False)
def plot_mask (mesh_path, circumpolar):

    # Plotting parameters
    if circumpolar:
        lat_max = -30 + 90

    # Build triangular patches for each element
    elements, patches = make_patches(mesh_path, circumpolar)

    # Read mask
    f = open(mesh_path + 'ssrestoring_mask.out', 'r')
    node_mask = []
    for line in f:
        node_mask.append(float(line))
    f.close()
    # Apply to elements
    mask = []
    for elm in elements:
        # For the purposes of plotting, unmask the element if at least 2 of the
        # 3 component nodes are unmasked
        if sum([node_mask[elm.nodes[0].id], node_mask[elm.nodes[1].id], node_mask[elm.nodes[2].id]]) >= 2:
            # Unmask this element
            mask.append(1)
        else:
            # Mask this element
            mask.append(0)

    # Set up figure
    fig = figure(figsize=(12, 9))
    if circumpolar:
        ax = fig.add_subplot(1,1,1, aspect='equal')
    else:
        ax = fig.add_subplot(1,1,1)
    # Set colours for patches and add them to plot
    img = PatchCollection(patches, cmap=jet)
    img.set_array(array(mask))
    img.set_edgecolor('face')
    ax.add_collection(img)

    # Configure plot
    if circumpolar:
        xlim([-lat_max, lat_max])
        ylim([-lat_max, lat_max])
    else:
        xlim([-180, 180])
        ylim([-90, 90])
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    title('Restoring mask')
    cbar = colorbar(img)
    img.set_clim(vmin=0, vmax=1)
    axis('off')

    show()


# Command-line interface
if __name__ == "__main__":

    mesh_path = input("Path to mesh directory: ")
    domain = input("Global (g) or circumpolar (c)? ")
    if domain == 'c':
        circumpolar = True
    elif domain == 'g':
        circumpolar = False
    plot_mask(mesh_path, circumpolar)
