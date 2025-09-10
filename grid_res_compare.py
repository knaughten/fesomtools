from patches import *
from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.patches import Polygon

# Make a 2x1 plot showing horizontal resolution (square root of the area of each
# triangular Element) for the low-res and high-res FESOM mesh.
# Input:
# flag = 1 (global lat-lon, colourbar limit 225 km) or 2 (circumpolar Antarctic,
#        colourbar limit 20 km)
# save = optional boolean indicating to save the figure, rather than display
# fig_name = if save=True, filename for figure
def grid_res_compare (flag, save=False, fig_name=None):

    # Paths to mesh directories
    directory_head = '/short/y99/kaa561/FESOM/mesh/'
    mesh_low = directory_head + 'low_res/'
    mesh_high = directory_head + 'high_res/'

    # Set plotting parameters based on flag
    if flag == 1:
        circumpolar = False
        res_max = 225        
    elif flag == 2:
        circumpolar = True
        lat_max = -63 + 90
        res_max = 20

    # Make plotting patches for each mesh
    elements_low, patches_low = make_patches(mesh_low, circumpolar)
    elements_high, patches_high = make_patches(mesh_high, circumpolar)

    # Calculate grid resolution for each element: square root of area,
    # convert to km
    values_low = []
    for elm in elements_low:
        values_low.append(sqrt(elm.area())*1e-3)
    values_high = []
    for elm in elements_high:
        values_high.append(sqrt(elm.area())*1e-3)

    # Plot
    fig = figure(figsize=(20,9))
    # Low-res
    if circumpolar:        
        ax = fig.add_subplot(1,2,1, aspect='equal')
    else:
        ax = fig.add_subplot(1,2,1, aspect=1.5)
    img = PatchCollection(patches_low, cmap=jet)
    img.set_array(array(values_low))
    img.set_clim(vmin=0, vmax=res_max)
    img.set_edgecolor('face')
    ax.add_collection(img)
    if circumpolar:
        xlim([-lat_max, lat_max])
        ylim([-lat_max, lat_max])
    else:
        xlim([-180, 180])
        ylim([-90, 90])
    axis('off')
    title('a) Low resolution', fontsize=24)
    # High-res
    if circumpolar:
        ax = fig.add_subplot(1,2,2, aspect='equal')
    else:
        ax = fig.add_subplot(1,2,2, aspect=1.5)
    img = PatchCollection(patches_high, cmap=jet)
    img.set_array(array(values_high))
    img.set_clim(vmin=0, vmax=res_max)
    img.set_edgecolor('face')
    ax.add_collection(img)
    if circumpolar:
        xlim([-lat_max, lat_max])
        ylim([-lat_max, lat_max])
    else:
        xlim([-180, 180])
        ylim([-90, 90])
    axis('off')
    title('b) High resolution', fontsize=24)
    cbaxes = fig.add_axes([0.92, 0.2, 0.01, 0.6])
    cbar = colorbar(img, cax=cbaxes, extend='max')
    cbar.ax.tick_params(labelsize=20)
    suptitle('Horizontal grid resolution (km)', fontsize=30)
    subplots_adjust(wspace=0.05)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    flag = int(input("Global (1) or circumpolar Antarctic (2)? "))
    action = input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = input("Filename for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    grid_res_compare(flag, save, fig_name)
