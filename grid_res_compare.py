from patches import *
from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.patches import Polygon

def grid_res_compare (flag, save=False, fig_name=None):

    directory_head = '/short/y99/kaa561/FESOM/mesh/'
    mesh_low = directory_head + 'low_res/'
    mesh_high = directory_head + 'high_res/'

    if flag == 1:
        circumpolar = False
        res_max = 225        
    elif flag == 2:
        circumpolar = True
        lat_max = -63 + 90
        res_max = 20

    elements_low, patches_low = make_patches(mesh_low, circumpolar)
    elements_high, patches_high = make_patches(mesh_high, circumpolar)

    values_low = []
    for elm in elements_low:
        values_low.append(sqrt(elm.area())*1e-3)
    values_high = []
    for elm in elements_high:
        values_high.append(sqrt(elm.area())*1e-3)

    fig = figure(figsize=(20,9))
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

    flag = int(raw_input("Global (1) or circumpolar Antarctic (2)? "))
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("Filename for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    grid_res_compare(flag, save, fig_name)
