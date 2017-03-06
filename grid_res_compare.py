from patches import *
from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.patches import Polygon

def grid_res_compare ():

    directory_head = '/short/y99/kaa561/FESOM/mesh/'
    mesh_low = directory_head + 'low_res/'
    mesh_high = directory_head + 'high_res/'
    circumpolar = True
    lat_max = -60 + 90
    res_max = 25

    elements_low, patches_low = make_patches(mesh_low, circumpolar)
    elements_high, patches_high = make_patches(mesh_high, circumpolar)

    values_low = []
    for elm in elements_low:
        values_low.append(sqrt(elm.area())*1e-3)
    values_high = []
    for elm in elements_high:
        values_high.append(sqrt(elm.area())*1e-3)

    fig = figure(figsize=(20,9))
    ax = fig.add_subplot(1,2,1, aspect='equal')
    img = PatchCollection(patches_low, cmap=jet)
    img.set_array(array(values_low))
    img.set_clim(vmin=0, vmax=res_max)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    axis('off')
    title('Low res', fontsize=24)
    ax = fig.add_subplot(1,2,2, aspect='equal')
    img = PatchCollection(patches_high, cmap=jet)
    img.set_array(array(values_high))
    img.set_clim(vmin=0, vmax=res_max)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    axis('off')
    title('High res', fontsize=24)
    cbaxes = fig.add_axes([0.92, 0.2, 0.01, 0.6])
    cbar = colorbar(img, cax=cbaxes, extend='max', ticks=arange(0,25+5,5))
    cbar.ax.tick_params(labelsize=20)
    suptitle('Horizontal grid resolution (km)', fontsize=30)

    #fig.show()
    fig.savefig('grid_res_compare.png')


# Command-line interface
if __name__ == "__main__":

    grid_res_compare()
