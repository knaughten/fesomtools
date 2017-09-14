from patches import *
from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.patches import Polygon

def rcp_resolution ():

    # Path to mesh directory
    mesh_path = '/short/y99/kaa561/FESOM/mesh/high_res/'
    # Colourbar limit for circumpolar view
    circumpolar_max = 15
    # Northern boundary for circumpolar view: 63S
    nbdry = -63 + 90

    # Build 2 versions of the mesh: global and circumpolar
    elements_global, patches_global = make_patches(mesh_path, circumpolar=False)
    elements_circumpolar, patches_circumpolar = make_patches(mesh_path, circumpolar=True)

    # Calculate grid resolution for each element: square root of area, convert
    # to km
    res_global = []
    for elm in elements_global:
        res_global.append(sqrt(elm.area())*1e-3)
    res_circumpolar = []
    for elm in elements_circumpolar:
        res_circumpolar.append(sqrt(elm.area())*1e-3)

    # Plot
    fig = figure(figsize=(20,9))
    # Global
    gs_a = GridSpec(1,1)
    gs_a.update(left=0.1, right=0.54, bottom=0.1, top=0.85)
    ax = subplot(gs_a[0,0])
    img = PatchCollection(patches_global, cmap='jet')
    img.set_array(array(res_global))
    img.set_clim(vmin=0, vmax=amax(res_global))
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-180, 180])
    ylim([-90, 90])
    ax.set_xticks([])
    ax.set_yticks([])
    title('a) Global', fontsize=28)
    # Colourbar on the left
    cbaxes = fig.add_axes([0.04, 0.3, 0.02, 0.4])
    cbar = colorbar(img, cax=cbaxes, ticks=arange(0,200+50,50))
    cbar.ax.tick_params(labelsize=20)
    # Circumpolar
    gs_b = GridSpec(1,1)
    gs_b.update(left=0.55, right=0.9, bottom=0.1, top=0.85)
    ax = subplot(gs_b[0,0], aspect='equal')
    img = PatchCollection(patches_circumpolar, cmap='jet')
    img.set_array(array(res_circumpolar))
    img.set_clim(vmin=0, vmax=circumpolar_max)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    ax.set_xticks([])
    ax.set_yticks([])
    title('b) Antarctic', fontsize=28)
    # Colourbar on the right
    cbaxes = fig.add_axes([0.92, 0.3, 0.02, 0.4])
    cbar = colorbar(img, cax=cbaxes, ticks=arange(0,15+3,3), extend='max')
    cbar.ax.tick_params(labelsize=20)
    suptitle('Horizontal resolution (km)', fontsize=32)
    fig.show()
    fig.savefig('resolution.png')


# Command-line interface
if __name__ == "__main__":

    rcp_resolution()
    
    
