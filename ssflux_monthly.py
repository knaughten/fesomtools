from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.patches import Polygon
from patches import *
from monthly_avg import *

def ssflux_monthly (elements, patches, file_path, month, bound=None, save=False, fig_name=None):

    # Month names for plot titles
    month_name = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
    ndays_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    lat_max = -30+90

    # Get monthly average of the FESOM output
    ssflux = monthly_avg(file_path, 'virtual_salt', month)
    ssflux *= 1e6

    # Build an array of data values corresponding to each Element
    values = []
    for elm in elements:
        values.append(mean([ssflux[elm.nodes[0].id], ssflux[elm.nodes[1].id], ssflux[elm.nodes[2].id]]))
    if bound is None:
        bound = amax(abs(array(values)))

    fig = figure(figsize=(16,12))
    ax = fig.add_subplot(1,1,1,aspect='equal')
    img = PatchCollection(patches, cmap='RdYlBu_r')
    img.set_array(array(values))
    img.set_edgecolor('face')
    img.set_clim(vmin=-bound, vmax=bound)
    ax.add_collection(img)
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    axis('off')
    title(month_name[month] + r' surface salinity flux (10$^{-6}$ kg/m$^2$/s)', fontsize=30)
    cbar = colorbar(img, extend='both')
    cbar.ax.tick_params(labelsize=20)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


if __name__ == "__main__":

    circumpolar=True
    mask_cavities=False

    mesh_path = raw_input("Path to mesh directory: ")
    file_path = raw_input("Path to FESOM forcing.diag.nc output file: ")
    month = int(raw_input("Month number (1-12): "))-1
    bound = None
    get_bound = raw_input("Set bounds on colour scale (y/n)? ")
    if get_bound == 'y':
        bound = float(raw_input("Maximum absolute value (1e-6 kg/m^2/s: "))
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None

    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)
    ssflux_monthly(elements, patches, file_path, month, bound, save, fig_name)

    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            update_mesh = False
            while True:
                changes = raw_input("Enter a parameter to change: (1) mesh path, (2) file path, (3) month, (4) colour bounds, (5) save/display; or enter to continue: ")
                if len(changes) == 0:
                    break
                else:
                    if int(changes) == 1:
                        update_mesh = True
                        mesh_path = raw_input("Path to mesh directory: ")
                        file_path = raw_input("Path to FESOM forcing.diag.nc output file: ")
                    elif int(changes) == 2:
                        file_path = raw_input("Path to FESOM forcing.diag.nc output file: ")
                    elif int(changes) == 3:
                        month = int(raw_input("Month number (1-12): "))-1
                    elif int(changes) == 4:
                        bound = None
                        get_bound = raw_input("Set bounds on colour scale (y/n)? ")
                        if get_bound == 'y':
                            bound = float(raw_input("Maximum absolute value (1e-6 kg/m^2/s: "))
                    elif int(changes) == 5:
                        save = not save
            if save:
                fig_name = raw_input("File name for figure: ")
            if update_mesh:
                elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)
            ssflux_monthly(elements, patches, file_path, month, bound, save, fig_name)
        else:
            break
        
    
