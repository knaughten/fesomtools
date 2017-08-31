from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from patches import *

def thdgr_peninsula_res (mesh_path_lr, seasonal_file_lr, mesh_path_hr, seasonal_file_hr):

    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Bounds on plot (for polar coordinate transformation)
    x_min = -22.5
    x_max = -12
    y_min = 6
    y_max = 15
    # Plotting parameters
    circumpolar = True
    mask_cavities = True
    # Season names for plot titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    # Colour bounds
    bound = 4

    print 'Processing low-res FESOM'
    # Build mesh
    elements_lr, patches_lr = make_patches(mesh_path_lr, circumpolar, mask_cavities)
    # Read data
    id = Dataset(seasonal_file_lr, 'r')
    thdgr_nodes_lr = id.variables['thdgr'][:,:]*1e7
    id.close()
    # Count the number of elements not in ice shelf cavities
    num_elm_lr = 0
    for elm in elements_lr:
        if not elm.cavity:
            num_elm_lr += 1
    # Set up array for element-averages for each season
    thdgr_lr = zeros([4, num_elm_lr])
    # Loop over elements to fill this in
    i = 0
    for elm in elements_lr:
        if not elm.cavity:
            # Average over 3 component nodes
            thdgr_lr[:,i] = (thdgr_nodes_lr[:,elm.nodes[0].id] + thdgr_nodes_lr[:,elm.nodes[1].id] + thdgr_nodes_lr[:,elm.nodes[2].id])/3
            i += 1

    print 'Processing high-res FESOM'
    elements_hr, patches_hr = make_patches(mesh_path_hr, circumpolar, mask_cavities)
    id = Dataset(seasonal_file_hr, 'r')
    thdgr_nodes_hr = id.variables['thdgr'][:,:]*1e7
    id.close()
    num_elm_hr = 0
    for elm in elements_hr:
        if not elm.cavity:
            num_elm_hr += 1
    thdgr_hr = zeros([4, num_elm_hr])
    i = 0
    for elm in elements_hr:
        if not elm.cavity:
            thdgr_hr[:,i] = (thdgr_nodes_hr[:,elm.nodes[0].id] + thdgr_nodes_hr[:,elm.nodes[1].id] + thdgr_nodes_hr[:,elm.nodes[2].id])/3
            i += 1

    print 'Plotting'
    fig = figure(figsize=(19,9))
    for season in range(4):
        # Low-res
        ax = fig.add_subplot(2, 4, season+1, aspect='equal')
        img = PatchCollection(patches_lr, cmap='RdBu_r')
        img.set_array(thdgr_lr[season,:])
        img.set_clim(vmin=-bound, vmax=bound)
        img.set_edgecolor('face')
        ax.add_collection(img)
        title(season_names[season], fontsize=24)
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        axis('off')
        if season == 0:
            text(-24, 10, 'low-res', fontsize=24, ha='right')
        # High-res
        ax = fig.add_subplot(2, 4, season+5, aspect='equal')
        img = PatchCollection(patches_hr, cmap='RdBu_r')
        img.set_array(thdgr_hr[season,:])
        img.set_clim(vmin=-bound, vmax=bound)
        img.set_edgecolor('face')
        ax.add_collection(img)
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        axis('off')
        if season == 0:
            text(-24, 10, 'high-res', fontsize=24, ha='right')
    cbaxes = fig.add_axes([0.35, 0.04, 0.3, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='max', ticks=arange(-bound, bound+2, 2))
    cbar.ax.tick_params(labelsize=20)
    suptitle(r'FESOM sea ice thermodynamic growth rate (10$^{-7}$ m/s), 1992-2016 average', fontsize=30)
    subplots_adjust(wspace=0.025,hspace=0.025)
    fig.show()
    fig.savefig('thdgr_peninsula_res.png')


# Command-line interface
if __name__ == "__main__":

    mesh_path_lr = raw_input("Path to low-res mesh directory: ")
    seasonal_file_lr = raw_input("Path to low-res seasonal climatology file containing thdgr: ")
    mesh_path_hr = raw_input("Path to high-res mesh directory: ")
    seasonal_file_hr = raw_input("Path to high-res seasonal climatology file containing thdgr: ")
    thdgr_peninsula_res(mesh_path_lr, seasonal_file_lr, mesh_path_hr, seasonal_file_hr)
            

