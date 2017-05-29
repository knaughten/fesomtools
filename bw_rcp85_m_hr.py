from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.patches import Polygon
from patches import *

def bw_rcp85_m_hr ():

    mesh_path = '../FESOM/mesh/high_res/'
    file_path = '../FESOM/rcp85_M_highres/output/annual_avg.oce.mean.diff.nc'
    lat_max = -63 + 90
    circumpolar = True
    mask_cavities = False
    temp_lim = 2.0
    salt_lim = 0.5

    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)

    id = Dataset(file_path, 'r')
    temp = id.variables['temp'][0,:]
    salt = id.variables['salt'][0,:]
    id.close()

    temp_vals = []
    salt_vals = []
    for elm in elements:
        temp_vals_tmp = []
        salt_vals_tmp = []
        for node in elm.nodes:
            id = node.find_bottom().id
            temp_vals_tmp.append(temp[id])
            salt_vals_tmp.append(salt[id])
        temp_vals.append(mean(temp_vals_tmp))
        salt_vals.append(mean(salt_vals_tmp))

    fig = figure(figsize=(20,10))
    ax1 = fig.add_subplot(1,2,1,aspect='equal')
    img1 = PatchCollection(patches, cmap='RdBu_r')
    img1.set_array(array(temp_vals))
    img1.set_edgecolor('face')
    img1.set_clim(vmin=-temp_lim, vmax=temp_lim)
    ax1.add_collection(img1)
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    axis('off')
    title(r'Bottom water temperature ($^{\circ}$C)', fontsize=24)
    cbaxes1 = fig.add_axes([0.05, 0.15, 0.02, 0.7])
    cbar1 = colorbar(img1, cax=cbaxes1, extend='both')
    cbar1.ax.tick_params(labelsize=14)
    ax2 = fig.add_subplot(1,2,2,aspect='equal')
    img2 = PatchCollection(patches, cmap='RdBu_r')
    img2.set_array(array(salt_vals))
    img2.set_edgecolor('face')
    img2.set_clim(vmin=-salt_lim, vmax=salt_lim)
    ax2.add_collection(img2)
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    axis('off')
    title('Bottom water salinity (psu)', fontsize=24)
    cbaxes2 = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar2 = colorbar(img2, cax=cbaxes2, extend='both')
    cbar2.ax.tick_params(labelsize=16)
    suptitle('RCP 8.5 M HR: 2091-2100 minus 2006-2015', fontsize=30)

    fig.savefig('bw_rcp85m_hr.png')


if __name__ == "__main__":

    bw_rcp85_m_hr()
