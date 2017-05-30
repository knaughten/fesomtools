from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.patches import Polygon
from patches import *

# Make a 2x1 plot showing bottom water temperature and salinity anomalies
# for the RCP 8.5 M HR experiment (2091-2100 average minus 2006-2015 average).
def bw_rcp85_m_hr ():

    # Path to FESOM mesh directory
    mesh_path = '../FESOM/mesh/high_res/'
    # Path to FESOM oce.mean file containing 2091-2100 average minus 2006-2015
    # average (single annual timestep, not climatology)
    file_path = '../FESOM/rcp85_M_highres/output/annual_avg.oce.mean.diff.nc'
    # Boundary on plot
    lat_max = -63 + 90
    # Plotting parameters
    circumpolar = True
    mask_cavities = False
    # Limits on colourbars
    temp_lim = 2.0
    salt_lim = 0.5

    # Make plotting patches
    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)

    # Read temperature and salinity anomalies at each node
    id = Dataset(file_path, 'r')
    temp = id.variables['temp'][0,:]
    salt = id.variables['salt'][0,:]
    id.close()

    # Calculate values at each element
    temp_vals = []
    salt_vals = []
    for elm in elements:
        temp_vals_tmp = []
        salt_vals_tmp = []
        # Find the bottommost nodes
        for node in elm.nodes:
            id = node.find_bottom().id
            temp_vals_tmp.append(temp[id])
            salt_vals_tmp.append(salt[id])
        # Average over the 3 corners
        temp_vals.append(mean(temp_vals_tmp))
        salt_vals.append(mean(salt_vals_tmp))

    # Plot
    fig = figure(figsize=(20,10))
    # Temperature
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
    # Colourbar on the left
    cbaxes1 = fig.add_axes([0.05, 0.15, 0.02, 0.7])
    cbar1 = colorbar(img1, cax=cbaxes1, extend='both')
    cbar1.ax.tick_params(labelsize=14)
    # Salinity
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
    # Colourbar on the right
    cbaxes2 = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar2 = colorbar(img2, cax=cbaxes2, extend='both')
    cbar2.ax.tick_params(labelsize=16)
    # Main title
    suptitle('RCP 8.5 M HR: 2091-2100 minus 2006-2015', fontsize=30)

    fig.savefig('bw_rcp85m_hr.png')


# Command-line interface
if __name__ == "__main__":

    bw_rcp85_m_hr()
