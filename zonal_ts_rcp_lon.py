from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from fesom_grid import *
from fesom_sidegrid import *

# Create a 2x2 plot showing zonal slices of temperature (top) and salinity
# (bottom) through the given longitude, comparing the first 10 years of the
# given RCP (left) to the last 10 years (right).
# Input:
# lon0 = longitude to slice through (-180 to 180)
# lat_min, lat_max = bounds on latitude to plot (-90 to 90)
# mesh_path = path to FESOM mesh directory
# output_path = path to output directory for the given RCP, containing the
#               files annual_avg.oce.mean.2006.2015.nc and
#               annual_avg.oce.mean.nc.2091.2100.nc, which have 3D temperature
#               and salinity averaged over the first 10 years and last 10 years
#               respectively
# save = optional boolean indicating to save the figure to a file, rather than
#        display it on the screen
# fig_name = if save=True, filename for figure
def zonal_ts_rcp_lon (lon0, lat_min, lat_max, mesh_path, output_path, save=False, fig_name=None):

    file_name_beg = 'annual_avg.oce.mean.2006.2015.nc'
    file_name_end = 'annual_avg.oce.mean.2091.2100.nc'

    print 'Building FESOM mesh'
    elm2D = fesom_grid(mesh_path)
    print 'Reading temperature and salinity data'
    id = Dataset(output_path + file_name_beg, 'r')
    temp_nodes_beg = id.variables['temp'][0,:]
    salt_nodes_beg = id.variables['salt'][0,:]
    id.close()
    id = Dataset(output_path + file_name_end, 'r')
    temp_nodes_end = id.variables['temp'][0,:]
    salt_nodes_end = id.variables['salt'][0,:]
    id.close()

    # Figure out what to write on the title about longitude
    if lon0 < 0:
        lon_string = str(-lon0)+r'$^{\circ}$W'
    else:
        lon_string = str(lon0)+r'$^{\circ}$E'

    # Build arrays of SideElements making up zonal slices
    selements_temp_beg = fesom_sidegrid(elm2D, temp_nodes_beg, lon0, lat_max)
    selements_salt_beg = fesom_sidegrid(elm2D, salt_nodes_beg, lon0, lat_max)
    selements_temp_end = fesom_sidegrid(elm2D, temp_nodes_end, lon0, lat_max)
    selements_salt_end = fesom_sidegrid(elm2D, salt_nodes_end, lon0, lat_max)

    # Build array of quadrilateral patches for the plots, and data values
    # corresponding to each SideElement
    patches = []
    temp_beg = []
    for selm in selements_temp_beg:
        # Make patch
        coord = transpose(vstack((selm.y, selm.z)))
        patches.append(Polygon(coord, True, linewidth=0.))
        # Save data value
        temp_beg.append(selm.var)
    temp_beg = array(temp_beg)
    # Other variables have same patches but different values
    salt_beg = []
    for selm in selements_salt_beg:
        salt_beg.append(selm.var)
    salt_beg = array(salt_beg)
    temp_end = []
    for selm in selements_temp_end:
        temp_end.append(selm.var)
    temp_end = array(temp_end)
    salt_end = []
    for selm in selements_salt_end:
        salt_end.append(selm.var)
    salt_end = array(salt_end)
    # Find bounds on each variable
    temp_min = min(amin(temp_beg), amin(temp_end))
    temp_max = max(amax(temp_beg), amax(temp_end))
    salt_min = min(amin(salt_beg), amin(salt_beg))
    salt_max = max(amax(salt_end), amax(salt_end))
    # Find deepest depth
    # Start with 0
    depth_min = 0
    # Modify with patches
    for selm in selements_temp_beg:
        depth_min = min(depth_min, amin(selm.z))
    # Round down to nearest 50 metres
    depth_min = floor(depth_min/50)*50
    # Plot
    fig = figure(figsize=(18,12))
    # Temperature, first 10 years
    ax = fig.add_subplot(2, 2, 1)
    img1 = PatchCollection(patches, cmap='jet')
    img1.set_array(temp_beg)
    img1.set_edgecolor('face')
    img1.set_clim(vmin=temp_min, vmax=temp_max)
    ax.add_collection(img1)
    xlim([lat_min, lat_max])
    ylim([depth_min, 0])
    title(r'Temperature ($^{\circ}$C), 2006-2015', fontsize=20)
    ylabel('Depth (m)', fontsize=16)
    # Temperature, last 10 years
    ax = fig.add_subplot(2, 2, 2)
    img2 = PatchCollection(patches, cmap='jet')
    img2.set_array(temp_end)
    img2.set_edgecolor('face')
    img2.set_clim(vmin=temp_min, vmax=temp_max)
    ax.add_collection(img2)
    xlim([lat_min, lat_max])
    ylim([depth_min, 0])
    title(r'Temperature ($^{\circ}$C), 2091-2100', fontsize=20)
    # Add colorbar for temperature
    cbaxes_temp = fig.add_axes([0.92, 0.575, 0.01, 0.3])
    cbar_temp = colorbar(img2, cax=cbaxes_temp)
    cbar_temp.ax.tick_params(labelsize=16)
    # Salinity, first 10 years
    ax = fig.add_subplot(2, 2, 3)
    img3 = PatchCollection(patches, cmap='jet')
    img3.set_array(salt_beg)
    img3.set_edgecolor('face')
    img3.set_clim(vmin=salt_min, vmax=salt_max)
    ax.add_collection(img3)
    xlim([lat_min, lat_max])
    ylim([depth_min, 0])
    title('Salinity (psu), 2006-2015', fontsize=20)
    xlabel('Latitude', fontsize=16)
    ylabel('Depth (m)', fontsize=16)
    # Salinity, last 10 years
    ax = fig.add_subplot(2, 2, 4)
    img4 = PatchCollection(patches, cmap='jet')
    img4.set_array(salt_end)
    img4.set_edgecolor('face')
    img4.set_clim(vmin=salt_min, vmax=salt_max)
    ax.add_collection(img4)
    xlim([lat_min, lat_max])
    ylim([depth_min, 0])
    title('Salinity (psu), 2091-2100', fontsize=20)
    xlabel('Latitude', fontsize=16)
    # Add colorbar for salinity
    cbaxes_salt = fig.add_axes([0.92, 0.125, 0.01, 0.3])
    cbar_salt = colorbar(img4, cax=cbaxes_salt)
    cbar_salt.ax.tick_params(labelsize=16)
    # Main title
    suptitle(lon_string, fontsize=28)
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    lon0 = float(raw_input('Longitude to plot (-180 to 180): '))
    lat_min = float(raw_input('Minimum latitude to plot (-90 to 90): '))
    lat_max = float(raw_input('Maximum latitude to plot (-90 to 90): '))
    mesh_path = raw_input("Path to FESOM mesh directory: ")
    output_path = raw_input("Path to output directory for RCP: ")
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input('Filename for figure: ')
    elif action == 'd':
        save = False
        fig_name = None
    zonal_ts_rcp_lon(lon0, lat_min, lat_max, mesh_path, output_path, save, fig_name)
