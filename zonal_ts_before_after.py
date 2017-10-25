from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from fesom_grid import *
from fesom_sidegrid import *

def zonal_ts_before_after (lon0, lat_min, lat_max, rcp, model, save=False, fig_name=None):

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/meshB/'
    file_beg = '/short/y99/kaa561/FESOM/highres_spinup/annual_avg.oce.mean.1996.2005.nc'
    file_end = '/short/y99/kaa561/FESOM/rcp'+rcp+'_'+model+'/annual_avg.oce.mean.2091.2100.nc'
    # Figure out what to write on the title about longitude
    if lon0 < 0:
        lon_string = str(int(round(-lon0)))+r'$^{\circ}$W'
    else:
        lon_string = str(int(round(lon0)))+r'$^{\circ}$E'

    print 'Building FESOM mesh'
    elm2D = fesom_grid(mesh_path)
    print 'Reading temperature and salinity data'
    id = Dataset(file_beg, 'r')
    temp_nodes_beg = id.variables['temp'][0,:]
    salt_nodes_beg = id.variables['salt'][0,:]
    id.close()
    id = Dataset(file_end, 'r')
    temp_nodes_end = id.variables['temp'][0,:]
    salt_nodes_end = id.variables['salt'][0,:]

    print 'Interpolating to ' + str(lon0)
    # Build arrays of SideElements making up zonal slices
    # Start with beginning
    selements_temp_beg = fesom_sidegrid(elm2D, temp_nodes_beg, lon0, lat_max)
    selements_salt_beg = fesom_sidegrid(elm2D, salt_nodes_beg, lon0, lat_max)
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
    # Salinity has same patches but different values
    salt_beg = []
    for selm in selements_salt_beg:
        salt_beg.append(selm.var)
    salt_beg = array(salt_beg)
    # Repeat for end
    selements_temp_end = fesom_sidegrid(elm2D, temp_nodes_end, lon0, lat_max)
    selements_salt_end = fesom_sidegrid(elm2D, salt_nodes_end, lon0, lat_max)
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
    salt_min = min(amin(salt_beg), amin(salt_end))
    salt_max = max(amax(salt_beg), amax(salt_end))
    # Find deepest depth
    # Start with 0
    depth_min = 0
    # Modify with patches
    for selm in selements_temp_beg:
        depth_min = min(depth_min, amin(selm.z))
    # Round down to nearest 50 metres
    depth_min = floor(depth_min/50)*50

    print 'Plotting'
    fig = figure(figsize=(16,10))    
    # Temperature
    gs_temp = GridSpec(1,2)
    gs_temp.update(left=0.11, right=0.9, bottom=0.5, top=0.9, wspace=0.05, hspace=0.5)
    # Beginning
    ax = subplot(gs_temp[0,0])
    img = PatchCollection(patches, cmap='jet')
    img.set_array(temp_beg)
    img.set_edgecolor('face')
    img.set_clim(vmin=temp_min, vmax=temp_max)
    ax.add_collection(img)
    xlim([lat_min, lat_max])
    ylim([depth_min, 0])
    title(r'Temperature ($^{\circ}$C), 1996-2005', fontsize=24)
    ax.set_xticklabels([])
    ylabel('Depth (m)', fontsize=18)
    # End
    ax = subplot(gs_temp[0,1])
    img = PatchCollection(patches, cmap='jet')
    img.set_array(temp_end)
    img.set_edgecolor('face')
    img.set_clim(vmin=temp_min, vmax=temp_max)
    ax.add_collection(img)
    xlim([lat_min, lat_max])
    ylim([depth_min, 0])
    title(r'Temperature ($^{\circ}$C), 2091-2100', fontsize=24)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    # Add a colorbar on the right
    cbaxes = fig.add_axes([0.92, 0.55, 0.02, 0.3])
    cbar = colorbar(img, cax=cbaxes, extend='both')
    cbar.ax.tick_params(labelsize=16)
    # Salinity
    gs_salt = GridSpec(1,2)
    gs_salt.update(left=0.11, right=0.9, bottom=0.05, top=0.45, wspace=0.05, hspace=0.5)
    # Beginning
    ax = subplot(gs_salt[0,0])
    img = PatchCollection(patches, cmap='jet')
    img.set_array(salt_beg)
    img.set_edgecolor('face')
    img.set_clim(vmin=salt_min, vmax=salt_max)
    ax.add_collection(img)
    xlim([lat_min, lat_max])
    ylim([depth_min, 0])
    title('Salinity (psu), 1996-2005', fontsize=24)
    xlabel('Latitude', fontsize=18)
    ylabel('Depth (m)', fontsize=18)
    # End
    ax = subplot(gs_salt[0,1])
    img = PatchCollection(patches, cmap='jet')
    img.set_array(salt_end)
    img.set_edgecolor('face')
    img.set_clim(vmin=salt_min, vmax=salt_max)
    ax.add_collection(img)
    xlim([lat_min, lat_max])
    ylim([depth_min, 0])
    title('Salinity (psu), 2091-2100', fontsize=24)
    xlabel('Latitude', fontsize=18)
    ax.set_yticklabels([])
    # Add a colorbar on the right
    cbaxes = fig.add_axes([0.92, 0.1, 0.02, 0.3])
    cbar = colorbar(img, cax=cbaxes, extend='both')
    cbar.ax.tick_params(labelsize=16)
    # Main title
    suptitle('RCP ' + rcp[0] + '.' + rcp[1] + ' ' + model + ', ' + lon_string, fontsize=28)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    lon0 = float(raw_input('Longitude to plot (-180 to 180): '))
    lat_min = float(raw_input('Minimum latitude to plot (-90 to 90): '))
    lat_max = float(raw_input('Maximum latitude to plot (-90 to 90): '))
    key = int(raw_input('RCP 4.5 (4) or 8.5 (8)? '))
    if key == 4:
        rcp = '45'
    elif key == 8:
        rcp = '85'
    key = int(raw_input('Multi-model mean (1) or ACCESS 1.0 (2)? '))
    if key == 1:
        model = 'M'
    elif key == 2:
        model = 'A'
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input('Filename for figure: ')
    elif action == 'd':
        save = False
        fig_name = None
    zonal_ts_before_after (lon0, lat_min, lat_max, rcp, model, save=False, fig_name=None)
