from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from fesom_grid import *
from fesom_sidegrid import *

# Input:
# lon0 = longitude to slice through (-180 to 180)
# lat_min, lat_max = bounds on latitude to plot (-90 to 90)
# save = optional boolean indicating to save the figure to a file, rather than
#        display it on the screen
# fig_name = if save=True, filename for figure
def zonal_ts_rcp_lon (lon0, lat_min, lat_max, save=False, fig_name=None):

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/high_res/'
    directory_beg = '/short/y99/kaa561/FESOM/highres_spinup/'
    directories = ['/short/y99/kaa561/FESOM/rcp45_M_highres/output/', '/short/y99/kaa561/FESOM/rcp45_A_highres/output/', '/short/y99/kaa561/FESOM/rcp85_M_highres/output/', '/short/y99/kaa561/FESOM/rcp85_A_highres/output/', '/short/y99/kaa561/FESOM/highres_spinup/']
    file_beg = 'annual_avg.oce.mean.1996.2005.nc'
    file_end = 'annual_avg.oce.mean.2091.2100.nc'
    # Titles for plotting
    expt_names = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A', 'CONTROL']
    num_expts = len(directories)
    # Start and end years for each period
    beg_years = [1996, 2005]
    end_years = [2091, 2100]
    # Figure out what to write on the title about longitude
    if lon0 < 0:
        lon_string = str(-lon0)+r'$^{\circ}$W'
    else:
        lon_string = str(lon0)+r'$^{\circ}$E'

    print 'Building FESOM mesh'
    elm2D = fesom_grid(mesh_path)
    print 'Reading temperature and salinity data'
    id = Dataset(directory_beg + file_beg, 'r')
    temp_nodes_beg = id.variables['temp'][0,:]
    salt_nodes_beg = id.variables['salt'][0,:]
    id.close()
    n3d = size(temp_nodes_beg)
    # Get anomalies for each experiment
    temp_nodes_diff = zeros([num_expts, n3d])
    salt_nodes_diff = zeros([num_expts, n3d])
    for expt in range(num_expts):
        id = Dataset(directories[expt] + file_end, 'r')
        temp_nodes_diff[expt,:] = id.variables['temp'][0,:] - temp_nodes_beg
        salt_nodes_diff[expt,:] = id.variables['salt'][0,:] - salt_nodes_beg

    print 'Interpolating to ' + lon_string
    # Build arrays of SideElements making up zonal slices
    # Start with beginning
    print '...'+str(beg_years[0])+'-'+str(beg_years[1])
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
    # Repeat for anomalies for each experiment
    temp_diff = zeros([num_expts, len(patches)])
    salt_diff = zeros([num_expts, len(patches)])
    for expt in range(num_expts):
        print '...'+expt_names[expt]
        selements_temp_diff = fesom_sidegrid(elm2D, temp_nodes_diff[expt,:], lon0, lat_max)
        selements_salt_diff = fesom_sidegrid(elm2D, salt_nodes_diff[expt,:], lon0, lat_max)
        i = 0
        for selm in selements_temp_diff:
            temp_diff[expt,i] = selm.var
            i += 1
        i = 0
        for selm in selements_salt_diff:
            salt_diff[expt,i] = selm.var
            i += 1
    # Find bounds on each variable
    temp_min = amin(temp_beg)
    temp_max = amax(temp_beg)
    salt_min = amin(salt_beg)
    salt_max = amax(salt_beg)
    temp_max_diff = 2.0
    salt_max_diff = 0.5
    # Find deepest depth
    # Start with 0
    depth_min = 0
    # Modify with patches
    for selm in selements_temp_beg:
        depth_min = min(depth_min, amin(selm.z))
    # Round down to nearest 50 metres
    depth_min = floor(depth_min/50)*50

    print 'Plotting'
    fig = figure(figsize=(20,16))    
    # Temperature
    gs_temp = GridSpec(2,3)
    gs_temp.update(left=0.11, right=0.9, bottom=0.51, top=0.9, wspace=0.05, hspace=0.5)
    # Beginning
    ax = subplot(gs_temp[0,0])
    img = PatchCollection(patches, cmap='jet')
    img.set_array(temp_beg)
    img.set_edgecolor('face')
    img.set_clim(vmin=temp_min, vmax=temp_max)
    ax.add_collection(img)
    xlim([lat_min, lat_max])
    ylim([depth_min, 0])
    title(str(beg_years[0])+'-'+str(beg_years[1]), fontsize=20)
    xlabel('Latitude', fontsize=16)
    ylabel('Depth (m)', fontsize=16)
    # Add a colorbar on the left
    cbaxes = fig.add_axes([0.02, 0.75, 0.015, 0.15])
    cbar = colorbar(img, cax=cbaxes, extend='both')
    # Anomalies for each experiment
    for expt in range(num_expts):
        if expt < 2:
            ax = subplot(gs_temp[0,expt+1])
        else:
            ax = subplot(gs_temp[1,expt-2])
        img = PatchCollection(patches, cmap='RdBu_r')
        img.set_array(temp_diff[expt,:])
        img.set_edgecolor('face')
        img.set_clim(vmin=-temp_max_diff, vmax=temp_max_diff)
        ax.add_collection(img)
        xlim([lat_min, lat_max])
        ylim([depth_min, 0])
        title(expt_names[expt], fontsize=20)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if expt == 0:
            xlabel(str(end_years[0])+'-'+str(end_years[1])+' anomalies', fontsize=20)
            text(-69.5, 1200, r'Temperature ($^{\circ}$C) at ' + lon_string, fontsize=28)
        if expt == num_expts-1:
            # Add a colorbar on the right
            cbaxes = fig.add_axes([0.92, 0.75, 0.015, 0.15])
            cbar = colorbar(img, cax=cbaxes, extend='both')
    # Salinity
    gs_salt = GridSpec(2,3)
    gs_salt.update(left=0.11, right=0.9, bottom=0.02, top=0.41, wspace=0.05, hspace=0.5)
    # Beginning
    ax = subplot(gs_salt[0,0])
    img = PatchCollection(patches, cmap='jet')
    img.set_array(salt_beg)
    img.set_edgecolor('face')
    img.set_clim(vmin=salt_min, vmax=salt_max)
    ax.add_collection(img)
    xlim([lat_min, lat_max])
    ylim([depth_min, 0])
    title(str(beg_years[0])+'-'+str(beg_years[1]), fontsize=20)
    xlabel('Latitude', fontsize=16)
    ylabel('Depth (m)', fontsize=16)
    # Add a colorbar on the left
    cbaxes = fig.add_axes([0.02, 0.26, 0.015, 0.15])
    cbar = colorbar(img, cax=cbaxes, extend='both')
    # Anomalies for each experiment
    for expt in range(num_expts):
        if expt < 2:
            ax = subplot(gs_salt[0,expt+1])
        else:
            ax = subplot(gs_salt[1,expt-2])
        img = PatchCollection(patches, cmap='RdBu_r')
        img.set_array(salt_diff[expt,:])
        img.set_edgecolor('face')
        img.set_clim(vmin=-salt_max_diff, vmax=salt_max_diff)
        ax.add_collection(img)
        xlim([lat_min, lat_max])
        ylim([depth_min, 0])
        title(expt_names[expt], fontsize=20)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if expt == 0:
            xlabel(str(end_years[0])+'-'+str(end_years[1])+' anomalies', fontsize=20)
            text(-68, 1200, 'Salinity (psu) at ' + lon_string, fontsize=28)
        if expt == num_expts-1:
            # Add a colorbar on the right
            cbaxes = fig.add_axes([0.92, 0.26, 0.015, 0.15])
            cbar = colorbar(img, cax=cbaxes, extend='both')
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    lon0 = float(raw_input('Longitude to plot (-180 to 180): '))
    lat_min = float(raw_input('Minimum latitude to plot (-90 to 90): '))
    lat_max = float(raw_input('Maximum latitude to plot (-90 to 90): '))
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input('Filename for figure: ')
    elif action == 'd':
        save = False
        fig_name = None
    zonal_ts_rcp_lon(lon0, lat_min, lat_max, save, fig_name)
