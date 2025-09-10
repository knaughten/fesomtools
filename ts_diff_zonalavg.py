from numpy import *
from netCDF4 import Dataset
from matplotlib.pyplot import *
from matplotlib.cm import *

def ts_diff_zonalavg ():

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/high_res/'
    directory_beg = '/short/y99/kaa561/FESOM/highres_spinup/'
    directories = ['/short/y99/kaa561/FESOM/rcp45_M_highres/output/', '/short/y99/kaa561/FESOM/rcp45_A_highres/output/', '/short/y99/kaa561/FESOM/rcp85_M_highres/output/', '/short/y99/kaa561/FESOM/rcp85_A_highres/output/', '/short/y99/kaa561/FESOM/highres_spinup/']
    file_beg = 'zonal_avg.1996.2005.nc'
    file_end = 'zonal_avg.2091.2100.nc'
    # Titles for plotting
    expt_names = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A', 'CONTROL']
    num_expts = len(directories)
    # Start and end years for each period
    beg_years = [1996, 2005]
    end_years = [2091, 2100]
    # Bounds on colour scale
    temp_min = -2.5
    temp_max = 7.5
    salt_min = 33.8
    salt_max = 34.8
    temp_max_diff = 3
    salt_max_diff = 0.4

    # Read grid and data at the beginning
    id = Dataset(directory_beg + file_beg, 'r')
    lat_vals = id.variables['latitude'][:]
    depth_vals = id.variables['depth'][:]
    temp_beg = id.variables['temp'][:,:]
    salt_beg = id.variables['salt'][:,:]
    id.close()

    # Get anomalies for each experiment
    temp_diff = ma.empty([num_expts, size(temp_beg,0), size(temp_beg,1)])
    salt_diff = ma.empty([num_expts, size(salt_beg,0), size(salt_beg,1)])
    for expt in range(num_expts):
        id = Dataset(directories[expt] + file_end, 'r')
        temp_diff[expt,:] = id.variables['temp'][:,:] - temp_beg
        salt_diff[expt,:] = id.variables['salt'][:,:] - salt_beg
        id.close()

    print('Plotting')
    fig = figure(figsize=(20,16))
    # Temperature
    gs_temp = GridSpec(2,3)
    gs_temp.update(left=0.11, right=0.9, bottom=0.51, top=0.9, wspace=0.05, hspace=0.5)
    # Beginning
    ax = subplot(gs_temp[0,0])
    img = ax.pcolorfast(lat_vals, depth_vals, temp_beg, cmap='jet')
    img.set_clim(vmin=temp_min, vmax=temp_max)
    xlim([amin(lat_vals), amax(lat_vals)])
    ylim([amin(depth_vals), amax(depth_vals)])
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
        img = ax.pcolorfast(lat_vals, depth_vals, temp_diff[expt], cmap='RdBu_r')
        img.set_clim(vmin=-temp_max_diff, vmax=temp_max_diff)
        xlim([amin(lat_vals), amax(lat_vals)])
        ylim([amin(depth_vals), amax(depth_vals)])
        title(expt_names[expt], fontsize=20)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if expt == 0:
            xlabel(str(end_years[0])+'-'+str(end_years[1])+' anomalies', fontsize=20)
            text(-100, 1800, r'Temperature ($^{\circ}$C), zonally averaged', fontsize=28)
        if expt == num_expts-1:
            # Add a colorbar on the right
            cbaxes = fig.add_axes([0.92, 0.75, 0.015, 0.15])
            cbar = colorbar(img, cax=cbaxes, extend='both')
    # Salinity
    gs_salt = GridSpec(2,3)
    gs_salt.update(left=0.11, right=0.9, bottom=0.02, top=0.41, wspace=0.05, hspace=0.5)
    # Beginning
    ax = subplot(gs_salt[0,0])
    img = ax.pcolorfast(lat_vals, depth_vals, salt_beg, cmap='jet')
    img.set_clim(vmin=salt_min, vmax=salt_max)
    xlim([amin(lat_vals), amax(lat_vals)])
    ylim([amin(depth_vals), amax(depth_vals)])
    title(str(beg_years[0])+'-'+str(beg_years[1]), fontsize=20)
    xlabel('Latitude', fontsize=16)
    ylabel('Depth', fontsize=16)
    # Add a colorbar on the left
    cbaxes = fig.add_axes([0.02, 0.26, 0.015, 0.15])
    cbar = colorbar(img, cax=cbaxes, extend='both')
    # Anomalies for each experiment
    for expt in range(num_expts):
        if expt < 2:
            ax = subplot(gs_salt[0,expt+1])
        else:
            ax = subplot(gs_salt[1,expt-2])
        img = ax.pcolorfast(lat_vals, depth_vals, salt_diff[expt], cmap='RdBu_r')
        img.set_clim(vmin=-salt_max_diff, vmax=salt_max_diff)
        xlim([amin(lat_vals), amax(lat_vals)])
        ylim([amin(depth_vals), amax(depth_vals)])
        title(expt_names[expt], fontsize=20)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if expt == 0:
            xlabel(str(end_years[0])+'-'+str(end_years[1])+' anomalies', fontsize=20)
            text(-95, 1800, 'Salinity (psu), zonally averaged', fontsize=28)
        if expt == num_expts-1:
            # Add a colorbar on the right
            cbaxes = fig.add_axes([0.92, 0.26, 0.015, 0.15])
            cbar = colorbar(img, cax=cbaxes, extend='both')
    fig.show()
    fig.savefig('ts_diff_zonalavg.png')
        

# Command-line interface
if __name__ == "__main__":

    ts_diff_zonalavg()
    
    
        

    

    
