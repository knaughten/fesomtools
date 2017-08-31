from netCDF4 import Dataset
from numpy import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from fesom_grid import *
from fesom_sidegrid import *
from temp_salt_seasonal import interp_lon_fesom 

def amundsen_slices (rcp, model, fig_name):

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/high_res/'
    file_beg = '/short/y99/kaa561/FESOM/highres_spinup/seasonal_climatology_oce_1996_2005.nc'
    file_end = '/short/y99/kaa561/FESOM/rcp'+rcp+'_'+model+'_highres/output/seasonal_climatology_oce_2091_2100.nc'
    # Longitude to interpolate to
    lon0 = -104
    # Latitude range to plot
    lat_min = -76
    lat_max = -70
    # Deepest depth to plot
    depth_min = -1000
    # Bounds on colour scales for temperature and salinity
    temp_min = -2
    temp_max = 2
    temp_tick = 1
    salt_min = 33.6
    salt_max = 34.8
    salt_tick = 0.2
    # Season names for titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    # Longitude format for title
    lon_string = str(int(round(-lon0))) + r'$^{\circ}$W'
    # RCP format for title
    if rcp == '45':
        rcp_title = 'RCP 4.5 ' + model
    elif rcp == '85':
        rcp_title = 'RCP 8.5 ' + model

    print 'Building FESOM mesh'
    elements = fesom_grid(mesh_path)

    print 'Reading data'
    id = Dataset(file_beg, 'r')
    temp_nodes_beg = id.variables['temp'][:,:]
    salt_nodes_beg = id.variables['salt'][:,:]
    id.close()
    id = Dataset(file_end, 'r')
    temp_nodes_end = id.variables['temp'][:,:]
    salt_nodes_end = id.variables['salt'][:,:]

    # Set up plot
    fig = figure(figsize=(20,18))
    gs_temp = GridSpec(2,4)
    gs_temp.update(left=0.05, right=0.9, bottom=0.54, top=0.91, wspace=0.05, hspace=0.2)
    gs_salt = GridSpec(2,4)
    gs_salt.update(left=0.05, right=0.9, bottom=0.06, top=0.44, wspace=0.05, hspace=0.2)
    print 'Temperature, 1996-2005'
    for season in range(4):
        print '...' + season_names[season]
        # Interpolate to lon0
        patches, values, tmp = interp_lon_fesom(elements, lat_max, lon0, temp_nodes_beg[season,:])
        ax = subplot(gs_temp[0,season])
        img = PatchCollection(patches, cmap='jet')
        img.set_array(array(values))
        img.set_edgecolor('face')
        img.set_clim(vmin=temp_min, vmax=temp_max)
        ax.add_collection(img)
        xlim([lat_min, lat_max])
        ylim([depth_min, 0])
        if season == 0:
            xlabel('Latitude', fontsize=14)
            ylabel('Depth (m)', fontsize=14)
        else:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        if season == 1:
            # Titles
            text(-69, 250, r'Temperature ($^{\circ}$C) at ' + lon_string, ha='center', fontsize=30)
            text(-69, 50, '1996-2005', ha='center', fontsize=24)
        if season == 3:
            # Colourbar
            cbaxes = fig.add_axes([0.93, 0.6, 0.015, 0.25])
            cbar = colorbar(img, cax=cbaxes, ticks=arange(temp_min, temp_max+temp_tick, temp_tick), extend='both')
            cbar.ax.tick_params(labelsize=16)
    print 'Temperature, 2091-2100'
    for season in range(4):
        print '...' + season_names[season]
        patches, values, tmp = interp_lon_fesom(elements, lat_max, lon0, temp_nodes_end[season,:])
        ax = subplot(gs_temp[1,season])
        img = PatchCollection(patches, cmap='jet')
        img.set_array(array(values))
        img.set_edgecolor('face')
        img.set_clim(vmin=temp_min, vmax=temp_max)
        ax.add_collection(img)
        xlim([lat_min, lat_max])
        ylim([depth_min, 0])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if season == 1:
            # Title
            text(-69, 50, rcp_title, ha='center', fontsize=24)
    print 'Salinity, 1996-2005'
    for season in range(4):
        print '...' + season_names[season]
        patches, values, tmp = interp_lon_fesom(elements, lat_max, lon0, salt_nodes_beg[season,:])
        ax = subplot(gs_salt[0,season])
        img = PatchCollection(patches, cmap='jet')
        img.set_array(array(values))
        img.set_edgecolor('face')
        img.set_clim(vmin=salt_min, vmax=salt_max)
        ax.add_collection(img)
        xlim([lat_min, lat_max])
        ylim([depth_min, 0])
        if season == 0:
            xlabel('Latitude', fontsize=14)
            ylabel('Depth (m)', fontsize=14)
        else:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        if season == 1:
            # Titles
            text(-69, 250, 'Salinity (psu) at ' + lon_string, ha='center', fontsize=30)
            text(-69, 50, '1996-2005', ha='center', fontsize=24)
        if season == 3:
            # Colourbar
            cbaxes = fig.add_axes([0.93, 0.13, 0.015, 0.25])
            cbar = colorbar(img, cax=cbaxes, ticks=arange(salt_min, salt_max+salt_tick, salt_tick), extend='both')
            cbar.ax.tick_params(labelsize=16)
    print 'Salinity, 2091-2100'
    for season in range(4):
        print '...' + season_names[season]
        patches, values, tmp = interp_lon_fesom(elements, lat_max, lon0, salt_nodes_end[season,:])
        ax = subplot(gs_salt[1,season])
        img = PatchCollection(patches, cmap='jet')
        img.set_array(array(values))
        img.set_edgecolor('face')
        img.set_clim(vmin=salt_min, vmax=salt_max)
        ax.add_collection(img)
        xlim([lat_min, lat_max])
        ylim([depth_min, 0])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        xlabel(season_names[season], fontsize=26)
        if season == 1:
            # Title
            text(-69, 50, rcp_title, ha='center', fontsize=24)
    fig.savefig(fig_name)


# Command-line interface
if __name__ == "__main__":

    tmp = int(raw_input("RCP 4.5 (1) or 8.5 (2)? "))
    if tmp == 1:
        rcp = '45'
    elif tmp == 2:
        rcp = '85'
    tmp = int(raw_input("Multi-model mean (1) or ACCESS 1.0 (2)? "))
    if tmp == 1:
        model = 'M'
    elif tmp == 2:
        model = 'A'
    fig_name = raw_input("Filename for figure: ")
    amundsen_slices(rcp, model, fig_name)
    
            
        
        
    
    
    
