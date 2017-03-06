from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *

def zonal_wind_scenarios ():

    era_dir = '/short/y99/kaa561/FESOM/ERA_Interim/'
    rcp_dir_head = '/short/y99/kaa561/FESOM/RCP_forcing/'
    expt_dir = ['rcp45/MMM/', 'rcp85/MMM/', 'rcp45/ACCESS1-3/', 'rcp85/ACCESS1-3/']
    file_head = 'uwind_'
    file_tails = ['00.nc', '06.nc', '12.nc', '18.nc']
    era_start_year = 1992
    era_end_year = 2005
    rcp_start_year = [2006, 2091]
    rcp_end_year = [2015, 2100]
    era_colour = 'black'
    era_title = 'ERA-Interim (1992-2005)'
    rcp_colours = ['red', 'blue', 'magenta', 'green']
    rcp_titles = ['RCP 4.5', 'RCP 8.5']
    year_titles = [' (2006-2015)', ' (2091-2100)']
    lat_start = 159

    id = Dataset(era_dir + file_head + file_tails[0], 'r')
    lon = id.variables['longitude'][:]
    lat = id.variables['latitude'][lat_start:]
    id.close()

    era_start = (era_start_year-1979)*365 + int(floor((era_start_year-1977)/4.0))
    era_end = (era_end_year+1-1979)*365 + int(floor((era_end_year+1-1977)/4.0))
    rcp_start1 = (rcp_start_year[0]-2006)*365 + int(floor((rcp_start_year[0]-2005)/4.0))
    rcp_end1 = (rcp_end_year[0]+1-2006)*365 + int(floor((rcp_end_year[0]+1-2005)/4.0))
    rcp_start2 = (rcp_start_year[1]-2006)*365 + int(floor((rcp_start_year[1]-2005)/4.0))
    rcp_end2 = (rcp_end_year[1]+1-2006)*365 + int(floor((rcp_end_year[1]+1-2005)/4.0))

    print 'Reading ERA-Interim data'
    era_wind = zeros([size(lat), size(lon)])
    for td in range(4):
        id = Dataset(era_dir + file_head + file_tails[td], 'r')
        era_wind[:,:] = era_wind[:,:] + mean(id.variables['u10'][era_start:era_end,lat_start:,:], axis=0)
    era_wind[:,:] = era_wind[:,:]/4.0
    era_wind_zonalavg = mean(era_wind, axis=1)

    print 'Reading RCP data, 2006-2015'
    rcp_wind1 = zeros([4, size(lat), size(lon)])
    for expt in range(4):
        print '...experiment ' + str(expt+1) + ' of 4'
        for td in range(4):
            id = Dataset(rcp_dir_head + expt_dir[expt] + file_head + file_tails[td], 'r')
            rcp_wind1[expt,:,:] = rcp_wind1[expt,:,:] + mean(id.variables['u10'][rcp_start1:rcp_end1,lat_start:,:], axis=0)
        rcp_wind1[expt,:,:] = rcp_wind1[expt,:,:]/4.0
    rcp_wind1_zonalavg = mean(rcp_wind1, axis=2)

    print 'Reading RCP data, 2091-2100'
    rcp_wind2 = zeros([4, size(lat), size(lon)])
    for expt in range(4):
        print '...experiment ' + str(expt+1) + ' of 4'
        for td in range(4):
            id = Dataset(rcp_dir_head + expt_dir[expt] + file_head + file_tails[td], 'r')
            rcp_wind2[expt,:,:] = rcp_wind2[expt,:,:] + mean(id.variables['u10'][rcp_start2:rcp_end2,lat_start:,:], axis=0)
        rcp_wind2[expt,:,:] = rcp_wind2[expt,:,:]/4.0
    rcp_wind2_zonalavg = mean(rcp_wind2, axis=2)

    fig1, ax1 = subplots(figsize=(10,6))
    ax1.plot(era_wind_zonalavg, lat, label=era_title, color=era_colour, linewidth=2, linestyle='solid')
    for expt in range(2):
        ax1.plot(rcp_wind1_zonalavg[expt,:], lat, label=rcp_titles[expt] + year_titles[0], color=rcp_colours[expt], linewidth=2, linestyle='solid')
        ax1.plot(rcp_wind2_zonalavg[expt,:], lat, label=rcp_titles[expt] + year_titles[1], color=rcp_colours[expt], linewidth=2, linestyle='dashed')
    title('Zonally averaged u-wind (MMM)', fontsize=18)
    xlabel('m/s', fontsize=14)
    ylabel('Latitude', fontsize=14)
    xlim([-3, 8.5])
    ylim([amin(lat), amax(lat)])
    grid(True)
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width*0.7, box.height])
    ax1.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig1.savefig('uwind_m.png')

    fig2, ax2 = subplots(figsize=(10,6))
    ax2.plot(era_wind_zonalavg, lat, label=era_title, color=era_colour, linewidth=2, linestyle='solid')
    for expt in range(2,4):
        ax2.plot(rcp_wind1_zonalavg[expt,:], lat, label=rcp_titles[expt-2] + year_titles[0], color=rcp_colours[expt], linewidth=2, linestyle='solid')
        ax2.plot(rcp_wind2_zonalavg[expt,:], lat, label=rcp_titles[expt-2] + year_titles[1], color=rcp_colours[expt], linewidth=2, linestyle='dashed')
    title('Zonally averaged u-wind (ACCESS1-3)', fontsize=18)
    xlabel('m/s', fontsize=14)
    ylabel('Latitude', fontsize=14)
    xlim([-3, 8.5])
    ylim([amin(lat), amax(lat)])
    grid(True)
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0, box.width*0.7, box.height])
    ax2.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig2.savefig('uwind_a.png')


if __name__ == "__main__":

    zonal_wind_scenarios()
    
    

    
                

    

    
    
    
