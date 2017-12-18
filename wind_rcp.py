from numpy import *
from netCDF4 import Dataset
from matplotlib.pyplot import *

def wind_rcp ():

    directory = '/short/y99/kaa561/CMIP5_forcing_monthly/'
    files = ['ERA-Interim.nc', 'rcp45_M.nc' , 'rcp45_A.nc', 'rcp85_M.nc', 'rcp85_A.nc']
    num_expts = len(files)
    expt_names = ['1992-2005', 'RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A']
    colours = ['black', 'blue', 'cyan', 'green', 'magenta']
    lat_min = -90
    lat_max = -30

    id = Dataset(directory + files[0], 'r')
    latitude = id.variables['latitude'][:]
    id.close()
    num_lat = size(latitude)
    uwind_avg = zeros([num_expts, num_lat])
    for expt in range(num_expts):
        id = Dataset(directory + files[expt], 'r')
        uwind = mean(id.variables['u10'][-12:,:,:], axis=0)
        uwind_avg[expt,:] = mean(uwind, axis=1)
        id.close()

    fig, ax = subplots(figsize=(10,6))
    for expt in range(num_expts):
        ax.plot(uwind_avg[expt,:], latitude, label=expt_names[expt], color=colours[expt], linewidth=2)
    title('Zonally averaged eastward wind (2091-2100)', fontsize=18)
    xlabel('m/s', fontsize=14)
    ylabel('latitude', fontsize=14)
    ylim([lat_min, lat_max])
    grid(True)
    # Move plot over to make room for legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    # Make legend
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.show()
    fig.savefig('wind_rcp.png')


# Command-line interface
if __name__ == "__main__":

    wind_rcp()

    
