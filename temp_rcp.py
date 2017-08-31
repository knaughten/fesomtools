from numpy import *
from netCDF4 import Dataset
from matplotlib.pyplot import *

def temp_rcp ():

    directory = '/short/y99/kaa561/CMIP5_forcing_monthly/'
    beg_file = 'ERA-Interim.nc'
    files = ['rcp45_M.nc' , 'rcp45_A.nc', 'rcp85_M.nc', 'rcp85_A.nc']
    num_expts = len(files)
    expt_names = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A']
    colours = ['blue', 'cyan', 'green', 'magenta']
    lat_min = -90
    lat_max = -30

    id = Dataset(directory + beg_file, 'r')
    latitude = id.variables['latitude'][:]
    tair = mean(id.variables['t2m'][-12:,:,:], axis=0)
    id.close()
    tair_beg = mean(tair, axis=1)
    num_lat = size(latitude)
    tair_diff = zeros([num_expts, num_lat])
    for expt in range(num_expts):
        id = Dataset(directory + files[expt], 'r')
        tair = mean(id.variables['t2m'][-12:,:,:], axis=0)
        id.close()
        tair_diff[expt,:] = mean(tair, axis=1) - tair_beg

    fig, ax = subplots(figsize=(10,6))
    for expt in range(num_expts):
        ax.plot(tair_diff[expt,:], latitude, label=expt_names[expt], color=colours[expt], linewidth=2)
    title('Zonally averaged surface air temperature anomaly\n2091-2100 minus 1992-2005', fontsize=16)
    xlabel(r'$^{\circ}$C', fontsize=14)
    ylabel('latitude', fontsize=14)
    xlim([0, 6])
    ylim([lat_min, lat_max])
    grid(True)
    # Move plot over to make room for legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    # Make legend
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('temp_rcp.png')


# Command-line interface
if __name__ == "__main__":

    temp_rcp()

    
