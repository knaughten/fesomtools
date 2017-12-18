from numpy import *
from netCDF4 import Dataset
from matplotlib.pyplot import *

def wind_curl_rcp ():

    directory = '/short/y99/kaa561/CMIP5_forcing_monthly/'
    files = ['ERA-Interim.nc', 'rcp45_M.nc' , 'rcp45_A.nc', 'rcp85_M.nc', 'rcp85_A.nc']
    num_expts = len(files)
    expt_names = ['1992-2005', 'RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A']
    file_tails = ['tmp', 'rcp45_m', 'rcp45_a', 'rcp85_m', 'rcp85_a']
    colours = ['black', 'blue', 'cyan', 'green', 'magenta']
    lat_min = -90
    lat_max = -30
    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0

    # Read ERA-Interim grid
    id = Dataset(directory + files[0], 'r')
    longitude = id.variables['longitude'][:]
    latitude = id.variables['latitude'][1:-1]  # Ignore the first and last values where cos(lat)=0
    id.close()
    num_lat = size(latitude)
    num_lon = size(longitude)
    # Tile latitude into a 2D array
    lat_2d = transpose(tile(latitude, (num_lon,1)))
    # Constant grid spacing in degrees
    dlon = longitude[1] - longitude[0]
    dlat = latitude[1] - latitude[0]
    # Cartesian differentials
    dx = r*cos(lat_2d*deg2rad)*dlon*deg2rad  # 2D array
    dy = r*dlat*deg2rad  # Scalar

    # Calculate curl of wind velocity in each experiment
    wind_curl = zeros([num_expts, num_lat, num_lon])
    for expt in range(num_expts):
        id = Dataset(directory + files[expt], 'r')
        uwind = mean(id.variables['u10'][-12:,1:-1,:], axis=0)
        vwind = mean(id.variables['v10'][-12:,1:-1,:], axis=0)
        id.close()
        # Calculate the two derivatives
        dv_dx = zeros(shape(uwind))
        du_dy = zeros(shape(uwind))
        # Forward difference approximation
        dv_dx[:,:-1] = (vwind[:,1:] - vwind[:,:-1])/dx[:,:-1]
        du_dy[:-1,:] = (uwind[1:,:] - uwind[:-1,:])/dy
        # Backward difference for the last row
        dv_dx[:,-1] = (vwind[:,-1] - vwind[:,-2])/dx[:,-1]
        du_dy[-1,:] = (vwind[-1,:] - vwind[-2,:])/dy
        wind_curl[expt,:,:] = dv_dx - du_dy
    # Convert to 1e-6 s-1 for readability
    wind_curl *= 1e6

    # Calculate zonal averages
    wind_curl_avg = mean(wind_curl, axis=2)
    # Plot zonal averages
    fig, ax = subplots(figsize=(10,6))
    for expt in range(num_expts):
        ax.plot(wind_curl_avg[expt,:], latitude, label=expt_names[expt], color=colours[expt], linewidth=2)
    title('Curl of wind velocity (2091-2100)', fontsize=18)
    xlabel(r'10$^{-6}$ s$^{-1}$', fontsize=14)
    ylabel('latitude', fontsize=14)
    ylim([lat_min, lat_max])
    grid(True)
    # Move plot over to make room for legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    # Make legend
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.show()
    fig.savefig('wind_curl_rcp.png')

    # Plot anomalies in zonal averages
    fig, ax = subplots(figsize=(10,6))
    for expt in range(1,num_expts):
        ax.plot(wind_curl_avg[expt,:]-wind_curl_avg[0,:], latitude, label=expt_names[expt], color=colours[expt], linewidth=2)
    title('Anomalies in curl of wind velocity (2091-2100 minus 1992-2005)', fontsize=18)
    xlabel(r'10$^{-6}$ s$^{-1}$', fontsize=14)
    ylabel('latitude', fontsize=14)
    ylim([lat_min, lat_max])
    grid(True)
    # Move plot over to make room for legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    # Make legend
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.show()
    fig.savefig('wind_curl_diff_rcp.png')

    # Plot initial 2D wind curl field
    bound = 100
    lev = linspace(-bound, bound, num=50)
    fig, ax = subplots(figsize=(10,6))
    img = ax.contourf(longitude, latitude, wind_curl[0,:,:], lev, cmap='RdBu_r', extend='both')
    xlabel('Longitude')
    ylabel('Latitude')
    ylim([lat_min, lat_max])
    title(r'Curl of wind velocity, 1992-2005 (s$^{-1}$)', fontsize=18)
    colorbar(img)
    fig.show()
    fig.savefig('wind_curl_2d.png')

    # Plot anomalies in 2D wind curl field
    for expt in range(1, num_expts):
        wind_curl_diff = wind_curl[expt,:,:] - wind_curl[0,:,:]
        bound = 10
        lev = linspace(-bound, bound, num=50)
        fig, ax = subplots(figsize=(10,6))
        img = ax.contourf(longitude, latitude, wind_curl_diff, lev, cmap='RdBu_r', extend='both')
        xlabel('Longitude')
        ylabel('Latitude')
        ylim([lat_min, lat_max])
        title('Anomalies in curl of wind velocity, ' + expt_names[expt], fontsize=18)
        colorbar(img)
        fig.show()
        fig.savefig('wind_curl_diff_2d_' + file_tails[expt] + '.png')


# Command-line interface
if __name__ == "__main__":

    wind_curl_rcp ()
        
        
