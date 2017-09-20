from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.colors import ListedColormap

def rtopo_plots ():

    data_file = "../FESOM_mesh/RTopo105/RTopo105_data.nc"
    aux_file = "../FESOM_mesh/RTopo105/RTopo105_aux.nc"
    deg2rad = pi/180
    nbdry = -60 + 90
    max_j = 2000

    id = Dataset(data_file, 'r')
    lon = id.variables['lon'][:]
    lat = id.variables['lat'][:max_j]
    bathy = -1*id.variables['bathy'][:max_j,:]
    draft = -1*id.variables['draft'][:max_j,:]
    id.close()

    id = Dataset(aux_file, 'r')
    amask = id.variables['amask'][:max_j,:]
    id.close()

    omask = zeros(shape(amask))
    index = amask == 0
    omask[index] = 1
    index = amask == 2
    omask[index] = 1
    imask = zeros(shape(amask))
    index = amask == 2
    imask[index] = 1
    land_zice = ma.masked_where(amask==0, amask)

    bathy = ma.masked_where(omask==0, bathy)
    draft = ma.masked_where(imask==0, draft)

    lon_2d, lat_2d = meshgrid(lon, lat)
    x = -(lat_2d+90)*cos(lon_2d*deg2rad+pi/2)
    y = (lat_2d+90)*sin(lon_2d*deg2rad+pi/2)

    fig = figure(figsize=(128,96))
    fig.add_subplot(1,1,1, aspect='equal')
    img = pcolor(x, y, bathy, vmin=0, vmax=2800)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    cbar = colorbar(img, extend='max')
    cbar.ax.tick_params(labelsize=160)
    title('Bathymetry (m)', fontsize=240)
    axis('off')
    fig.savefig('rtopo_depth.png')

    fig = figure(figsize=(128,96))
    fig.add_subplot(1,1,1, aspect='equal')
    fig.patch.set_facecolor('white')
    grey_cmap = ListedColormap([(0.6, 0.6, 0.6)])
    pcolor(x, y, land_zice, cmap=grey_cmap)    
    img = pcolor(x, y, draft, vmin=0, vmax=2300)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    cbar = colorbar(img, extend='max')
    cbar.ax.tick_params(labelsize=160)
    title('Ice shelf draft (m)', fontsize=240)
    axis('off')
    fig.savefig('rtopo_shelf.png')


if __name__ == "__main__":

    rtopo_plots()

    
    
    
    
