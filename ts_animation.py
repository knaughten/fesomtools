from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from fesom_grid import *
from unesco import *

def ts_animation (mesh_path, directory, start_year, end_year, fig_dir):

    nbdry = -50
    num_bins = 1000
    circumpolar = False
    cross_180 = False
    min_salt = 31.8
    max_salt = 35.2
    min_temp = -3
    max_temp = 12
    min_vol = 18
    max_vol = 33
    file_head = 'MK44005.'
    file_tail = '.oce.mean.nc'

    temp_bins = linspace(min_temp, max_temp, num=num_bins)
    temp_centres = 0.5*(temp_bins[:-1] + temp_bins[1:])
    salt_bins = linspace(min_salt, max_salt, num=num_bins)
    salt_centres = 0.5*(salt_bins[:-1] + salt_bins[1:])

    freezing_pt = -0.0575*salt_centres + 1.7105e-3*sqrt(salt_centres**3) - 2.155e-4*salt_centres**2
    salt_2d, temp_2d = meshgrid(salt_centres, temp_centres)
    density = unesco(temp_2d, salt_2d, zeros(shape(temp_centres)))-1000
    density_lev = arange(24.4, 28.4, 0.2)

    elements = fesom_grid(mesh_path, circumpolar, cross_180)

    for year in range(start_year, end_year+1):
        print 'Processing ' + str(year)
        id = Dataset(directory + file_head + str(year) + file_tail, 'r')
        temp = mean(id.variables['temp'][:,:], axis=0)
        salt = mean(id.variables['salt'][:,:], axis=0)
        id.close()
        ts_vals = zeros([size(temp_centres), size(salt_centres)])
        for elm in elements:
            if all(elm.lat < nbdry):
                area = elm.area()
                nodes = [elm.nodes[0], elm.nodes[1], elm.nodes[2]]
                while True:
                    if nodes[0].below is None or nodes[1].below is None or nodes[2].below is None:
                        break
                    temp_vals = []
                    salt_vals = []
                    dz = []
                    for i in range(3):
                        temp_vals.append(temp[nodes[i].id])
                        temp_vals.append(temp[nodes[i].below.id])
                        salt_vals.append(salt[nodes[i].id])
                        salt_vals.append(salt[nodes[i].below.id])
                        dz.append(abs(nodes[i].depth - nodes[i].below.depth))
                        nodes[i] = nodes[i].below
                    temp_elm = mean(array(temp_vals))
                    salt_elm = mean(array(salt_vals))
                    volume = area*mean(array(dz))
                    temp_index = nonzero(temp_bins > temp_elm)[0][0] - 1
                    salt_index = nonzero(salt_bins > salt_elm)[0][0] - 1
                    ts_vals[temp_index, salt_index] += volume
        ts_vals = ma.masked_where(ts_vals==0, ts_vals)
        fig = figure(figsize=(12,12))
        img=pcolor(salt_centres, temp_centres, log(ts_vals), vmin=min_vol, vmax=max_vol, cmap='jet')
        plot(salt_centres, freezing_pt, color='black', linestyle='dashed')
        cs=contour(salt_centres, temp_centres, density, density_lev, colors=(0.6,0.6,0.6), linestyles='dotted')
        manual_locations = [(32,11.4),(32.3,11.4),(32.5,11.4),(32.8,11.4),(33.1,11.4),(33.3,11.4),(33.5,11.3),(33.8,11.3),(34.1,11.3),(34.3,11.3),(34.6,11.3),(34.8,11.4),(35,10.8),(35,9.9),(35,8.1),(35,7.5),(35,6),(35,4.4),(35,2.6),(35.1,0)]
        clabel(cs, inline=1, fontsize=12, color=(0.6,0.6,0.6), fmt='%1.1f', manual=manual_locations)
        xlim([min_salt, max_salt])
        ylim([min_temp, max_temp])
        xlabel('Salinity (psu)', fontsize=16)
        ylabel(r'Temperature ($^{\circ}$C)', fontsize=16)
        title('Water masses south of ' + str(-nbdry) + r'$^{\circ}$S, log(volume)', fontsize=24)
        colorbar(img)
        text(35.8, -4, str(year), fontsize=30)

        fig.savefig(fig_dir + str(year) + '.png')


if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    directory = raw_input("Path to FESOM output directory: ")
    start_year = int(raw_input("First year: "))
    end_year = int(raw_input("Last year: "))
    fig_dir = raw_input("Directory to store figures: ")

    ts_animation(mesh_path, directory, start_year, end_year, fig_dir)
    
