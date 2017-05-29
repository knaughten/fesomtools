from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from fesom_grid import *
from unesco import *

def ts_distribution (mesh_path, file_path, tstep, save=False, fig_name=None):

    nbdry = -50
    num_bins = 1000
    circumpolar = False
    cross_180 = False
    min_salt = 31.8
    max_salt = 35.2
    min_temp = -3
    max_temp = 12

    id = Dataset(file_path, 'r')
    temp = id.variables['temp'][tstep-1,:]
    salt = id.variables['salt'][tstep-1,:]
    id.close()

    temp_bins = linspace(min_temp, max_temp, num=num_bins)
    temp_centres = 0.5*(temp_bins[:-1] + temp_bins[1:])
    salt_bins = linspace(min_salt, max_salt, num=num_bins)
    salt_centres = 0.5*(salt_bins[:-1] + salt_bins[1:])
    ts_vals = zeros([size(temp_centres), size(salt_centres)])

    freezing_pt = -0.0575*salt_centres + 1.7105e-3*sqrt(salt_centres**3) - 2.155e-4*salt_centres**2
    salt_2d, temp_2d = meshgrid(salt_centres, temp_centres)
    density = unesco(temp_2d, salt_2d, zeros(shape(temp_centres)))-1000
    density_lev = arange(24.4, 28.4, 0.2)

    elements = fesom_grid(mesh_path, circumpolar, cross_180)

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
    img=pcolor(salt_centres, temp_centres, log(ts_vals), cmap='jet')
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

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    file_path = raw_input("Path to FESOM oce.mean.nc file: ")
    tstep = int(raw_input("Time index to plot (starting at 1): "))
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None

    ts_distribution(mesh_path, file_path, tstep, save, fig_name)
    
