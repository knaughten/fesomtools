from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from fesom_grid import *

def ts_distribution (mesh_path, file_path, tstep, save=False, fig_name=None):

    nbdry = -50
    num_bins = 1000
    circumpolar = False
    cross_180 = False

    id = Dataset(file_path, 'r')
    temp = id.variables['temp'][tstep-1,:]
    salt = id.variables['salt'][tstep-1,:]
    id.close()

    f = open(mesh_path + 'nod3d.out', 'r')
    f.readline()
    lat = []
    for line in f:
        tmp = line.split()
        lat.append(float(tmp[2]))
    f.close()

    min_temp = temp[0]
    max_temp = temp[0]
    min_salt = salt[0]
    max_salt = salt[0]
    for i in range(1, len(lat)):
        if lat[i] < nbdry:
            if temp[i] < min_temp:
                min_temp = temp[i]
            if temp[i] > max_temp:
                max_temp = temp[i]
            if salt[i] < min_salt:
                min_salt = salt[i]
            if salt[i] > max_salt:
                max_salt = salt[i]

    temp_bins = linspace(min_temp, max_temp, num=num_bins)
    temp_centres = 0.5*(temp_bins[:-1] + temp_bins[1:])
    salt_bins = linspace(min_salt, max_salt, num=num_bins)
    salt_centres = 0.5*(salt_bins[:-1] + salt_bins[1:])
    ts_vals = zeros([size(temp_centres), size(salt_centres)])

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
    pcolor(salt_centres, temp_centres, log(ts_vals), cmap='jet')
    xlim([31.8, 35.2])
    ylim([-3, 12])
    xlabel('Salinity (psu)', fontsize=16)
    ylabel(r'Temperature ($^{\circ}$C)', fontsize=16)
    title('Water masses south of ' + str(-nbdry) + r'$^{\circ}$S, log(volume)', fontsize=24)
    colorbar()

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
    
