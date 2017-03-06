from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.pyplot import *
from matplotlib.colors import LinearSegmentedColormap
from patches import *

def massloss_percent_change ():

    directory_head = '/short/y99/kaa561/FESOM/'
    mesh_path = directory_head + 'mesh/low_res/'
    expt_paths = ['rcp45_M', 'rcp85_M'] #'rcp45_A', 'rcp85_M', 'rcp85_A']
    expt_titles = ['RCP 4.5', 'RCP 8.5'] #['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A']
    year_skip = 2005-1992+1
    peryear = 365/5

    lon_min = [-62.67, -65.5, -79.17, -85, -104.17, -102.5, -108.33, -114.5, -135.67, -149.17, -155, 144, 115, 94.17, 80.83, 65, 33.83, 19, 12.9, 9.33, -10.05, -28.33, -181, 158.33]
    lon_max = [-59.33, -60, -66.67, -28.33, -88.83, -99.17, -103.33, -111.5, -114.33, -140, -145, 146.62, 123.33, 102.5, 89.17, 75, 37.67, 33.33, 16.17, 12.88, 7.6, -10.33, -146.67, 181]
    lat_min = [-73.03, -69.35, -74.17, -83.5, -73.28, -75.5, -75.5, -75.33, -74.9, -76.42, -78, -67.83, -67.17, -66.67, -67.83, -73.67, -69.83, -71.67, -70.5, -70.75, -71.83, -76.33, -85, -84.5]
    lat_max = [-69.37, -66.13, -69.5, -74.67, -71.67, -74.17, -74.67, -73.67, -73, -75.17, -76.41, -66.67, -66.5, -64.83, -66.17, -68.33, -68.67, -68.33, -69.33, -69.83, -69.33, -71.5, -77.77, -77]
    num_shelves = len(lon_min)-1

    max_lat_plot = -63+90
    mask_cavities = True
    circumpolar = True

    elements, mask_patches = make_patches(mesh_path, circumpolar, mask_cavities, only_major=True)
    patches = iceshelf_mask(elements, only_major=True)

    values = empty([2, len(patches)]) #4, len(patches)])
    for expt in range(2): #4):
        f = open(directory_head + expt_paths[expt] + '/massloss.log', 'r')
        f.readline()
        for line in f:
            try:
                tmp = float(line)
            except(ValueError):
                break
        percent_change = empty(num_shelves)
        for index in range(num_shelves):
            massloss_tmp = []
            for line in f:
                try:
                    massloss_tmp.append(float(line))
                except(ValueError):
                    break
            massloss_first10y = mean(array(massloss_tmp[year_skip*peryear:(year_skip+10)*peryear]))
            massloss_last10y = mean(array(massloss_tmp[-10*peryear:]))
            percent_change[index] = (massloss_last10y - massloss_first10y)/massloss_first10y*100

        values_tmp = []
        for elm in elements:
            if elm.cavity:
                keep = False
                for index in range(num_shelves):
                    if all(elm.lon >= lon_min[index]) and all(elm.lon <= lon_max[index]) and all(elm.lat >= lat_min[index]) and all(elm.lat <= lat_max[index]):
                        keep = True
                        tmp = percent_change[index]
                    if index == num_shelves-1:
                        if all(elm.lon >= lon_min[index+1]) and all(elm.lon <= lon_max[index+1]) and all(elm.lat >= lat_min[index+1]) and all(elm.lat <= lat_max[index+1]):
                            keep = True
                            tmp = percent_change[index]
                if keep:
                    values_tmp.append(tmp)
        values[expt,:] = values_tmp

    x_reg, y_reg = meshgrid(linspace(-max_lat_plot, max_lat_plot, num=100), linspace(-max_lat_plot, max_lat_plot, num=100))
    land_square = zeros(shape(x_reg))

    contour_lines = []
    for elm in elements:
        if count_nonzero(elm.cavity_nodes) == 2:
            coast_tmp = []
            x_tmp = []
            y_tmp = []
            for i in range(3):
                if elm.cavity_nodes[i]:
                    coast_tmp.append(elm.coast_nodes[i])
                    x_tmp.append(elm.x[i])
                    y_tmp.append(elm.y[i])
            if count_nonzero(coast_tmp) < 2:
                contour_lines.append([(x_tmp[0], y_tmp[0]), (x_tmp[1], y_tmp[1])])

    cmap_vals = array([0, 35, 70, amax(values)])
    cmap_colours = [(1,1,1),(1,0.9,0.4),(0.99,0.59,0.18),(0.5,0.0,0.08)]
    cmap_vals_norm = cmap_vals/amax(values)
    cmap_list = []
    for i in range(size(cmap_vals)):
        cmap_list.append((cmap_vals_norm[i], cmap_colours[i]))
    my_cmap = LinearSegmentedColormap.from_list('my_cmap',cmap_list)

    fig = figure(figsize=(16,8)) #12))
    for expt in range(2): #4):
        ax = fig.add_subplot(1, 2, expt+1, aspect='equal') #2, 2, expt+1, aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(patches, cmap=my_cmap)
        img.set_array(values[expt,:])
        img.set_edgecolor('face')
        ax.add_collection(img)
        overlay = PatchCollection(mask_patches, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
        ax.add_collection(contours)
        xlim([-max_lat_plot, max_lat_plot])
        ylim([-max_lat_plot, max_lat_plot])
        axis('off')
        title(expt_titles[expt], fontsize=24)
        if expt == 0:
            cbaxes = fig.add_axes([0.06, 0.2, 0.03, 0.6]) #0.55, 0.02, 0.3])
            cbar = colorbar(img, cax=cbaxes)
            cbar.ax.tick_params(labelsize=16)
    suptitle('% Change in Ice Shelf Mass Loss (2091-2100 vs 2006-2015)', fontsize=30)
    subplots_adjust(wspace=0.025,hspace=0.1)

    #fig.show()
    fig.savefig('massloss_percent_change.png')
    

# Command-line interface
if __name__ == "__main__":

    massloss_percent_change()

    


    
