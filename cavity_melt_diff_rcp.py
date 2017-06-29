from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.colors import LinearSegmentedColormap
from fesom_grid import *
from patches import *
from unrotate_vector import *
from unrotate_grid import *

# For each major ice shelf, plot the change in ice shelf melt rate during the
# RCP (last 10 years minus first 10 years) zoomed into the region of interest.
# Input:
# mesh_path = path to FESOM mesh directory
# output_path = path to output directory for that RCP, containing the filenames
#               specified in file_name_beg and file_name_end (either
#               forcing.diag.nc or oce.mean.nc averaged over the first and last
#               10 years of the RCP)
# fig_dir = path to directory for saving figures (must end in '/' unless it is
#           the current directory, i.e. empty string)
def cavity_melt_diff_rcp (mesh_path, output_path, fig_dir=''):

    file_name_beg = 'annual_avg.forcing.diag.2006.2015.nc'
    file_name_end = 'annual_avg.forcing.diag.2091.2100.nc'

    # Name of each ice shelf
    shelf_names = ['Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    # Beginnings of filenames for figures
    fig_heads = ['larsen_d', 'larsen_c', 'wilkins_georgevi_stange', 'ronne_filchner', 'abbot', 'pig', 'thwaites', 'dotson', 'getz', 'nickerson', 'sulzberger', 'mertz', 'totten_moscowuni', 'shackleton', 'west', 'amery', 'prince_harald', 'baudouin_borchgrevink', 'lazarev', 'nivl', 'fimbul_jelbart_ekstrom', 'brunt_riiser_larsen', 'ross']
    # Limits on longitude and latitude for each ice shelf
    # Note Ross crosses 180W=180E
    lon_min = [-62.67, -65.5, -79.17, -85, -104.17, -102.5, -108.33, -114.5, -135.67, -149.17, -155, 144, 115, 94.17, 80.83, 65, 33.83, 19, 12.9, 9.33, -10.05, -28.33, 158.33]
    lon_max = [-59.33, -60, -66.67, -28.33, -88.83, -99.17, -103.33, -111.5, -114.33, -140, -145, 146.62, 123.33, 102.5, 89.17, 75, 37.67, 33.33, 16.17, 12.88, 7.6, -10.33, -146.67]
    lat_min = [-73.03, -69.35, -74.17, -83.5, -73.28, -75.5, -75.5, -75.33, -74.9, -76.42, -78, -67.83, -67.17, -66.67, -67.83, -73.67, -69.83, -71.67, -70.5, -70.75, -71.83, -76.33, -85]
    lat_max = [-69.37, -66.13, -69.5, -74.67, -71.67, -74.17, -74.67, -73.67, -73, -75.17, -76.41, -66.67, -66.5, -64.83, -66.17, -68.33, -68.67, -68.33, -69.33, -69.83, -69.33, -71.5, -77]
    num_shelves = len(shelf_names)

    # Constants
    sec_per_year = 365*24*3600
    deg2rad = pi/180.0

    print 'Building FESOM mesh'
    # Mask open ocean
    elements, mask_patches = make_patches(mesh_path, circumpolar=True, mask_cavities=True)
    # Unmask ice shelves
    patches = iceshelf_mask(elements)
    print 'Reading data'
    id_beg = Dataset(output_path + file_name_beg, 'r')
    node_data_beg = id_beg.variables['wnet'][0,:]*sec_per_year
    id_beg.close()
    id_end = Dataset(output_path + file_name_end, 'r')    
    node_data_end = id_end.variables['wnet'][0,:]*sec_per_year    
    id_end.close()
    # Calculate the difference
    node_data = node_data_end - node_data_beg
    # Calculate value at each element
    data = []
    for elm in elements:
        # For each element in an ice shelf cavity, append the mean value
        # for the 3 component Nodes
        if elm.cavity:
            data.append(mean([node_data[elm.nodes[0].id], node_data[elm.nodes[1].id], node_data[elm.nodes[2].id]]))

    # Loop over ice shelves
    for index in range(num_shelves):
        print 'Processing ' + shelf_names[index]
        # Convert lat/lon bounds to polar coordinates for plotting
        x1 = -(lat_min[index]+90)*cos(lon_min[index]*deg2rad+pi/2)
        y1 = (lat_min[index]+90)*sin(lon_min[index]*deg2rad+pi/2)
        x2 = -(lat_min[index]+90)*cos(lon_max[index]*deg2rad+pi/2)
        y2 = (lat_min[index]+90)*sin(lon_max[index]*deg2rad+pi/2)
        x3 = -(lat_max[index]+90)*cos(lon_min[index]*deg2rad+pi/2)
        y3 = (lat_max[index]+90)*sin(lon_min[index]*deg2rad+pi/2)
        x4 = -(lat_max[index]+90)*cos(lon_max[index]*deg2rad+pi/2)
        y4 = (lat_max[index]+90)*sin(lon_max[index]*deg2rad+pi/2)
        # Find the new bounds on x and y
        x_min = amin(array([x1, x2, x3, x4]))
        x_max = amax(array([x1, x2, x3, x4]))
        y_min = amin(array([y1, y2, y3, y4]))
        y_max = amax(array([y1, y2, y3, y4]))
        # Now make the plot square: enlarge the smaller of delta_x and delta_y
        # so they are equal
        delta_x = x_max - x_min
        delta_y = y_max - y_min
        if delta_x > delta_y:
            diff = 0.5*(delta_x - delta_y)
            y_min -= diff
            y_max += diff
        elif delta_y > delta_x:
            diff = 0.5*(delta_y - delta_x)
            x_min -= diff
            x_max += diff
        # Set up a grey square to fill the background with land
        x_reg, y_reg = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
        land_square = zeros(shape(x_reg))
        # Find bounds on variables in this region
        # Start with something impossible
        var_min = amax(data)
        var_max = amin(data)
        # Modify
        i = 0
        for elm in elements:
            if elm.cavity:
                if any(elm.x >= x_min) and any(elm.x <= x_max) and any(elm.y >= y_min) and any(elm.y <= y_max):
                    if data[i] < var_min:
                        var_min = data[i]
                    if data[i] > var_max:
                        var_max = data[i]
                i += 1
        # Set up colour map: centered on zero, but truncated
        bound = max(abs(var_min), abs(var_max))
        min_colour = (var_min+bound)/(2.0*bound)
        max_colour = (var_max+bound)/(2.0*bound)
        new_cmap = truncate_colormap(get_cmap('RdBu_r'), min_colour, max_colour)
        # Plot
        fig = figure(figsize=(15,12))
        ax = fig.add_subplot(1,1,1, aspect='equal')
        # Start with land background
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        # Add ice shelf elements
        img = PatchCollection(patches, cmap=new_cmap)
        img.set_array(array(data))
        img.set_edgecolor('face')
        img.set_clim(vmin=var_min, vmax=var_max)
        ax.add_collection(img)
        # Mask out the open ocean in white
        overlay = PatchCollection(mask_patches, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        axis('off')
        title('Change in melt rate (m/y)\n2091-2100 minus 2006-2015', fontsize=24)
        cbar = colorbar(img)
        cbar.ax.tick_params(labelsize=20)
        #fig.show()
        fig.savefig(fig_dir + fig_heads[index] + '.png')


# Truncate colourmap function from https://stackoverflow.com/questions/40929467/how-to-use-and-plot-only-a-part-of-a-colorbar-in-matplotlib
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=-1):
    if n== -1:
        n = cmap.N
    new_cmap = LinearSegmentedColormap.from_list('trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval), cmap(linspace(minval, maxval, n)))
    return new_cmap


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    file_path = raw_input("Path to output directory for RCP: ")
    fig_dir = raw_input("Path to directory to save figures into: ")
    cavity_melt_diff_rcp(mesh_path, file_path, fig_dir)
            
                    
            

    
