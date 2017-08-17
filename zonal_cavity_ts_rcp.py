from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from fesom_grid import *
from fesom_sidegrid import *

# Create one 3x2 plot for each major ice shelf: zonal slices of temperature
# (top) and salinity (bottom), comparing the first 10 years of the RCP (left),
# the last 10 years (middle), and the difference (right). Longitudes to slice
# through, and latitude bounds, are pre-determined.
# Input:
# mesh_path = path to FESOM mesh directory
# output_path = path to directory containing the files
#               annual_avg.oce.mean.2006.2015.nc and
#               annual_avg.oce.mean.2091.2100.nc, which contain 3D temperature
#               and salinity averaged over 2006-2015 and 2091-2100 respectively
# fig_dir = optional string containing path to directory to store figures in
def zonal_cavity_ts_rcp (mesh_path, output_path, fig_dir=''):

    file_name_beg = 'annual_avg.oce.mean.2006.2015.nc'
    file_name_end = 'annual_avg.oce.mean.2091.2100.nc'

    # Name of each ice shelf
    shelf_names = ['Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    # Beginnings of filenames for figures
    fig_heads = ['larsen_d', 'larsen_c', 'wilkins_georgevi_stange', 'ronne_filchner', 'abbot', 'pig', 'thwaites', 'dotson', 'getz', 'nickerson', 'sulzberger', 'mertz', 'totten_moscowuni', 'shackleton', 'west', 'amery', 'prince_harald', 'baudouin_borchgrevink', 'lazarev', 'nivl', 'fimbul_jelbart_ekstrom', 'brunt_riiser_larsen', 'ross']
    # Longitudes intersecting each ice shelf
    lon0 = [-60, -62, -68, -55, -93, -101, -106, -113, -120, -145, -150, 145, 116, 96, 85, 71, 36, 25, 15, 11, -1, -20, 180]
    # Latitude bounds for each ice shelf
    lat_min = [-73.1, -69.35, -73.1, -82.6, -73.28, -75.4, -75.5, -75, -74.9, -75.9, -77.8, -67.7, -67.17, -66.67, -67.25, -72, -69.7, -71, -70.4, -70.75, -71.83, -75.6, -84.6]
    lat_max = [-72, -66.13, -70, -75.5, -72.3, -74.4, -74.67, -74, -73.5, -75.3, -76.41, -67, -66.5, -64.83, -66.25, -68.5, -68.7, -69.9, -69.33, -69.83, -69.33, -72.9, -77]
    num_shelves = len(shelf_names)

    print 'Building FESOM mesh'
    elm2D = fesom_grid(mesh_path)
    print 'Reading temperature and salinity data'
    id = Dataset(output_path + file_name_beg, 'r')
    temp_nodes_beg = id.variables['temp'][0,:]
    salt_nodes_beg = id.variables['salt'][0,:]
    id.close()
    id = Dataset(output_path + file_name_end, 'r')
    temp_nodes_end = id.variables['temp'][0,:]
    salt_nodes_end = id.variables['salt'][0,:]
    id.close()
    temp_nodes_diff = temp_nodes_end - temp_nodes_beg
    salt_nodes_diff = salt_nodes_end - salt_nodes_beg

    # Loop over ice shelves
    for index in range(num_shelves):
        print 'Processing ' + shelf_names[index]
        # Figure out what to write on the title about longitude
        if lon0[index] < 0:
            lon_string = ' ('+str(-lon0[index])+r'$^{\circ}$W)'
        else:
            lon_string = ' ('+str(lon0[index])+r'$^{\circ}$E)'
        # Build arrays of SideElements making up zonal slices
        selements_temp_beg = fesom_sidegrid(elm2D, temp_nodes_beg, lon0[index], lat_max[index])
        selements_salt_beg = fesom_sidegrid(elm2D, salt_nodes_beg, lon0[index], lat_max[index])
        selements_temp_end = fesom_sidegrid(elm2D, temp_nodes_end, lon0[index], lat_max[index])
        selements_salt_end = fesom_sidegrid(elm2D, salt_nodes_end, lon0[index], lat_max[index])
        selements_temp_diff = fesom_sidegrid(elm2D, temp_nodes_diff, lon0[index], lat_max[index])
        selements_salt_diff = fesom_sidegrid(elm2D, salt_nodes_diff, lon0[index], lat_max[index])
        # Build array of quadrilateral patches for the plots, and data values
        # corresponding to each SideElement
        patches = []
        temp_beg = []
        for selm in selements_temp_beg:
            # Make patch
            coord = transpose(vstack((selm.y, selm.z)))
            patches.append(Polygon(coord, True, linewidth=0.))
            # Save data value
            temp_beg.append(selm.var)
        temp_beg = array(temp_beg)
        # Other variables have same patches but different values
        salt_beg = []
        for selm in selements_salt_beg:
            salt_beg.append(selm.var)
        salt_beg = array(salt_beg)
        temp_end = []
        for selm in selements_temp_end:
            temp_end.append(selm.var)
        temp_end = array(temp_end)
        salt_end = []
        for selm in selements_salt_end:
            salt_end.append(selm.var)
        salt_end = array(salt_end)
        temp_diff = []
        for selm in selements_temp_diff:
            temp_diff.append(selm.var)
        temp_diff = array(temp_diff)
        salt_diff = []
        for selm in selements_salt_diff:
            salt_diff.append(selm.var)
        salt_diff = array(salt_diff)
        # Find bounds on each variable
        temp_min = min(amin(temp_beg), amin(temp_end))
        temp_max = max(amax(temp_beg), amax(temp_end))
        temp_max_diff = amax(abs(temp_diff))
        salt_min = min(amin(salt_beg), amin(salt_end))
        salt_max = max(amax(salt_beg), amax(salt_end))
        salt_max_diff = amax(abs(salt_diff))
        # Find deepest depth
        depth_min = 0
        for selm in selements_temp_beg:
            depth_min = min(depth_min, amin(selm.z))
        # Round down to nearest 50 metres
        depth_min = floor(depth_min/50)*50
        # Plot
        fig = figure(figsize=(24,12))
        # Temperature (beginning)
        ax = fig.add_subplot(2, 3, 1)
        img = PatchCollection(patches, cmap='jet')
        img.set_array(temp_beg)
        img.set_edgecolor('face')
        img.set_clim(vmin=temp_min, vmax=temp_max)
        ax.add_collection(img)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        title(r'Temperature ($^{\circ}$C), 2006-2015', fontsize=20)
        ylabel('Depth (m)', fontsize=16)
        # Add colorbar for absolute temperature
        cbaxes_temp = fig.add_axes([0.05, 0.575, 0.01, 0.3])
        cbar_temp = colorbar(img, cax=cbaxes_temp)
        cbar_temp.ax.tick_params(labelsize=16)
        # Temperature (end)
        ax = fig.add_subplot(2, 3, 2)
        img = PatchCollection(patches, cmap='jet')
        img.set_array(temp_end)
        img.set_edgecolor('face')
        img.set_clim(vmin=temp_min, vmax=temp_max)
        ax.add_collection(img)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        title(r'Temperature ($^{\circ}$C), 2091-2100', fontsize=20)
        ylabel('Depth (m)', fontsize=16)
        # Temperature (difference)
        ax = fig.add_subplot(2, 3, 3)
        img = PatchCollection(patches, cmap='RdBu_r')
        img.set_array(temp_diff)
        img.set_edgecolor('face')
        img.set_clim(vmin=-temp_max_diff, vmax=temp_max_diff)
        ax.add_collection(img)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        title(r'Temperature ($^{\circ}$C), change', fontsize=20)
        ylabel('Depth (m)', fontsize=16)
        # Add colorbar for temperature difference
        cbaxes_dtemp = fig.add_axes([0.92, 0.575, 0.01, 0.3])
        cbar_dtemp = colorbar(img, cax=cbaxes_dtemp)
        cbar_dtemp.ax.tick_params(labelsize=16)
        # Salinity (beginning)
        ax = fig.add_subplot(2, 3, 4)
        img = PatchCollection(patches, cmap='jet')
        img.set_array(salt_beg)
        img.set_edgecolor('face')
        img.set_clim(vmin=salt_min, vmax=salt_max)
        ax.add_collection(img)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        title('Salinity (psu), 2006-2015', fontsize=20)
        ylabel('Depth (m)', fontsize=16)
        # Add colorbar for absolute salinity
        cbaxes_salt = fig.add_axes([0.05, 0.125, 0.01, 0.3])
        cbar_salt = colorbar(img, cax=cbaxes_salt)
        cbar_salt.ax.tick_params(labelsize=16)
        # Salinity (end)
        ax = fig.add_subplot(2, 3, 5)
        img = PatchCollection(patches, cmap='jet')
        img.set_array(salt_end)
        img.set_edgecolor('face')
        img.set_clim(vmin=salt_min, vmax=salt_max)
        ax.add_collection(img)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        title('Salinity (psu), 2091-2100', fontsize=20)
        ylabel('Depth (m)', fontsize=16)
        # Salinity (difference)
        ax = fig.add_subplot(2, 3, 6)
        img = PatchCollection(patches, cmap='RdBu_r')
        img.set_array(salt_diff)
        img.set_edgecolor('face')
        img.set_clim(vmin=-salt_max_diff, vmax=salt_max_diff)
        ax.add_collection(img)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        title('Salinity (psu), change', fontsize=20)
        ylabel('Depth (m)', fontsize=16)
        # Add colorbar for salinity difference
        cbaxes_dsalt = fig.add_axes([0.92, 0.125, 0.01, 0.3])
        cbar_dsalt = colorbar(img, cax=cbaxes_dsalt)
        cbar_dsalt.ax.tick_params(labelsize=16)
        # Main title
        suptitle(shelf_names[index] + lon_string, fontsize=28)
        #fig.show()
        fig.savefig(fig_dir + fig_heads[index] + '_zonal_ts.png')
        

# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    output_path = raw_input("Path to output directory for RCP: ")
    fig_dir = raw_input("Directory to store figures (blank for current directory, otherwise ending in /): ")
    zonal_cavity_ts_rcp(mesh_path, output_path, fig_dir)
