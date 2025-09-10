from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from fesom_grid import *
from fesom_sidegrid import *

# Create one 2x2 plot for each major ice shelf: zonal slices of temperature
# (top) and salinity (bottom), comparing the low-res mesh (left) and high-res 
# (right). Longitudes to slice through, and latitude bounds, are pre-determined.
def zonal_cavity_ts_res ():

    # Paths to mesh directories
    mesh_path_low = '../FESOM/mesh/low_res/'
    mesh_path_high = '../FESOM/mesh/high_res/'
    # Paths to output files
    output_path_low = '../FESOM/lowres_spinup/rep3/'
    output_path_high = '../FESOM/highres_spinup/rep3/'
    file_name = 'annual_avg.oce.mean.nc'

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

    print('Building FESOM mesh')
    elm2D_low = fesom_grid(mesh_path_low)
    elm2D_high = fesom_grid(mesh_path_high)
    print('Reading temperature and salinity data')
    id = Dataset(output_path_low + file_name, 'r')
    temp_nodes_low = id.variables['temp'][0,:]
    salt_nodes_low = id.variables['salt'][0,:]
    id.close()
    id = Dataset(output_path_high + file_name, 'r')
    temp_nodes_high = id.variables['temp'][0,:]
    salt_nodes_high = id.variables['salt'][0,:]
    id.close()

    # Loop over ice shelves
    for index in range(num_shelves):
        print('Processing ' + shelf_names[index])
        # Figure out what to write on the title about longitude
        if lon0[index] < 0:
            lon_string = ' ('+str(-lon0[index])+r'$^{\circ}$W)'
        else:
            lon_string = ' ('+str(lon0[index])+r'$^{\circ}$E)'
        # Build arrays of SideElements making up zonal slices
        selements_temp_low = fesom_sidegrid(elm2D_low, temp_nodes_low, lon0[index], lat_max[index])
        selements_salt_low = fesom_sidegrid(elm2D_low, salt_nodes_low, lon0[index], lat_max[index])
        selements_temp_high = fesom_sidegrid(elm2D_high, temp_nodes_high, lon0[index], lat_max[index])
        selements_salt_high = fesom_sidegrid(elm2D_high, salt_nodes_high, lon0[index], lat_max[index])
        # Build array of quadrilateral patches for the plots, and data values
        # corresponding to each SideElement
        patches_low = []
        temp_low = []
        for selm in selements_temp_low:
            # Make patch
            coord = transpose(vstack((selm.y, selm.z)))
            patches_low.append(Polygon(coord, True, linewidth=0.))
            # Save data value
            temp_low.append(selm.var)
        temp_low = array(temp_low)
        # Salinity has same patches but different values
        salt_low = []
        for selm in selements_salt_low:
            salt_low.append(selm.var)
        salt_low = array(salt_low)
        # Repeat for high-res
        patches_high = []
        temp_high = []
        for selm in selements_temp_high:
            coord = transpose(vstack((selm.y, selm.z)))
            patches_high.append(Polygon(coord, True, linewidth=0.))
            temp_high.append(selm.var)
        temp_high = array(temp_high)
        salt_high = []
        for selm in selements_salt_high:
            salt_high.append(selm.var)
        salt_high = array(salt_high)
        # Find bounds on each variable
        temp_min = min(amin(temp_low), amin(temp_high))
        temp_max = max(amax(temp_low), amax(temp_high))
        salt_min = min(amin(salt_low), amin(salt_high))
        salt_max = max(amax(salt_low), amax(salt_high))
        # Find deepest depth
        # Start with 0
        depth_min = 0
        # Modify with low-res patches
        for selm in selements_temp_low:
            depth_min = min(depth_min, amin(selm.z))
        # Modify with high-res patches
        for selm in selements_temp_high:
            depth_min = min(depth_min, amin(selm.z))
        # Round down to nearest 50 metres
        depth_min = floor(depth_min/50)*50
        # Plot
        fig = figure(figsize=(18,12))
        # Low-res temperature
        ax = fig.add_subplot(2, 2, 1)
        img1 = PatchCollection(patches_low, cmap='jet')
        img1.set_array(temp_low)
        img1.set_edgecolor('face')
        img1.set_clim(vmin=temp_min, vmax=temp_max)
        ax.add_collection(img1)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        title(r'Low-res temperature ($^{\circ}$C)', fontsize=20)
        ylabel('Depth (m)', fontsize=16)
        # High-res temperature
        ax = fig.add_subplot(2, 2, 2)
        img2 = PatchCollection(patches_high, cmap='jet')
        img2.set_array(temp_high)
        img2.set_edgecolor('face')
        img2.set_clim(vmin=temp_min, vmax=temp_max)
        ax.add_collection(img2)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        title(r'High-res temperature ($^{\circ}$C)', fontsize=20)
        # Add colorbar for temperature
        cbaxes_temp = fig.add_axes([0.92, 0.575, 0.01, 0.3])
        cbar_temp = colorbar(img2, cax=cbaxes_temp)
        cbar_temp.ax.tick_params(labelsize=16)
        # Low-res salinity
        ax = fig.add_subplot(2, 2, 3)
        img3 = PatchCollection(patches_low, cmap='jet')
        img3.set_array(salt_low)
        img3.set_edgecolor('face')
        img3.set_clim(vmin=salt_min, vmax=salt_max)
        ax.add_collection(img3)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        title('Low-res salinity (psu)', fontsize=20)
        xlabel('Latitude', fontsize=16)
        ylabel('Depth (m)', fontsize=16)
        # High-res salinity
        ax = fig.add_subplot(2, 2, 4)
        img4 = PatchCollection(patches_high, cmap='jet')
        img4.set_array(salt_high)
        img4.set_edgecolor('face')
        img4.set_clim(vmin=salt_min, vmax=salt_max)
        ax.add_collection(img4)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        title('High-res salinity (psu)', fontsize=20)
        xlabel('Latitude', fontsize=16)
        # Add colorbar for salinity
        cbaxes_salt = fig.add_axes([0.92, 0.125, 0.01, 0.3])
        cbar_salt = colorbar(img4, cax=cbaxes_salt)
        cbar_salt.ax.tick_params(labelsize=16)
        # Main title
        suptitle(shelf_names[index] + lon_string, fontsize=28)
        #fig.show()
        fig.savefig(fig_heads[index] + '_zonal_ts.png')
        

# Command-line interface
if __name__ == "__main__":

    zonal_cavity_ts_res()
