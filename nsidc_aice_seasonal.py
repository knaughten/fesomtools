from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from patches import *

# Plot seasonal averages of sea ice concentration over the last year of
# simulation with FESOM, compared with NSIDC satellite data for 1995.
# Input:
# mesh_path = path to FESOM mesh directory
# file_path1 = path to a FESOM output file containing one year of 5-day
#              averages for sea ice variables (we will just use December)
# file_path2 = path to a FESOM output file containing the following year of
#              5-day averages for sea ice variables (we will use January
#              through November)
# save = optional boolean indicating to save the figure to a file, rather than
#        display it on the screen
# fig_name = if save=True, filename for figure
def nsidc_aice_seasonal (mesh_path, file_path1, file_path2, save=False, fig_name=None):

    # FESOM parameters
    circumpolar=True
    mask_cavities=True
    # Season names for plot titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    # NSIDC file paths
    nsidc_head = '/short/m68/kaa561/nsidc_aice/seaice_conc_monthly_sh'
    nsidc_head_0 = nsidc_head + '_f11_'
    nsidc_head_1 = nsidc_head + '_f13_'
    nsidc_tail = '_v02r00.nc'
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Number of days per month (just for CICE)
    ndays_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    # Build FESOM mesh
    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)

    # Get seasonal averages of the FESOM output
    # This is hard-coded and ugly
    id = Dataset(file_path1, 'r')
    n2d = id.variables['area'].shape[1]
    fesom_data = zeros([4, n2d])
    # DJF: 1/5 of index 67 (1-based) and indices 68-73 in file1; indices 1-11
    # and 4/5 of index 12 in file2; 90 days in total
    fesom_data[0,:] = id.variables['area'][66,:] + sum(id.variables['area'][67:73,:]*5, axis=0)
    id.close()
    id = Dataset(file_path2, 'r')
    fesom_data[0,:] += sum(id.variables['area'][0:11,:]*5, axis=0) + id.variables['area'][11,:]*4
    fesom_data[0,:] /= 90
    # MAM: 1/5 of index 12, indices 13-30, and 1/5 of index 31 in file2;
    # 92 days in total
    fesom_data[1,:] = id.variables['area'][11,:] + sum(id.variables['area'][12:30,:]*5, axis=0) + id.variables['area'][30,:]
    fesom_data[1,:] /= 92
    # JJA: 4/5 of index 31, indices 32-48, and 3/5 of index 49 in file2;
    # 92 days in total
    fesom_data[2,:] = id.variables['area'][30,:]*4 + sum(id.variables['area'][31:48]*5, axis=0) + id.variables['area'][48,:]*3
    fesom_data[2,:] /= 92
    # SON: 2/5 of index 49, indices 50-66, and 4/5 of index 67 in file2;
    # 91 days in total
    fesom_data[3,:] = id.variables['area'][48,:]*2 + sum(id.variables['area'][49:66,:]*5, axis=0) + id.variables['area'][66,:]*4
    fesom_data[3,:] /= 91
    id.close()

    # Read NSIDC grid from the January file
    id = Dataset(nsidc_head_0 + '199501' + nsidc_tail, 'r')
    nsidc_lon = id.variables['longitude'][:,:]
    nsidc_lat = id.variables['latitude'][:,:]
    id.close()

    # Initialise seasonal averages of NSIDC data
    nsidc_data = ma.empty([4, size(nsidc_lon,0), size(nsidc_lon,1)])
    nsidc_data[:,:] = 0.0
    # Process one season at a time
    for season in range(4):
        # Figure out which months we care about
        if season == 0:
            # DJF
            nsidc_months = [12, 1, 2]
        elif season == 1:
            # MAM
            nsidc_months = [3, 4, 5]
        elif season == 2:
            # JJA
            nsidc_months = [6, 7, 8]
        elif season == 3:
            # SON
            nsidc_months = [9, 10, 11]
        season_days = 0 # Number of days in season; this will be incremented

        # Process one month at a time
        for month in nsidc_months:
            # Construct NSIDC file path
            if month < 10:
                nsidc_file = nsidc_head_0 + '19950' + str(month) + nsidc_tail
            else:
                nsidc_file = nsidc_head_1 + '1995' + str(month) + nsidc_tail
            # Read concentration data
            id = Dataset(nsidc_file, 'r')
            nsidc_data_raw = id.variables['seaice_conc_monthly_cdr'][0,:,:]
            # Read std just for the mask
            nsidc_mask = id.variables['stdev_of_seaice_conc_monthly_cdr'][0,:,:]
            id.close()
            # Set land mask
            nsidc_data_tmp = ma.empty(shape(nsidc_data_raw))
            nsidc_data_tmp[:,:] = 0.0
            nsidc_data_tmp[~nsidc_mask.mask] = nsidc_data_raw[~nsidc_mask.mask]
            nsidc_data_tmp[nsidc_mask.mask] = ma.masked
            # Accumulate master array, weighted with number of days per month
            nsidc_data[season,:,:] += nsidc_data_tmp*ndays_month[month-1]
            season_days += ndays_month[month-1]

        # Convert from sum to average
        nsidc_data[season,:,:] /= season_days

    # Convert to spherical coordinates
    nsidc_x = -(nsidc_lat+90)*cos(nsidc_lon*deg2rad+pi/2)
    nsidc_y = (nsidc_lat+90)*sin(nsidc_lon*deg2rad+pi/2)

    # Find boundaries for each side of plot based on extent of NSIDC grid
    bdry1 = amax(nsidc_x[:,0])
    bdry2 = amin(nsidc_x[:,-1])
    bdry3 = amin(nsidc_y[:,0])
    bdry4 = amax(nsidc_y[:,-1])

    # Set consistent colour levels
    lev = linspace(0, 1, num=50)

    # Plot
    fig = figure(figsize=(20,9))
    # Loop over seasons
    for season in range(4):
        # NSIDC
        ax = fig.add_subplot(2, 4, season+1, aspect='equal')
        contourf(nsidc_x, nsidc_y, nsidc_data[season,:,:], lev)
        if season == 0:
            text(-39, 0, 'NSIDC', fontsize=24, ha='right')
        title(season_names[season], fontsize=24)
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
        # Build an array of FESOM data values corresponding to each Element
        values = []
        for elm in elements:
            # For each element not in an ice shelf cavity, append the mean
            # value for the 3 component Nodes
            if not elm.cavity:
                values.append(mean([fesom_data[season,elm.nodes[0].id], fesom_data[season,elm.nodes[1].id], fesom_data[season,elm.nodes[2].id]]))
        # Plot FESOM data
        ax = fig.add_subplot(2, 4, season+5, aspect='equal')
        img = PatchCollection(patches, cmap=jet)
        img.set_array(array(values))
        img.set_clim(vmin=0, vmax=1)
        img.set_edgecolor('face')
        ax.add_collection(img)
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
        if season == 0:
            text(-39, 0, 'FESOM', fontsize=24, ha='right')
    # Add a horizontal colorbar at the bottom
    cbaxes = fig.add_axes([0.25, 0.04, 0.5, 0.02])
    cbar = colorbar(img, orientation='horizontal', ticks=arange(0,1+0.25,0.25), cax=cbaxes)
    cbar.ax.tick_params(labelsize=16)
    # Add the main title
    suptitle('Sea ice concentration', fontsize=30)
    # Decrease space between plots
    subplots_adjust(wspace=0.025,hspace=0.025)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()
        

# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    file_path1 = raw_input("Path to output ice.mean.nc containing one year of 5-day averages (December will be used): ")
    file_path2 = raw_input("Path to the following ice.mean.nc containing 5-day averages for the next year (January through November will be used): ")
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None

    nsidc_aice_seasonal(mesh_path, file_path1, file_path2, save, fig_name)

    

    

    
