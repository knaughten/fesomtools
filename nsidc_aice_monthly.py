from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from patches import *
from monthly_avg import *

# Plot monthly averages of sea ice concentration over the last year of
# simulation with FESOM, compared with NSIDC satellite data for 1995.
# Input:
# elements, patches = arrays of grid Elements and plotting patches created using
#                     patches.py
# file_path = path to FESOM output file containing one year of 5-day averages
#             for sea ice variables
# month = index for month to plot (0-based)
# save = optional boolean indicating to save the figure to a file, rather than
#        display it on the screen
# fig_name = if save=True, filename for figure
def nsidc_aice_monthly (elements, patches, file_path, month, save=False, fig_name=None):

    # Month names for plot titles
    month_name = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
    # NSIDC file paths
    nsidc_head = '/short/m68/kaa561/nsidc_aice/seaice_conc_monthly_sh'
    nsidc_head_0 = nsidc_head + '_f11_'
    nsidc_head_1 = nsidc_head + '_f13_'
    nsidc_tail = '_v02r00.nc'
    # Degrees to radians conversion factor
    deg2rad = pi/180.0

    # Get monthly average of the FESOM output
    fesom_data = monthly_avg(file_path, 'area', month)

    # Build an array of FESOM data values corresponding to each Element
    values = []
    for elm in elements:
        # For each element not in an ice shelf cavity, append the mean value
        # for the 3 component Nodes
        if not elm.cavity:
            values.append(mean([fesom_data[elm.nodes[0].id], fesom_data[elm.nodes[1].id], fesom_data[elm.nodes[2].id]]))

    # Construct NSIDC file path
    if month+1 < 10:
        nsidc_file = nsidc_head_0 + '19950' + str(month+1) + nsidc_tail
    else:
        nsidc_file = nsidc_head_1 + '1995' + str(month+1) + nsidc_tail

    # Read NSIDC grid and monthly data
    id = Dataset(nsidc_file, 'r')
    nsidc_lon = id.variables['longitude'][:,:]
    nsidc_lat = id.variables['latitude'][:,:]
    nsidc_data_tmp = id.variables['seaice_conc_monthly_cdr'][0,:,:]
    # Read std just for the land mask
    nsidc_mask = id.variables['stdev_of_seaice_conc_monthly_cdr'][0,:,:]
    id.close()

    # Set land mask on NSIDC sea ice concentration
    nsidc_data = ma.empty(shape(nsidc_data_tmp))
    nsidc_data[:,:] = 0.0
    nsidc_data[~nsidc_mask.mask] = nsidc_data_tmp[~nsidc_mask.mask]
    nsidc_data[nsidc_mask.mask] = ma.masked

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
    # NSIDC
    ax = fig.add_subplot(1,2,1, aspect='equal')
    contourf(nsidc_x, nsidc_y, nsidc_data, lev)
    title('NSIDC', fontsize=24)
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    # FESOM
    ax = fig.add_subplot(1,2,2, aspect='equal')
    img = PatchCollection(patches, cmap=jet)
    img.set_array(array(values))
    img.set_clim(vmin=0, vmax=1)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    title('FESOM', fontsize=24)
    # Add a horizontal colorbar at the bottom
    cbaxes = fig.add_axes([0.35, 0.04, 0.3, 0.04])
    cbar = colorbar(img, orientation='horizontal', ticks=arange(0,1+0.25,0.25), cax=cbaxes)
    cbar.ax.tick_params(labelsize=20)
    # Add the main title
    suptitle(month_name[month] + ' sea ice concentration', fontsize=30)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    mesh_path = input("Path to FESOM mesh directory: ")
    file_path = input("Path to output ice.mean.nc containing one year of 5-day averages: ")
    month = int(input("Month number (1-12): ")) - 1
    action = input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None

    # Pre-compute the FESOM mesh (so that it doesn't have to be recomputed
    # if the user wants to repeat)
    circumpolar = True
    mask_cavities = True
    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)
    # Make the plot
    nsidc_aice_monthly(elements, patches, file_path, month, save, fig_name)

    while True:
        # Repeat until the user wants to exit
        repeat = input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to parameters until the user is done
                changes = input("Enter a parameter to change: (1) file path, (2) month number, (3) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New FESOM file
                        file_path = input("Path to one year of 5-day averages for sea ice variables: ")
                    elif int(changes) == 2:
                        # New month
                        month = int(input("Month number (1-12): ")) - 1
                    elif int(changes) == 3:
                        # Change from save to display, or vice versa
                        save = not save
            if save:
                # Get a new figure name
                fig_name = input("File name for figure: ")
            # Make the plot
            nsidc_aice_monthly(elements, patches, file_path, month, save, fig_name)
        else:
            break
                
    
