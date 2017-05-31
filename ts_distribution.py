from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from fesom_grid import *
from unesco import *

# Make a temperature-salinity distribution plot, showing the volume of each
# water mass with a log scale, including the surface freezing point line and
# density contours.
# Input:
# mesh_path = path to FESOM mesh directory
# file_path = path to FESOM oce.mean.nc file
# tstep = time index in file_path to plot (1-based)
# save = optional boolean indicating to save the figure rather than display
# fig_name = if save=True, filename for figure
def ts_distribution (mesh_path, file_path, tstep, save=False, fig_name=None):

    # Northern boundary of water masses to consider
    nbdry = -50
    # Number of temperature and salinity bins
    num_bins = 1000
    # Plotting parameters
    circumpolar = False
    cross_180 = False
    # Bounds on temperature and salinity bins (pre-computed, change if needed)
    min_salt = 31.8
    max_salt = 35.2
    min_temp = -3
    max_temp = 12

    # Read temperature and salinity at each 3D node
    id = Dataset(file_path, 'r')
    temp = id.variables['temp'][tstep-1,:]
    salt = id.variables['salt'][tstep-1,:]
    id.close()

    # Calculate boundaries of temperature bins
    temp_bins = linspace(min_temp, max_temp, num=num_bins)
    # Calculate centres of temperature bins (for plotting)
    temp_centres = 0.5*(temp_bins[:-1] + temp_bins[1:])
    # Repeat for salinity
    salt_bins = linspace(min_salt, max_salt, num=num_bins)
    salt_centres = 0.5*(salt_bins[:-1] + salt_bins[1:])
    # Set up a 2D array of temperature bins x salinity bins to increment with
    # volume of water masses
    ts_vals = zeros([size(temp_centres), size(salt_centres)])

    # Calculate surface freezing point as a function of salinity: this is the
    # equation the FESOM sea ice code uses
    freezing_pt = -0.0575*salt_centres + 1.7105e-3*sqrt(salt_centres**3) - 2.155e-4*salt_centres**2
    # Get 2D versions of the temperature and salinity bins
    salt_2d, temp_2d = meshgrid(salt_centres, temp_centres)
    # Calculate potential density of each combination of temperature and
    # salinity bins
    density = unesco(temp_2d, salt_2d, zeros(shape(temp_centres)))-1000
    # Density contours to plot
    density_lev = arange(24.4, 28.4, 0.2)

    # Make FESOM grid elements
    elements = fesom_grid(mesh_path, circumpolar, cross_180)

    # Loop over elements
    for elm in elements:
        # See if we're in the region of interest
        if all(elm.lat < nbdry):
            # Get area of 2D triangle
            area = elm.area()
            nodes = [elm.nodes[0], elm.nodes[1], elm.nodes[2]]
            # Loop downward
            while True:
                if nodes[0].below is None or nodes[1].below is None or nodes[2].below is None:
                    # We've reached the bottom
                    break
                # Calculate average temperature, salinity, and layer thickness
                # over this 3D triangular prism
                temp_vals = []
                salt_vals = []
                dz = []
                for i in range(3):
                    # Average temperature over 6 nodes
                    temp_vals.append(temp[nodes[i].id])
                    temp_vals.append(temp[nodes[i].below.id])
                    # Average salinity over 6 nodes
                    salt_vals.append(salt[nodes[i].id])
                    salt_vals.append(salt[nodes[i].below.id])
                    # Average dz over 3 vertical edges
                    dz.append(abs(nodes[i].depth - nodes[i].below.depth))
                    # Get ready for next repetition of loop
                    nodes[i] = nodes[i].below
                temp_elm = mean(array(temp_vals))
                salt_elm = mean(array(salt_vals))
                # Calculate volume of 3D triangular prism
                volume = area*mean(array(dz))
                # Figure out which bins this falls into
                temp_index = nonzero(temp_bins > temp_elm)[0][0] - 1
                salt_index = nonzero(salt_bins > salt_elm)[0][0] - 1
                # Increment bins with volume
                ts_vals[temp_index, salt_index] += volume

    # Mask bins with zero volume
    ts_vals = ma.masked_where(ts_vals==0, ts_vals)

    # Plot
    fig = figure(figsize=(12,12))
    # Log scale is more visible
    img=pcolor(salt_centres, temp_centres, log(ts_vals), cmap='jet')
    # Add surface freezing point line
    plot(salt_centres, freezing_pt, color='black', linestyle='dashed')
    # Add density contours
    cs=contour(salt_centres, temp_centres, density, density_lev, colors=(0.6,0.6,0.6), linestyles='dotted')
    # Label density contours
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


# Command-line interface
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
    
