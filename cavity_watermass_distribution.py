from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from fesom_grid import *

def cavity_watermass_distribution():

    # Path to mesh directory
    mesh_path = '/short/y99/kaa561/FESOM/mesh/high_res/'
    # File containing temperature and salinity averaged over rep3 1992-2005
    ts_file = '/short/y99/kaa561/FESOM/highres_spinup/rep3/annual_avg.oce.mean.nc'
    # Number of temperature and salinity bins
    num_bins = 1000
    # Mesh parameters
    circumpolar = True
    cross_180 = False
    # Bounds on temperature and salinity bins (pre-computed, change if needed)
    min_salt = 32.8
    max_salt = 35
    min_temp = -3
    max_temp = 0.2

    # Read temperature and salinity at each 3D node
    id = Dataset(ts_file, 'r')
    temp = id.variables['temp'][0,:]
    salt = id.variables['salt'][0,:]
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

    # Make FESOM mesh elements
    elements = fesom_grid(mesh_path, circumpolar, cross_180)
    # Loop over elements
    for elm in elements:
        # Only consider ice shelf cavities
        if elm.cavity:
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
    fig = figure(figsize=(8,6))
    ax = fig.add_subplot(1,1,1)
    # Plot log of volume
    img = pcolor(salt_centres, temp_centres, log(ts_vals), cmap='jet')
    # Add surface freezing point line
    plot(salt_centres, freezing_pt, color='black', linestyle='dashed', linewidth=2)
    # Add dividing line at 34 psu
    tmp = -0.0575*34 + 1.7105e-3*sqrt(34**3) - 2.155e-4*34**2
    plot([34, 34], [tmp, 0.2], color='black', linestyle='dashed', linewidth=2)
    # Add dividing line at 34.5 psu
    tmp = -0.0575*34.5 + 1.7105e-3*sqrt(34.5**3) - 2.155e-4*34.5**2
    plot([34.5, 34.5], [tmp, -1], color='black', linestyle='dashed', linewidth=2)
    # Add dividing line at -1 C
    plot([34, 35], [-1, -1], color='black', linestyle='dashed', linewidth=2)
    # Label water masses
    text(33.25, -2.5, 'ISW', fontsize=20)
    text(33.5, -0.25, 'AASW', fontsize=20)
    text(34.12, -1.2, 'WW', fontsize=20)
    text(34.66, -1.2, 'HSSW', fontsize=20)
    text(34.5, -0.5, 'MCDW', fontsize=20)
    # Configure plot
    xlim([33, max_salt])
    ylim([min_temp, max_temp])
    xlabel('Salinity (psu)', fontsize=14)
    ylabel(r'Temperature ($^{\circ}$C)', fontsize=14)
    title('T-S distribution in ice shelf cavities, 1992-2005', fontsize=18)
    colorbar(img)
    # Label colourbar units
    text(35.5, -1, 'log of volume', fontsize=16, rotation=90)
    fig.show()    
    fig.savefig('watermass_key.png')


# Command-line interface
if __name__ == "__main__":

    cavity_watermass_distribution()
    
