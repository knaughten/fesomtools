from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from fesom_grid import *
from unesco import *

# Calculate the meridional overturning streamfunction in latitude-density space.
# Input:
# mesh_path = path to FESOM mesh directory
# file_path = path to FESOM output oce.mean.nc file, containing 1 year of output
# save = optional boolean indicating to save the figure, rather than display
# fig_name = if save=True, filename for figure
def moc_lat_density (mesh_path, file_path, save=False, fig_name=None):

    # Options for grid objects
    circumpolar = False
    cross_180 = False

    # Read vertical velocity, temperature, and salinity at every node
    id = Dataset(file_path, 'r')
    w = mean(id.variables['w'][:,:], axis=0)
    temp = mean(id.variables['temp'][:,:], axis=0)
    salt = mean(id.variables['salt'][:,:], axis=0)
    id.close()

    # Calculate potential density (depth 0) at every node
    density = unesco(temp, salt, zeros(shape(temp)))-1000

    # Build FESOM grid
    elements = fesom_grid(mesh_path, circumpolar, cross_180)

    # Set up arrays of vertical transport, latitude, upstream density, and
    # downstream density at every interface between vertical layers of elements
    transport_all = []
    lat_all = []
    density_us_all = []
    density_ds_all = []
    # Loop over 2D elements
    for elm in elements:
        # Get area and latitude (average over 3 nodes)
        area = elm.area()
        lat = mean(elm.lat)
        nodes_above = [elm.nodes[0], elm.nodes[1], elm.nodes[2]]
        nodes = [nodes_above[0].below, nodes_above[1].below, nodes_above[2].below]
        # Loop from the second layer from the surface, down to the second layer
        # from the bottom
        while True:
            nodes_below = [nodes[0].below, nodes[1].below, nodes[2].below]
            if None in nodes_below:
                # Reached the bottom
                break
            # Vertical velocity average over 3 nodes
            w_avg = mean([w[nodes[0].id], w[nodes[1].id], w[nodes[2].id]])
            # Vertical transport through this triangular interface
            transport = abs(w_avg)*area*1e-6
            # Density average over 3D triangular prism above
            density_above = mean([density[nodes[0].id], density[nodes[1].id], density[nodes[2].id], density[nodes_above[0].id], density[nodes_above[1].id], density[nodes_above[2].id]])
            # Density average over 3D triangular prism below
            density_below = mean([density[nodes[0].id], density[nodes[1].id], density[nodes[2].id], density[nodes_below[0].id], density[nodes_below[1].id], density[nodes_below[2].id]])
            # Figure out which is triangular prism upstream and which is
            # downstream; save the density values correspondingly
            if w_avg > 0:
                density_us = density_below
                density_ds = density_above
            else:
                density_us = density_above
                density_ds = density_below
            # Save vertical transport, latitude, upstream and downstream
            # densities for this interface
            transport_all.append(transport)
            lat_all.append(lat)
            density_us_all.append(density_us)
            density_ds_all.append(density_ds)
            # Get ready for next layer down
            nodes_above = nodes
            nodes = nodes_below

    # Get regular values of latitude and density
    lat_reg = linspace(-90, 90, num=50)
    density_reg = linspace(floor(amin(density)), ceil(amax(density)), num=25)
    # Set up array for overturning streamfunction
    moc = zeros([size(density_reg), size(lat_reg)])
    # Loop over latitude
    for j in range(size(lat_reg)):
        print 'Processing latitude ' + str(j+1) + ' of ' + str(size(lat_reg))
        # Make a flag which is 1 for interfaces south of the current latitude,
        # 0 otherwise
        flag_lat = zeros(shape(lat_all))
        index = lat_all <= lat_reg[j]
        flag_lat[index] = 1
        # Loop over density
        for k in range(size(density_reg)):
            # Make a flag which is 1 or -1 (depending on direction) for
            # interfaces where the upstream-downstream density gradient crosses
            # the current density, 0 otherwise
            flag_density = zeros(shape(density_us_all))
            index = (density_us_all <= density_reg[k])*(density_ds_all >= density_reg[k])
            flag_density[index] = 1
            index = (density_ds_all <= density_reg[k])*(density_us_all >= density_reg[k])
            flag_density[index] = -1
            # Calculate MOC
            moc[k,j] = sum(transport_all*flag_lat*flag_density)            

    # Make colour levels
    bound = amax(abs(moc))
    lev = linspace(-bound, bound, num=50)

    # Plot
    fig = figure()
    img = contourf(lat_reg, density_reg, moc, lev, cmap='RdBu_r')
    ylim([density_reg[-1], density_reg[0]])
    xlabel('Latitude')
    ylabel(r'Density (kg/m$^3$)')
    title('Meridional Overturning Streamfunction (Sv)')
    colorbar(img)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()
    

# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    file_path = raw_input("Path to FESOM oce.mean.nc file containing one year of data: ")
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("Filename for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    moc_lat_density(mesh_path, file_path, save, fig_name)
