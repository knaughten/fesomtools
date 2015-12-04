from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.patches import Polygon
from patches import *


# Create a plot of a specified variable on a lon-lat domain.
# Input:
# file_path = string containing path to FESOM output file
# var_name = string containing name of variable in file_path
# depth_key = int specifying whether to plot surface nodes (0), bottom nodes
#             (1), vertical average throughout the entire water column (2), 
#             a single specified depth (3), or vertical average between two
#             specified depths (4)
# depth = if depth_key==3, specified depth in m
# depth_bounds = if depth_key==4, array containing shallow and deep bounds for
#                vertical average
# tstep = int specifying index of time axis in file_path
# circumpolar = boolean flag indicating whether to use the circumpolar domain
#               or the global domain
# elements = array of Elements for the global grid (created using fesom_grid)
# patches = array of Polygon patches corresponding to elements (created using
#           make_patches)
# mask_cavities = optional boolean flag indicating whether to mask ice shelf
#                 cavities in grey
# save = optional boolean flag indicating whether to save plot to file
#        (otherwise will display on screen)
# fig_name = optional string containing name of figure file, if save = True
# set_limits = optional boolean flag indicating whether or not to set manual
#              limits on colourbar (otherwise limits determined automatically)
# limits = optional array containing min and max limits
def lonlat_plot (file_path, var_name, depth_key, depth, depth_bounds, tstep, circumpolar, elements, patches, mask_cavities=False, save=False, fig_name=None, set_limits=False, limits=None):

    # Set bounds for domain
    if circumpolar:
        # Northern boundary 60S
        lat_max = -60+90
    else:
        lon_min = -180
        lon_max = 180
        lat_min = -90
        lat_max = 90
        # Configure position of latitude and longitude labels
        lon_ticks = arange(-120, 120+1, 60)
        lat_ticks = arange(-60, 60+1, 30)
    # Set font sizes
    font_sizes = [18, 16, 12]
    # Seconds per year, for conversion of ice shelf melt rate
    sec_per_year = 365.25*24*3600

    # Read data
    file = Dataset(file_path, 'r')
    varid = file.variables[var_name]
    data = varid[tstep-1,:]
    # Set descriptive variable name and units for title
    if var_name == 'area':
        name = 'ice concentration'
        units = 'fraction'
    elif var_name == 'wnet':
        # Convert from m/s to m/y
        data = data*sec_per_year
        name = 'ice shelf melt rate'
        units = 'm/y'
    else:
        name = varid.getncattr('description')
        units = varid.getncattr('units')
    if depth_key == 0:
        if var_name in ['temp', 'salt', 'u', 'v']:
            depth_string = 'at surface'
        else:
            depth_string = ''
    elif depth_key == 1:
        depth_string = 'at bottom'
    elif depth_key == 2:
        depth_string = 'vertically averaged'
    elif depth_key == 3:
        depth_string = 'at ' + str(depth) + ' m'
    elif depth_key == 4:
        depth_string = 'vertically averaged between '+str(depth_bounds[0])+' and '+str(depth_bounds[1])+' m'

    # Build an array of data values corresponding to each Element
    values = []
    plot_patches = []
    for elm in elements:
        # If mask_cavities is true, only include elements which are not inside
        # an ice shelf cavity; otherwise, include all elements
        if (mask_cavities and not elm.cavity) or (not mask_cavities):

            if depth_key == 0:
                # Surface nodes; this is easy
                # Average the data value for each of the three component nodes
                values.append(mean([data[elm.nodes[0].id], data[elm.nodes[1].id], data[elm.nodes[2].id]]))

            elif depth_key == 1:
                # Bottom nodes
                values_tmp = []
                # For each of the three component nodes, find the id of the
                # bottom node beneath it
                for node in elm.nodes:
                    id = node.find_bottom().id
                    values_tmp.append(data[id])
                # Average over these three values
                values.append(mean(values_tmp))

            elif depth_key == 2:                
                # Vertical average throughout entire water column
                # First calculate volume: area of triangular face * water
                # column thickness (mean of three corners)
                wct = []
                for node in elm.nodes:
                    wct.append(abs(node.find_bottom().depth - node.depth))
                area = elm.area()
                volume = area*mean(array(wct))
                # Now integrate down
                integral = 0
                nodes = [elm.nodes[0], elm.nodes[1], elm.nodes[2]]
                while True:
                    if nodes[0].below is None or nodes[1].below is None or nodes[2].below is None:
                        break
                    # Calculate mean of data at six corners of this triangular
                    # prism, and mean depths at three edges
                    values_tmp = []
                    dz_tmp = []
                    # Must loop over indices of nodes array, not entries,
                    # because we want to change the contents of the nodes array
                    for i in range(3):
                        values_tmp.append(data[nodes[i].id])
                        values_tmp.append(data[nodes[i].below.id])
                        dz_tmp.append(abs(nodes[i].depth - nodes[i].below.depth))
                        # Get ready for next iteration of the while loop
                        nodes[i] = nodes[i].below
                    # Integrand is mean of data at corners * area of upper face
                    # * mean of depths at edges
                    integral += mean(array(values_tmp))*area*mean(array(dz_tmp))
                # All done; divide integral by volume to get the average
                values.append(integral/volume)    

            elif depth_key == 3:
                # Specified depth
                values_tmp = []
                # For each of the three component nodes, linearly interpolate
                # to the correct depth
                for i in range(3):
                    # Find the ids of the nodes above and below this depth,
                    # and the coefficients for the linear interpolation
                    id1, id2, coeff1, coeff2 = elm.nodes[i].find_depth(depth)
                    if any(isnan(array([id1, id2, coeff1, coeff2]))):
                        # No such node at this depth
                        values_tmp.append(NaN)
                    else:
                        values_tmp.append(coeff1*data[id1] + coeff2*data[id2])
                if any (isnan(array(values_tmp))):
                    pass
                else:
                    values.append(mean(values_tmp))
                    coord = transpose(vstack((elm.x, elm.y)))
                    # Make new patches for elements which exist at this depth
                    plot_patches.append(Polygon(coord, True, linewidth=0.))

            elif depth_key == 4:
                # Vertical average between two specified depths
                shallow_bound = depth_bounds[0]
                deep_bound = depth_bounds[1]
                area = elm.area()
                # Volume is area of triangular face * difference between
                # depth bounds
                volume = abs(deep_bound - shallow_bound)*area
                nodes = [elm.nodes[0], elm.nodes[1], elm.nodes[2]]
                # Check if we're already too deep
                if nodes[0].depth > deep_bound or nodes[1].depth > deep_bound or nodes[2].depth > deep_bound:
                    pass
                # Check if the seafloor is shallower than shallow_bound
                elif nodes[0].find_bottom().depth <= shallow_bound or nodes[1].find_bottom().depth <= shallow_bound or nodes[2].find_bottom().depth <= shallow_bound:
                    pass

                else:
                    # Here we go

                    # Find the first 3D element which is entirely below
                    # shallow_bound
                    while True:
                        if nodes[0].depth > shallow_bound and nodes[1].depth > shallow_bound and nodes[2].depth > shallow_bound:
                            # Save these nodes
                            first_nodes = []
                            for i in range(3):
                                first_nodes.append(nodes[i])
                            break
                        else:
                            for i in range(3):
                                nodes[i] = nodes[i].below

                    # Integrate downward until one of the next nodes hits the
                    # seafloor, or is deeper than deep_bound
                    integral = 0
                    while True:
                        if nodes[0].below is None or nodes[1].below is None or nodes[2].below is None:
                            # We've reached the seafloor
                            last_nodes = None
                            break
                        if nodes[0].below.depth > deep_bound or nodes[1].below.depth > deep_bound or nodes[2].below.depth > deep_bound:
                            # At least one of the next nodes will pass
                            # deep_bound; save the current nodes
                            last_nodes = []
                            for i in range(3):
                                last_nodes.append(nodes[i])
                            break
                        else:
                            # Calculate mean of data at six corners of this
                            # triangular prism, and mean depths at three edges
                            values_tmp = []
                            dz_tmp = []
                            for i in range(3):
                                values_tmp.append(data[nodes[i].id])
                                values_tmp.append(data[nodes[i].below.id])
                                dz_tmp.append(abs(nodes[i].depth - nodes[i].below.depth))
                                # Get ready for next iteration of while loop
                                nodes[i] = nodes[i].below
                            # Integrand is mean of data at corners *
                            # area of upper face * mean of depths at edges
                            integral += mean(array(values_tmp))*area*mean(array(dz_tmp))

                    # Now integrate from shallow_bound to first_nodes by
                    # linearly interpolating each node to shallow_bound
                    values_tmp = []
                    dz_tmp = []
                    for i in range(3):
                        values_tmp.append(data[first_nodes[i].id])
                        id1, id2, coeff1, coeff2 = elm.nodes[i].find_depth(shallow_bound)
                        if any(isnan(array([id1, id2, coeff1, coeff2]))):
                            # first_nodes[i] was the shallowest node, we can't
                            # interpolate above it
                            values_tmp.append(NaN)
                        else:
                            values_tmp.append(coeff1*data[id1] + coeff2*data[id2])
                        dz_tmp.append(abs(first_nodes[i].depth - shallow_bound))
                    if any(isnan(array(values_tmp))):
                        pass
                    else:
                        integral += mean(array(values_tmp))*area*mean(array(dz_tmp))  

                    # Now integrate from last_nodes to deep_bound by linearly
                    # interpolating each node to shallow_bound, unless we hit
                    # the seafloor earlier
                    if last_nodes is not None:
                        values_tmp = []
                        dz_tmp = []
                        for i in range(3):
                            values_tmp.append(data[last_nodes[i].id])
                            id1, id2, coeff1, coeff2 = elm.nodes[i].find_depth(deep_bound)
                            if any(isnan(array([id1, id2, coeff1, coeff2]))):
                                # last_nodes[i] was the deepest node, we can't
                                # interpolate below it
                                values_tmp.append(NaN)
                            else:
                                values_tmp.append(coeff1*data[id1] + coeff2*data[id2])
                            dz_tmp.append(abs(deep_bound - last_nodes[i].depth))
                        if any(isnan(array(values_tmp))):
                            pass
                        else:
                            integral += mean(array(values_tmp))*area*mean(array(dz_tmp))

                    # All done; divide integral by volume to get the average
                    values.append(integral/volume)                    
                    # Make new patches for elements which exist at this depth
                    coord = transpose(vstack((elm.x, elm.y)))
                    plot_patches.append(Polygon(coord, True, linewidth=0.))

    if depth_key < 3:
        # Use all patches
        plot_patches = patches[:]

    if mask_cavities:
        # Get mask array of patches for ice shelf cavity elements
        mask_patches = iceshelf_mask(elements)
        if var_name == 'wnet':
            # Swap with regular patches so that open ocean elements are masked,
            # ice shelf cavity nodes are not
            tmp = plot_patches
            plot_patches = mask_patches
            mask_patches = tmp            

    # Set up plot
    if circumpolar:
        fig = figure(figsize=(16,12))
        ax = fig.add_subplot(1,1,1, aspect='equal')
    else:
        fig = figure(figsize=(16,8))
        ax = fig.add_subplot(1,1,1)
    # Set colourmap for patches, and refer it to the values array
    img = PatchCollection(plot_patches, cmap=jet)
    img.set_array(array(values))
    img.set_edgecolor('face')
    # Add patches to plot
    ax.add_collection(img)
    if mask_cavities:
        # Set colour to light grey for patches in mask
        overlay = PatchCollection(mask_patches, facecolor=(0.6, 0.6, 0.6))
        overlay.set_edgecolor('face')
        # Add mask to plot
        ax.add_collection(overlay)

    # Configure plot
    if circumpolar:
        xlim([-lat_max, lat_max])
        ylim([-lat_max, lat_max])
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
    else:
        xlim([lon_min, lon_max])
        ylim([lat_min, lat_max])
        xticks(lon_ticks)
        yticks(lat_ticks)
        xlabel('Longitude', fontsize=font_sizes[1])
        ylabel('Latitude', fontsize=font_sizes[1])
        setp(ax.get_xticklabels(), fontsize=font_sizes[2])
        setp(ax.get_xticklabels(), fontsize=font_sizes[2])
    title(name + ' (' + units + ') ' + depth_string, fontsize=font_sizes[0])    
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=font_sizes[2])
    if set_limits:
        img.set_clim(vmin=limits[0], vmax=limits[1])

    # Plot specified points
    #problem_ids = [4415, 4431, 4432, 6130]
    #problem_x = []
    #problem_y = []
    #for elm in elements:
        #for i in range(3):
            #if elm.nodes[i].id in problem_ids:
                #problem_x.append(elm.x[i])
                #problem_y.append(elm.y[i])
                #problem_ids.remove(elm.nodes[i].id)
    #ax.plot(problem_x, problem_y, 'or')    

    if save:
        savefig(fig_name)
    else:
        show()
