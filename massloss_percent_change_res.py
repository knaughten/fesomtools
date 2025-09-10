from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.pyplot import *
from matplotlib.colors import LinearSegmentedColormap
from patches import *

# Make a circumpolar plot of the percentage change in mass loss for each major
# ice shelf in the high-res spinup (1992-2005 average, rep3) compared to
# low-res.
def massloss_percent_change_res ():

    # Paths to low-res mesh directory and both output directories
    directory_head = '/short/y99/kaa561/FESOM/'
    mesh_path = directory_head + 'mesh/low_res/'
    expt_paths = ['lowres_spinup/', 'highres_spinup/']
    # Skip the first two repetitions
    skipyears = 28
    num_years = 14
    # 5-day averages
    peryear = 365//5

    # Limits on longitude and latitude for each ice shelf
    # These depend on the source geometry, in this case RTopo 1.05
    # Note there is one extra index at the end of each array; this is because
    # the Ross region crosses the line 180W and therefore is split into two
    # We have -181 and 181 not -180 and 180 at this boundary so that
    # elements which cross the boundary are still counted
    lon_min = [-62.67, -65.5, -79.17, -85, -104.17, -102.5, -108.33, -114.5, -135.67, -149.17, -155, 144, 115, 94.17, 80.83, 65, 33.83, 19, 12.9, 9.33, -10.05, -28.33, -181, 158.33]
    lon_max = [-59.33, -60, -66.67, -28.33, -88.83, -99.17, -103.33, -111.5, -114.33, -140, -145, 146.62, 123.33, 102.5, 89.17, 75, 37.67, 33.33, 16.17, 12.88, 7.6, -10.33, -146.67, 181]
    lat_min = [-73.03, -69.35, -74.17, -83.5, -73.28, -75.5, -75.5, -75.33, -74.9, -76.42, -78, -67.83, -67.17, -66.67, -67.83, -73.67, -69.83, -71.67, -70.5, -70.75, -71.83, -76.33, -85, -84.5]
    lat_max = [-69.37, -66.13, -69.5, -74.67, -71.67, -74.17, -74.67, -73.67, -73, -75.17, -76.41, -66.67, -66.5, -64.83, -66.17, -68.33, -68.67, -68.33, -69.33, -69.83, -69.33, -71.5, -77.77, -77]
    num_shelves = len(lon_min)-1

    # Plotting parameters
    max_lat_plot = -63+90
    mask_cavities = True
    circumpolar = True

    # Build FESOM mesh, just for low res
    # Get separate patches for the open ocean and minor ice shelf elements
    # so we can mask them out
    elements, mask_patches = make_patches(mesh_path, circumpolar, mask_cavities, only_major=True)
    patches = iceshelf_mask(elements, only_major=True)

    # Set up array for mass loss at each ice shelf, 1992-2005 average, for each
    # simulation
    massloss = empty([2, num_shelves])
    for expt in range(2):
        # Read log file
        f = open(directory_head + expt_paths[expt] + '/massloss.log', 'r')
        f.readline()
        for line in f:
            try:
                tmp = float(line)
            except(ValueError):
                break
        for index in range(num_shelves):
            massloss_tmp = []
            for line in f:
                try:
                    massloss_tmp.append(float(line))
                except(ValueError):
                    break
            # Calculate average over third repetition of 1992-2005
            massloss_tmp = mean(array(massloss_tmp[skipyears*peryear:(skipyears+num_years)*peryear]))
            massloss[expt,index] = massloss_tmp

    # Calculate percent change
    percent_change = (massloss[1,:] - massloss[0,:])/massloss[0,:]*100
    # Save these values to the right Elements
    values = []
    for elm in elements:
        # Make sure we're actually in an ice shelf cavity
        if elm.cavity:
            keep = False
            # Loop over ice shelves
            for index in range(num_shelves):
                # Figure out whether or not this element is part of the given
                # ice shelf
                if all(elm.lon >= lon_min[index]) and all(elm.lon <= lon_max[index]) and all(elm.lat >= lat_min[index]) and all(elm.lat <= lat_max[index]):
                    keep = True
                    tmp = percent_change[index]
                if index == num_shelves-1:
                    # Ross region is split into two
                    if all(elm.lon >= lon_min[index+1]) and all(elm.lon <= lon_max[index+1]) and all(elm.lat >= lat_min[index+1]) and all(elm.lat <= lat_max[index+1]):
                        keep = True
                        tmp = percent_change[index]
            if keep:
                values.append(tmp)                

    # Set up a grey square covering the domain, anything that isn't covered
    # up later is land
    x_reg, y_reg = meshgrid(linspace(-max_lat_plot, max_lat_plot, num=100), linspace(-max_lat_plot, max_lat_plot, num=100))
    land_square = zeros(shape(x_reg))

    # Make contour lines for ice shelf front
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

    # Make custom colour scale
    values = array(values)
    cmap_vals = array([amin(values), 0, 40, 80, amax(values)])
    cmap_colours = [(0.26,0.45,0.86),(1,1,1),(1,0.9,0.4),(0.99,0.59,0.18),(0.5,0.0,0.08)]
    cmap_vals_norm = (cmap_vals-amin(values))/(amax(values)-amin(values))
    cmap_list = []
    for i in range(size(cmap_vals)):
        cmap_list.append((cmap_vals_norm[i], cmap_colours[i]))
    my_cmap = LinearSegmentedColormap.from_list('my_cmap',cmap_list)

    # Plot
    fig = figure(figsize=(16,12))
    ax = fig.add_subplot(1,1,1, aspect='equal')
    # Start with grey square background for land
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    # Add ice shelf elements
    img = PatchCollection(patches, cmap=my_cmap)
    img.set_array(values)
    img.set_edgecolor('face')
    ax.add_collection(img)
    # Mask out the open ocean in white
    overlay = PatchCollection(mask_patches, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    # Add ice shelf front contours
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax.add_collection(contours)
    xlim([-max_lat_plot, max_lat_plot])
    ylim([-max_lat_plot, max_lat_plot])
    axis('off')
    title('% Change in Ice Shelf Mass Loss (high res vs low res)', fontsize=30)
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=20)

    #fig.show()
    fig.savefig('massloss_percent_change_res.png')
    

# Command-line interface
if __name__ == "__main__":

    massloss_percent_change_res()

    


    
