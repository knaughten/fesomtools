from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.cm import *
from matplotlib.colors import LinearSegmentedColormap
from patches import *
from fesom_grid import *

# This script needs python 2.7.6


# For the given variable, find the min and max values across the 2 simulations
# in the given reigon.
# Input:
# data_old, data_new = data for each 2D element on the old and new meshes
#                      respectively
# x_min, x_max, y_min, y_max = bounds on x and y (using polar coordinate
#                              transformation) for the desired region
# cavity = optional boolean indicating to only consider values in ice shelf
#          cavities (default True)
# Output:
# var_min, var_max = min and max data values in this region across both
#                    simulations
def get_min_max (data_old, data_new, x_min, x_max, y_min, y_max, cavity=True):

    # Start with something impossible
    var_min = amax(data_old)
    var_max = amin(data_old)
    # Modify with old data
    i = 0
    for elm in elements_old:
        if (not cavity) or (cavity and elm.cavity):
            if any(elm.x >= x_min) and any(elm.x <= x_max) and any(elm.y >= y_min) and any(elm.y <= y_max):
                if data_old[i] < var_min:
                    var_min = data_old[i]
                if data_old[i] > var_max:
                    var_max = data_old[i]
            i += 1
    # Modify with new data
    i = 0
    for elm in elements_new:
        if (not cavity) or (cavity and elm.cavity):
            if any(elm.x >= x_min) and any(elm.x <= x_max) and any(elm.y >= y_min) and any(elm.y <= y_max):
                if data_new[i] < var_min:
                    var_min = data_new[i]
                if data_new[i] > var_max:
                    var_max = data_new[i]
            i += 1

    return var_min, var_max


# Plot bathymetry for both simulations for the given region, in the given axis.
# Input:
# x_min, x_max, y_min, y_max = bounds on x and y (using polar coordinate
#               transformation) for the desired region
# gs = GridSpec object of size 1x2 to plot in
# cbaxes = Axes object for location of colourbar
# cbar_ticks = 1D array containing values for ticks on colourbar
# letter = 'a', 'b', 'c', etc. to add before the bathymetry title, for use in a
#          figure showing multiple variables
# y0 = y-coordinate of model titles for the entire plot, assuming bathymetry is
#      at the top (i.e. letter='a'). Play around between 1.15 and 1.35.
def plot_bathy (x_min, x_max, y_min, y_max, gs, cbaxes, cbar_ticks, letter, y0):

    # Set up a grey square to fill the background with land
    x_reg, y_reg = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
    land_square = zeros(shape(x_reg))
    # Find bounds on variable in this region
    var_min, var_max = get_min_max(bathy_old, bathy_new, x_min, x_max, y_min, y_max, cavity=False)
    print('Bounds on bathymetry: ' + str(var_min) + ' ' + str(var_max))

    # Old simulation
    ax = subplot(gs[0,0], aspect='equal')
    # Start with land background
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    # Add ocean elements
    img = PatchCollection(patches_all_old, cmap='jet')
    img.set_array(array(bathy_old))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    # Add ice shelf front contour lines
    contours_old = LineCollection(contour_lines_old, edgecolor='black', linewidth=1)
    ax.add_collection(contours_old)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    # Model label
    text(0.5, y0, 'Old mesh', fontsize=16, horizontalalignment='center', transform=ax.transAxes)
    # Variable title
    title(letter + ') Bathymetry (m)', fontsize=18, loc='left')

    # New simulation
    ax = subplot(gs[0,1], aspect='equal')
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches_all_new, cmap='jet')
    img.set_array(array(bathy_new))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    contours_new = LineCollection(contour_lines_new, edgecolor='black', linewidth=1)
    ax.add_collection(contours_new)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    text(0.5, y0, 'New mesh', fontsize=16, horizontalalignment='center', transform=ax.transAxes)
    # Colourbar
    cbar = colorbar(img, cax=cbaxes, ticks=cbar_ticks)


# Plot water column thickness for both simulations for the given region, in the
# given axis.
# Input:
# x_min, x_max, y_min, y_max = bounds on x and y (using polar coordinate
#               transformation) for the desired region
# gs = GridSpec object of size 1x2 to plot in
# cbaxes = Axes object for location of colourbar
# cbar_ticks = 1D array containing values for ticks on colourbar
# letter = 'a', 'b', 'c', etc. to add before the wct title, for use in a figure
#          showing multiple variables
# y0 = y-coordinate of model titles for the entire plot, assuming wct is at the
#      top (i.e. letter='a'). Play around between 1.15 and 1.35.
def plot_wct (x_min, x_max, y_min, y_max, gs, cbaxes, cbar_ticks, letter, y0):

    # Set up a grey square to fill the background with land
    x_reg, y_reg = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
    land_square = zeros(shape(x_reg))
    # Find bounds on variable in this region
    var_min, var_max = get_min_max(wct_old, wct_new, x_min, x_max, y_min, y_max, cavity=True)
    print('Bounds on wct: ' + str(var_min) + ' ' + str(var_max))

    # Old simulation
    ax = subplot(gs[0,0], aspect='equal')
    # Start with land background
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    # Add ice shelf elements
    img = PatchCollection(patches_old, cmap='jet')
    img.set_array(array(wct_old))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    # Mask out the open ocean in white
    overlay = PatchCollection(mask_patches_old, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    # Model label
    text(0.5, y0, 'Old mesh', fontsize=16, horizontalalignment='center', transform=ax.transAxes)
    # Variable title
    title(letter + ') Water column thickness (m)', fontsize=18, loc='left')

    # New simulation
    ax = subplot(gs[0,1], aspect='equal')
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches_new, cmap='jet')
    img.set_array(array(wct_new))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    overlay = PatchCollection(mask_patches_new, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    text(0.5, y0, 'New mesh', fontsize=16, horizontalalignment='center', transform=ax.transAxes)
    # Colourbar
    cbar = colorbar(img, cax=cbaxes, ticks=cbar_ticks)


# Plot bottom water temperature for both simulations for the given region, in
# the given axis.
# Input:
# x_min, x_max, y_min, y_max = bounds on x and y (using polar coordinate
#               transformation) for the desired region
# gs = GridSpec object of size 1x2 to plot in
# cbaxes = Axes object for location of colourbar
# cbar_ticks = 1D array containing values for ticks on colourbar
# letter = 'a', 'b', 'c', etc. to add before the bwtemp title, for use in a
#          figure showing multiple variables
def plot_bwtemp (x_min, x_max, y_min, y_max, gs, cbaxes, cbar_ticks, letter):

    # Set up a grey square to fill the background with land
    x_reg, y_reg = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
    land_square = zeros(shape(x_reg))
    # Find bounds on variable in this region
    var_min, var_max = get_min_max(bwtemp_old, bwtemp_new, x_min, x_max, y_min, y_max, cavity=False)
    print('Bounds on bottom temperature: ' + str(var_min) + ' ' + str(var_max))

    # Old simulation
    ax = subplot(gs[0,0], aspect='equal')
    # Start with land background
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    # Add ocean elements
    img = PatchCollection(patches_all_old, cmap='jet')
    img.set_array(array(bwtemp_old))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    # Add ice shelf front contour lines
    contours_old = LineCollection(contour_lines_old, edgecolor='black', linewidth=1)
    ax.add_collection(contours_old)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    # Variable title
    title(letter + r') Bottom water temperature ($^{\circ}$C)', fontsize=18, loc='left')

    # New simulation
    ax = subplot(gs[0,1], aspect='equal')
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches_all_new, cmap='jet')
    img.set_array(array(bwtemp_new))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    contours_new = LineCollection(contour_lines_new, edgecolor='black', linewidth=1)
    ax.add_collection(contours_new)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    # Colourbar
    cbar = colorbar(img, cax=cbaxes, ticks=cbar_ticks)


# Plot ice shelf melt rate for both simulations for the given region, in the
# given axis.
# Input:
# x_min, x_max, y_min, y_max = bounds on x and y (using polar coordinate
#               transformation) for the desired region
# gs = GridSpec object of size 1x2 to plot in
# cbaxes = Axes object for location of colourbar
# cbar_ticks = 1D array containing values for ticks on colourbar
# change_points = list of size 3 containing values where the colourmap should
#                 transition (1) from yellow to orange, (2) from orange to red,
#                 (3) from red to magenta. Should not include the minimum value,
#                 0, or the maximum value. This way the custom colourmap can be 
#                 adjusted so that all melt rates are visible, particularly
#                 for ice shelves with strong spatial variations in melt.
# letter = 'a', 'b', 'c', etc. to add before the melt rate title, for use in a
#          figure showing multiple variables
def plot_melt (x_min, x_max, y_min, y_max, gs, cbaxes, cbar_ticks, change_points, letter):

    # Set up a grey square to fill the background with land
    x_reg, y_reg = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
    land_square = zeros(shape(x_reg))
    # Find bounds on variable in this region
    var_min, var_max = get_min_max(melt_old, melt_new, x_min, x_max, y_min, y_max, cavity=True)
    print('Bounds on melt: ' + str(var_min) + ' ' + str(var_max))
    # Special colour map
    if var_min < 0:
        # There is refreezing here; include blue for elements < 0
        cmap_vals = array([var_min, 0, change_points[0], change_points[1], change_points[2], var_max])
        cmap_colors = [(0.26, 0.45, 0.86), (1, 1, 1), (1, 0.9, 0.4), (0.99, 0.59, 0.18), (0.5, 0.0, 0.08), (0.96, 0.17, 0.89)]
        cmap_vals_norm = (cmap_vals - var_min)/(var_max - var_min)
        cmap_vals_norm[-1] = 1
        cmap_list = []
        for i in range(size(cmap_vals)):
            cmap_list.append((cmap_vals_norm[i], cmap_colors[i]))
        mf_cmap = LinearSegmentedColormap.from_list('melt_freeze', cmap_list)
    else:
        # No refreezing
        cmap_vals = array([0, change_points[0], change_points[1], change_points[2], var_max])
        cmap_colors = [(1, 1, 1), (1, 0.9, 0.4), (0.99, 0.59, 0.18), (0.5, 0.0, 0.08), (0.96, 0.17, 0.89)]
        cmap_vals_norm = cmap_vals/var_max
        cmap_vals_norm[-1] = 1
        cmap_list = []
        for i in range(size(cmap_vals)):
            cmap_list.append((cmap_vals_norm[i], cmap_colors[i]))
        mf_cmap = LinearSegmentedColormap.from_list('melt_freeze', cmap_list)

    # Old simulation
    ax = subplot(gs[0,0], aspect='equal')
    # Start with land background
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    # Add ice shelf elements
    img = PatchCollection(patches_old, cmap=mf_cmap)
    img.set_array(array(melt_old))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    # Mask out the open ocean in white
    overlay = PatchCollection(mask_patches_old, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    # Variable title
    title(letter + ') Ice shelf melt rate (m/y)', fontsize=18, loc='left')

    # New simulation
    ax = subplot(gs[0,1], aspect='equal')
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches_new, cmap=mf_cmap)
    img.set_array(array(melt_new))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    overlay = PatchCollection(mask_patches_new, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    # Colourbar
    cbar = colorbar(img, cax=cbaxes, ticks=cbar_ticks)


#***********MAIN PROCESSING***********g
mesh_path_old = '/short/y99/kaa561/FESOM/mesh/high_res/'
mesh_path_new = '/short/y99/kaa561/FESOM/mesh/meshB/'
ocn_file_old = '/short/y99/kaa561/FESOM/intercomparison_highres_old/oce_2002_2016_avg.nc'
ocn_file_new = '/short/y99/kaa561/FESOM/intercomparison_highres/output/oce_2002_2016_avg.nc'
ice_file_old = '/short/y99/kaa561/FESOM/intercomparison_highres_old/wnet_2002_2016_avg.nc'
ice_file_new = '/short/y99/kaa561/FESOM/intercomparison_highres/output/wnet_2002_2016_avg.nc'
# Constants
sec_per_year = 365.25*24*3600

print('Building old mesh')
# Mask open ocean
elements_old, mask_patches_old = make_patches(mesh_path_old, circumpolar=True, mask_cavities=True)
# Unmask ice shelves
patches_old = iceshelf_mask(elements_old)
# Also make a set of patches with open ocean unmasked
patches_all_old = []
for elm in elements_old:
    coord = transpose(vstack((elm.x, elm.y)))
    patches_all_old.append(Polygon(coord, True, linewidth=0.))
print('Building new mesh')
elements_new, mask_patches_new = make_patches(mesh_path_new, circumpolar=True, mask_cavities=True)
patches_new = iceshelf_mask(elements_new)
patches_all_new = []
for elm in elements_new:
    coord = transpose(vstack((elm.x, elm.y)))
    patches_all_new.append(Polygon(coord, True, linewidth=0.))

print('Building ice shelf front contours')
contour_lines_old = []
for elm in elements_old:
    # Select elements where exactly 2 of the 3 nodes are in a cavity
    if count_nonzero(elm.cavity_nodes) == 2:
        # Save the coastal flags and x- and y- coordinates of these 2
        coast_tmp = []
        x_tmp = []
        y_tmp = []
        for i in range(3):
            if elm.cavity_nodes[i]:
                coast_tmp.append(elm.coast_nodes[i])
                x_tmp.append(elm.x[i])
                y_tmp.append(elm.y[i])
        # Select elements where at most 1 of these 2 nodes are coastal
        if count_nonzero(coast_tmp) < 2:
            # Draw a line between the 2 nodes
            contour_lines_old.append([(x_tmp[0], y_tmp[0]), (x_tmp[1], y_tmp[1])])
contour_lines_new = []
for elm in elements_new:
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
            contour_lines_new.append([(x_tmp[0], y_tmp[0]), (x_tmp[1], y_tmp[1])])

print('Calculating bathymetry')
# Depth of bottom layer, averaged over 3 corners
bathy_old = []
for elm in elements_old:
    bathy_old.append(mean([elm.nodes[0].find_bottom().depth, elm.nodes[1].find_bottom().depth, elm.nodes[2].find_bottom().depth]))
bathy_new = []
for elm in elements_new:
    bathy_new.append(mean([elm.nodes[0].find_bottom().depth, elm.nodes[1].find_bottom().depth, elm.nodes[2].find_bottom().depth]))

print('Calculating water column thickness')
# Depth of bottom layer minus depth of surface layer in ice shelf cavities,
# averaged over 3 corners
wct_old = []
for elm in elements_old:
    if elm.cavity:
        wct_old.append(mean([elm.nodes[0].find_bottom().depth - elm.nodes[0].depth, elm.nodes[1].find_bottom().depth - elm.nodes[1].depth, elm.nodes[2].find_bottom().depth - elm.nodes[2].depth]))
wct_new = []
for elm in elements_new:
    if elm.cavity:
        wct_new.append(mean([elm.nodes[0].find_bottom().depth - elm.nodes[0].depth, elm.nodes[1].find_bottom().depth - elm.nodes[1].depth, elm.nodes[2].find_bottom().depth - elm.nodes[2].depth]))

print('Calculating bottom water temperature')
# Read full 3D field to start
id = Dataset(ocn_file_old, 'r')
node_temp_old = id.variables['temp'][0,:]
id.close()
# For each element, average over 3 corners
bwtemp_old = []
for elm in elements_old:
    bwtemp_old.append(mean([node_temp_old[elm.nodes[0].find_bottom().id], node_temp_old[elm.nodes[1].find_bottom().id], node_temp_old[elm.nodes[2].find_bottom().id]]))
id = Dataset(ocn_file_new, 'r')
node_temp_new = id.variables['temp'][0,:]
id.close()
bwtemp_new = []
for elm in elements_new:
    bwtemp_new.append(mean([node_temp_new[elm.nodes[0].find_bottom().id], node_temp_new[elm.nodes[1].find_bottom().id], node_temp_new[elm.nodes[2].find_bottom().id]]))

print('Calculating ice shelf melt rate')
# Read melt rate at 2D nodes and convert from m/s to m/y
id = Dataset(ice_file_old, 'r')
node_melt_old = id.variables['wnet'][0,:]*sec_per_year
id.close()
# For each element in an ice shelf cavity, average over 3 corners
melt_old = []
for elm in elements_old:
    if elm.cavity:
        melt_old.append(mean([node_melt_old[elm.nodes[0].id], node_melt_old[elm.nodes[1].id], node_melt_old[elm.nodes[2].id]]))
id = Dataset(ice_file_new, 'r')
node_melt_new = id.variables['wnet'][0,:]*sec_per_year
id.close()
melt_new = []
for elm in elements_new:
    if elm.cavity:
        melt_new.append(mean([node_melt_new[elm.nodes[0].id], node_melt_new[elm.nodes[1].id], node_melt_new[elm.nodes[2].id]]))

print('Plotting Pine Island Ice Shelf')
x_min_tmp = -15.9
x_max_tmp = -14.2
y_min_tmp = -3.8
y_max_tmp = -2.25
fig = figure(figsize=(6,8))
fig.patch.set_facecolor('white')
# Bathymetry
gs_a = GridSpec(1,2)
gs_a.update(left=0.05, right=0.87, bottom=0.63, top=0.88, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.88, 0.68, 0.025, 0.15])
cbar_ticks = arange(300, 900+200, 200)
plot_bathy(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_a, cbaxes_tmp, cbar_ticks, 'a', 1.15)
# Bottom water temperature
gs_b = GridSpec(1,2)
gs_b.update(left=0.05, right=0.87, bottom=0.33, top=0.58, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.88, 0.38, 0.025, 0.15])
cbar_ticks = arange(-1.8, -1.2+0.2, 0.2)
plot_bwtemp(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_b, cbaxes_tmp, cbar_ticks, 'b')
# Melt
gs_c = GridSpec(1,2)
gs_c.update(left=0.05, right=0.87, bottom=0.03, top=0.28, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.88, 0.08, 0.025, 0.15])
cbar_ticks = arange(0, 4+1, 1)
plot_melt(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_c, cbaxes_tmp, cbar_ticks, [1, 2, 4], 'c')
suptitle('Pine Island Ice Shelf', fontsize=24)
fig.show()
fig.savefig('bugs_fesom_pig.png')

print('Plotting Amery Ice Shelf')
x_min_tmp = 15.25
x_max_tmp = 20.5
y_min_tmp = 4.75
y_max_tmp = 8
fig = figure(figsize=(8,8))
fig.patch.set_facecolor('white')
# Water column thickness
gs_a = GridSpec(1,2)
gs_a.update(left=0.05, right=0.89, bottom=0.63, top=0.88, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.9, 0.68, 0.025, 0.15])
cbar_ticks = arange(300, 1200+300, 300)
plot_wct(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_a, cbaxes_tmp, cbar_ticks, 'a', 1.15)
# Bottom water temperature
gs_b = GridSpec(1,2)
gs_b.update(left=0.05, right=0.89, bottom=0.33, top=0.58, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.9, 0.38, 0.025, 0.15])
cbar_ticks = arange(-2.3, -1.4+0.3, 0.3)
plot_bwtemp(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_b, cbaxes_tmp, cbar_ticks, 'b')
# Melt
gs_c = GridSpec(1,2)
gs_c.update(left=0.05, right=0.89, bottom=0.03, top=0.28, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.9, 0.08, 0.025, 0.15])
cbar_ticks = arange(0, 12+4, 4)
plot_melt(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_c, cbaxes_tmp, cbar_ticks, [1, 3, 8], 'c')
suptitle('Amery Ice Shelf', fontsize=24)
fig.show()
fig.savefig('bugs_fesom_amery.png')
