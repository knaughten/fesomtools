from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.cm import *
from matplotlib.colors import LinearSegmentedColormap
from patches import *
from unrotate_grid import *
from in_triangle import *

# This script needs python 2.7.6

# Truncate colourmap function from https://stackoverflow.com/questions/40929467/how-to-use-and-plot-only-a-part-of-a-colorbar-in-matplotlib
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=-1):
    if n== -1:
        n = cmap.N
    new_cmap = LinearSegmentedColormap.from_list('trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval), cmap(linspace(minval, maxval, n)))
    return new_cmap


# For the given variable over the given region, find the min and max values
# at the beginning of the simulation, as well as the max absolute value of the
# anomalies at the end of each scenario.
# Input:
# data_beg = data for each 2D element on the FESOM mesh (if cavity=True, cavity
#            elements only) for the beginning of the simulation
# data_diff = array of size (num_expts) x (num_elm or num_cavity_elm) for
#             anomalies at each simulation
# x_min, x_max, y_min, y_max = bounds on x and y (using polar coordinate
#               transformation from above) for the desired region
# cavity = optional boolean indicating to only consider values in ice shelf
#          cavities (default True)
# Output:
# var_min, var_max = min and max data values in this region across data_beg
# diff_min, diff_max = min and max values of data_diff in this region for all
#                      scenarios
def get_min_max (data_beg, data_diff, x_min, x_max, y_min, y_max, cavity=True):

    # Start with min and max of beginning
    # Initialise with something impossible
    var_min = amax(array(data_beg))
    var_max = amin(array(data_beg))
    diff_max = 0
    # Modify as needed
    i = 0
    for elm in elements:
        if (not cavity) or (cavity and elm.cavity):
            if any(elm.x >= x_min) and any(elm.x <= x_max) and any(elm.y >= y_min) and any(elm.y <= y_max):
                if data_beg[i] < var_min:
                    var_min = data_beg[i]
                if data_beg[i] > var_max:
                    var_max = data_beg[i]
            i += 1
    # Now find min and max difference across all scenarios
    diff_min = 0
    diff_max = 0
    for expt in range(num_expts):
        i = 0
        for elm in elements:
            if (not cavity) or (cavity and elm.cavity):
                if any(elm.x >= x_min) and any(elm.x <= x_max) and any(elm.y >= y_min) and any(elm.y <= y_max):
                    if data_diff[expt,i] < diff_min:
                        diff_min = data_diff[expt,i]
                    if data_diff[expt,i] > diff_max:
                        diff_max = data_diff[expt,i]
                i += 1

    return var_min, var_max, diff_min, diff_max


# Interpolate bathymetry to a regular grid in the given region so that we can
# contour isobaths later.
# Input:
# lon_min, lon_max = bounds on longitude for regular grid (-180 to 180)
# lat_min, lat_max = bounds on latitude (-90 to 90)
# Output:
# x_lonlat_reg, y_lonlat_reg = polar coordinates of regular grid for plotting
# bathy_reg = bathymetry interpolated to the regular grid
def interpolate_bathy (lon_min, lon_max, lat_min, lat_max):

    # Size of regular grid
    num_lon = 500
    num_lat = 500

    # Set up regular grid
    # Start with boundaries
    lon_reg_edges = linspace(lon_min, lon_max, num_lon+1)
    lat_reg_edges = linspace(lat_min, lat_max, num_lat+1)
    # Now get centres
    lon_reg = 0.5*(lon_reg_edges[:-1] + lon_reg_edges[1:])
    lat_reg = 0.5*(lat_reg_edges[:-1] + lat_reg_edges[1:])
    # Make 2D versions
    lon_reg_2d, lat_reg_2d = meshgrid(lon_reg, lat_reg)
    # Calculate polar coordinates transformation for plotting
    x_lonlat_reg = -(lat_reg_2d+90)*cos(lon_reg_2d*deg2rad+pi/2)
    y_lonlat_reg = (lat_reg_2d+90)*sin(lon_reg_2d*deg2rad+pi/2)
    # Set up array for bathymetry on this regular grid
    bathy_reg = zeros([num_lat, num_lon])

    # For each element, check if a point on the regular lat-lon grid lies
    # within. If so, do barycentric interpolation to that point, of bottom
    # node depths (i.e. bathymetry).
    for elm in elements:
        # Ignore ice shelf cavities
        if not elm.cavity:
            # Check if we are within domain of regular grid
            if amax(elm.lon) > lon_min and amin(elm.lon) < lon_max and amax(elm.lat) > lat_min and amin(elm.lat) < lat_max:
                # Find largest regular longitude value west of Element
                tmp = nonzero(lon_reg > amin(elm.lon))[0]
                if len(tmp) == 0:
                    # Element crosses the western boundary
                    iW = 0
                else:
                    iW = tmp[0] - 1
                # Find smallest regular longitude value east of Element
                tmp = nonzero(lon_reg > amax(elm.lon))[0]
                if len(tmp) == 0:
                    # Element crosses the eastern boundary
                    iE = num_lon
                else:
                    iE = tmp[0]
                # Find largest regular latitude value south of Element
                tmp = nonzero(lat_reg > amin(elm.lat))[0]
                if len(tmp) == 0:
                    # Element crosses the southern boundary
                    jS = 0
                else:
                    jS = tmp[0] - 1
                # Find smallest regular latitude value north of Element
                tmp = nonzero(lat_reg > amax(elm.lat))[0]
                if len(tmp) == 0:
                    # Element crosses the northern boundary
                    jN = num_lat
                else:
                    jN = tmp[0]
                for i in range(iW+1,iE):
                    for j in range(jS+1,jN):
                        # There is a chance that the regular gridpoint at (i,j)
                        # lies within this element
                        lon0 = lon_reg[i]
                        lat0 = lat_reg[j]
                        if in_triangle(elm, lon0, lat0):
                            # Yes it does
                            # Get area of entire triangle
                            area = triangle_area(elm.lon, elm.lat)
                            # Get area of each sub-triangle formed by
                            # (lon0, lat0)
                            area0 = triangle_area([lon0, elm.lon[1], elm.lon[2]], [lat0, elm.lat[1], elm.lat[2]])
                            area1 = triangle_area([lon0, elm.lon[0], elm.lon[2]], [lat0, elm.lat[0], elm.lat[2]])
                            area2 = triangle_area([lon0, elm.lon[0], elm.lon[1]], [lat0, elm.lat[0], elm.lat[1]])
                            # Find fractional area of each
                            cff = [area0/area, area1/area, area2/area]
                            # Find bottom depth of each node
                            bathy_vals = []
                            for n in range(3):
                                bathy_vals.append(elm.nodes[n].find_bottom().depth)
                            # Barycentric interpolation to lon0, lat0
                            bathy_reg[j,i] = sum(array(cff)*array(bathy_vals))
    # Mask out points which are identically zero (land and ice shelves)
    bathy_reg = ma.masked_where(bathy_reg==0, bathy_reg)

    return x_lonlat_reg, y_lonlat_reg, bathy_reg


# Plot ice shelf melt rate in the given axes: absolute melt rate for the 
# beginning of the simulation, followed by anomalies at the end of each
# experiment.
# Input:
# x_min, x_max, y_min, y_max = bounds on x and y (using polar coordinate
#               transformation from above) for the desired region
# gs = GridSpec object of size 1x3 to plot in
# cbaxes1, cbaxes2 = Axes objects for location of colourbars (one for
#                    absolute melt, one for anomalies)
# ticks1, ticks2 = lists of corresponding values for colourbar ticks
# change_points = list of size 3 containing values where the colourmap should
#                 transition (1) from yellow to orange, (2) from orange to red,
#                 (3) from red to magenta. Should not include the minimum value,
#                 0, or the maximum value. This way the custom colourmap can be 
#                 adjusted so that all melt rates are visible, particularly
#                 for ice shelves with strong spatial variations in melt.
# letter = 'a', 'b', 'c', etc. to add before the ice shelf melt rate title, for
#          use in a figure showing multiple variables
# y0 = y-coordinate of model titles for the entire plot, assuming melt rate is
#      always at the top (i.e. letter='a'). Play around between 1.15 and 1.35.
# set_diff_max = optional float to use as maximum for anomaly colour scale
# lon_plot = optional list of longitudes to overlay as dotted lines on the
#            first panel (-180 to 180)
# lat_lines = optional list of latitudes to overlay (-90 to 90)
def plot_melt (x_min, x_max, y_min, y_max, gs, cbaxes1, ticks1, cbaxes2, ticks2, change_points, letter, y0, set_diff_max=None, lon_lines=None, lat_lines=None):

    # Set up a grey square to fill the background with land
    x_reg, y_reg = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
    land_square = zeros(shape(x_reg))
    # Also get longitude and latitude values on this regular x-y grid,
    # by inverting the polar coordinate transformation
    lat_xyreg = sqrt(x_reg**2 + y_reg**2) - 90
    lon_xyreg = (arctan2(y_reg, -1*x_reg) - pi/2)/deg2rad

    # Find bounds on melt in this region
    var_min, var_max, diff_min, diff_max = get_min_max(melt_beg, melt_diff, x_min, x_max, y_min, y_max)
    print 'Bounds on melt rate: ' + str(var_min) + ' ' + str(var_max) + ' ' + str(diff_min) + ' ' + str(diff_max)
    if set_diff_max is not None:
        diff_max = set_diff_max
    # Special colour map for absolute melt
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
    # Truncate difference colourmap; make sure we plot 0 though
    diff_bound = max(abs(diff_min), abs(diff_max))
    diff_min = min(diff_min, 0)
    diff_max = max(diff_max, 0)
    min_colour = (diff_min + diff_bound)/(2*diff_bound)
    max_colour = (diff_max + diff_bound)/(2*diff_bound)
    diff_cmap = truncate_colormap(get_cmap('RdBu_r'), min_colour, max_colour)

    # Plot absolute melt at the beginning of the simulation
    ax = subplot(gs[0,0], aspect='equal')
    # Start with land background
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    # Add ice shelf elements
    img = PatchCollection(patches, cmap=mf_cmap)
    img.set_array(array(melt_beg))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    # Mask out the open ocean in white
    overlay = PatchCollection(mask_patches, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    # Overlay longitudes
    contour(x_reg, y_reg, lon_xyreg, lon_lines, colors='black', linestyles='dotted')
    # Overlay latitudes
    contour(x_reg, y_reg, lat_xyreg, lat_lines, colors='black', linestyles='dotted')
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    # Add experiment title for top of plot (melt rate is always on top)
    text(0.5, y0, str(beg_years[0])+'-'+str(beg_years[1]), fontsize=18, horizontalalignment='center', transform=ax.transAxes)
    # Add variable title
    title(letter + ') Ice shelf melt rate (m/y)', loc='left', fontsize=20)
    # Colourbar on the left
    cbar = colorbar(img, cax=cbaxes1, ticks=ticks1)    

    # Plot anomalies for each experiment
    for expt in range(num_expts):
        ax = subplot(gs[0,expt+1], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(patches, cmap=diff_cmap)
        img.set_array(array(melt_diff[expt,:]))
        img.set_edgecolor('face')
        img.set_clim(vmin=diff_min, vmax=diff_max)
        ax.add_collection(img)
        overlay = PatchCollection(mask_patches, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        ax.set_xticks([])
        ax.set_yticks([])
        text(0.5, y0, expt_names[expt], fontsize=18, horizontalalignment='center', transform=ax.transAxes)        
        if expt == num_expts-1:
            # Add subtitle for anomalies
            title(str(end_years[0])+'-'+str(end_years[1])+' anomalies', loc='right', fontsize=18)            
        if expt == num_expts-1:
            # Colourbar on the right
            if set_diff_max is None:
                cbar = colorbar(img, cax=cbaxes2, ticks=ticks2)
            else:
                cbar = colorbar(img, cax=cbaxes2, extend='max', ticks=ticks2)


# Plot bottom water temperature in the given axes: absolute temperature for the 
# beginning of the simulation, followed by anomalies at the end of each
# experiment.
# Input:
# x_min, x_max, y_min, y_max = bounds on x and y (using polar coordinate
#               transformation from above) for the desired region
# gs = GridSpec object of size 1x3 to plot in
# cbaxes1, cbaxes2 = Axes objects for location of colourbars (one for
#                    absolute melt, one for anomalies)
# ticks1, ticks2 = lists of corresponding values for colourbar ticks
# letter = 'a', 'b', 'c', etc. to add before the bottom water temp title, for
#          use in a figure showing multiple variables
# bathy_contour = optional float containing single isobath to contour (positive,
#                 in metres; make sure you call interpolate_bathy first)
def plot_bwtemp (x_min, x_max, y_min, y_max, gs, cbaxes1, ticks1, cbaxes2, ticks2, letter, bathy_contour=None):

    # Set up a grey square to fill the background with land
    x_reg, y_reg = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
    land_square = zeros(shape(x_reg))
    # Find bounds on temperature in this region
    var_min, var_max, diff_min, diff_max = get_min_max(bwtemp_beg, bwtemp_diff, x_min, x_max, y_min, y_max, cavity=False)
    print 'Bounds on bottom water temp: ' + str(var_min) + ' ' + str(var_max) + ' ' + str(diff_min) + ' ' + str(diff_max)
    # Truncate difference colourmap; make sure we plot 0 though
    diff_bound = max(abs(diff_min), abs(diff_max))
    diff_min = min(diff_min, 0)
    diff_max = max(diff_max, 0)
    min_colour = (diff_min + diff_bound)/(2*diff_bound)
    max_colour = (diff_max + diff_bound)/(2*diff_bound)
    diff_cmap = truncate_colormap(get_cmap('RdBu_r'), min_colour, max_colour)

    # Plot absolute temperature at the beginning of the simulation
    ax = subplot(gs[0,0], aspect='equal')
    # Start with land background
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    # Add ocean elements
    img = PatchCollection(patches_all, cmap='jet')
    img.set_array(array(bwtemp_beg))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    # Add ice shelf front contour lines
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax.add_collection(contours)
    if bathy_contour is not None:
        # Overlay dashed contour on regular grid
        contour(x_lonlat_reg, y_lonlat_reg, bathy_reg, levels=[bathy_contour], colors=('black'), linestyles=('dashed'))
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    # Add variable title
    title(letter + r') Bottom water temperature ($^{\circ}$C)', loc='left', fontsize=20)
    # Colourbar on the left
    cbar = colorbar(img, cax=cbaxes1, ticks=ticks1)    

    # Plot anomalies for each experiment
    for expt in range(num_expts):
        ax = subplot(gs[0,expt+1], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(patches_all, cmap=diff_cmap)
        img.set_array(array(bwtemp_diff[expt,:]))
        img.set_edgecolor('face')
        img.set_clim(vmin=diff_min, vmax=diff_max)
        ax.add_collection(img)
        contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
        ax.add_collection(contours)
        if bathy_contour is not None:
            contour(x_lonlat_reg, y_lonlat_reg, bathy_reg, levels=[bathy_contour], colors=('black'), linestyles=('dashed'))
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        ax.set_xticks([])
        ax.set_yticks([])            
        if expt == num_expts-1:
            # Colourbar on the right
            cbar = colorbar(img, cax=cbaxes2, ticks=ticks2)


# Plot bottom water salinity in the given axes: absolute salinity for the 
# beginning of the simulation, followed by anomalies at the end of each
# experiment.
# Input:
# x_min, x_max, y_min, y_max = bounds on x and y (using polar coordinate
#               transformation from above) for the desired region
# gs = GridSpec object of size 1x3 to plot in
# cbaxes1, cbaxes2 = Axes objects for location of colourbars (one for
#                    absolute melt, one for anomalies)
# ticks1, ticks2 = lists of corresponding values for colourbar ticks
# letter = 'a', 'b', 'c', etc. to add before the bottom water salt title, for
#          use in a figure showing multiple variables
# bathy_contour = optional float containing single isobath to contour (positive,
#                 in metres; make sure you call interpolate_bathy first)
def plot_bwsalt (x_min, x_max, y_min, y_max, gs, cbaxes1, ticks1, cbaxes2, ticks2, letter, bathy_contour=None):

    # Set up a grey square to fill the background with land
    x_reg, y_reg = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
    land_square = zeros(shape(x_reg))
    # Find bounds on salinity in this region
    var_min, var_max, diff_min, diff_max = get_min_max(bwsalt_beg, bwsalt_diff, x_min, x_max, y_min, y_max, cavity=False)
    print 'Bounds on bottom water salt: ' + str(var_min) + ' ' + str(var_max) + ' ' + str(diff_min) + ' ' + str(diff_max)
    # Truncate difference colourmap; make sure we plot 0 though
    diff_bound = max(abs(diff_min), abs(diff_max))
    diff_min = min(diff_min, 0)
    diff_max = max(diff_max, 0)
    min_colour = (diff_min + diff_bound)/(2*diff_bound)
    max_colour = (diff_max + diff_bound)/(2*diff_bound)
    diff_cmap = truncate_colormap(get_cmap('RdBu_r'), min_colour, max_colour)

    # Plot absolute salinity at the beginning of the simulation
    ax = subplot(gs[0,0], aspect='equal')
    # Start with land background
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    # Add ocean elements
    img = PatchCollection(patches_all, cmap='jet')
    img.set_array(array(bwsalt_beg))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    # Add ice shelf front contour lines
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax.add_collection(contours)
    if bathy_contour is not None:
        # Overlay dashed contour on regular grid
        contour(x_lonlat_reg, y_lonlat_reg, bathy_reg, levels=[bathy_contour], colors=('black'), linestyles=('dashed'))
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    # Add variable title
    title(letter + ') Bottom water salinity (psu)', loc='left', fontsize=20)
    # Colourbar on the left
    cbar = colorbar(img, cax=cbaxes1, ticks=ticks1)    

    # Plot anomalies for each experiment
    for expt in range(num_expts):
        ax = subplot(gs[0,expt+1], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(patches_all, cmap=diff_cmap)
        img.set_array(array(bwsalt_diff[expt,:]))
        img.set_edgecolor('face')
        img.set_clim(vmin=diff_min, vmax=diff_max)
        ax.add_collection(img)
        contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
        ax.add_collection(contours)
        if bathy_contour is not None:
            contour(x_lonlat_reg, y_lonlat_reg, bathy_reg, levels=[bathy_contour], colors=('black'), linestyles=('dashed'))
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        ax.set_xticks([])
        ax.set_yticks([])            
        if expt == num_expts-1:
            # Colourbar on the right
            cbar = colorbar(img, cax=cbaxes2, ticks=ticks2)


# MAIN PROCESSING
# File paths
mesh_path = '/short/y99/kaa561/FESOM/mesh/meshB/'
directory_beg = '/short/y99/kaa561/FESOM/highres_spinup/'
directories = ['/short/y99/kaa561/FESOM/rcp45_M/', '/short/y99/kaa561/FESOM/rcp45_A/', '/short/y99/kaa561/FESOM/rcp85_M/', '/short/y99/kaa561/FESOM/rcp85_A/', '/short/y99/kaa561/FESOM/highres_spinup/']
forcing_file_beg = 'annual_avg.forcing.diag.1996.2005.nc'
forcing_file_end = 'annual_avg.forcing.diag.2091.2100.nc'
oce_file_beg = 'annual_avg.oce.mean.1996.2005.nc'
oce_file_end = 'annual_avg.oce.mean.2091.2100.nc'
# Titles for plotting
expt_names = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A', 'CONTROL']
num_expts = len(directories)
# Start and end years for each period
beg_years = [1996, 2005]
end_years = [2091, 2100]
# Constants
deg2rad = pi/180.0
sec_per_year = 365.25*24*3600

print 'Building mesh'
# Mask open ocean
elements, mask_patches = make_patches(mesh_path, circumpolar=True, mask_cavities=True)
# Unmask ice shelves
patches = iceshelf_mask(elements)
# Also make a set of patches with open ocean unmasked (for bottom water T/S)
patches_all = []
for elm in elements:
    coord = transpose(vstack((elm.x, elm.y)))
    patches_all.append(Polygon(coord, True, linewidth=0.))
# Count the number of elements in ice shelf cavities
num_cavity_elm = 0
for elm in elements:
    if elm.cavity:
        num_cavity_elm += 1
num_elm = len(elements)

print 'Building ice shelf front contours'
contour_lines = []
for elm in elements:
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
            contour_lines.append([(x_tmp[0], y_tmp[0]), (x_tmp[1], y_tmp[1])])

print 'Calculating ice shelf melt rate'
melt_beg = zeros(num_cavity_elm)
# First read melt rate at beginning
id = Dataset(directory_beg + forcing_file_beg, 'r')
node_melt_beg = id.variables['wnet'][0,:]*sec_per_year
id.close()
# For each element, calculate average over 3 corners
i = 0
for elm in elements:
    if elm.cavity:
        melt_beg[i] = mean([node_melt_beg[elm.nodes[0].id], node_melt_beg[elm.nodes[1].id], node_melt_beg[elm.nodes[2].id]])
        i += 1
# Find anomalies for all other experiments
melt_diff = zeros([num_expts, num_cavity_elm])
for expt in range(num_expts):
    # Read melt rate at end of simulation
    id = Dataset(directories[expt] + forcing_file_end, 'r')
    node_melt_end = id.variables['wnet'][0,:]*sec_per_year
    id.close()
    # Calculate difference
    node_melt_diff = node_melt_end - node_melt_beg
    i = 0
    for elm in elements:
        if elm.cavity:
            melt_diff[expt,i] = mean([node_melt_diff[elm.nodes[0].id], node_melt_diff[elm.nodes[1].id], node_melt_diff[elm.nodes[2].id]])
            i += 1

print 'Calculating bottom water temperature'
bwtemp_beg = zeros(num_elm)
# Read full 3D field to start
id = Dataset(directory_beg + oce_file_beg, 'r')
node_temp_beg = id.variables['temp'][0,:]
id.close()
# For each element, find the surface or bottom nodes and calculate average
i = 0
for elm in elements:
    bwtemp_beg[i] = mean([node_temp_beg[elm.nodes[0].find_bottom().id], node_temp_beg[elm.nodes[1].find_bottom().id], node_temp_beg[elm.nodes[2].find_bottom().id]])
    i += 1
bwtemp_diff = zeros([num_expts, num_elm])
for expt in range(num_expts):
    id = Dataset(directories[expt] + oce_file_end, 'r')
    node_temp_end = id.variables['temp'][0,:]
    id.close()
    node_temp_diff = node_temp_end - node_temp_beg
    i = 0
    for elm in elements:
        bwtemp_diff[expt,i] = mean([node_temp_diff[elm.nodes[0].find_bottom().id], node_temp_diff[elm.nodes[1].find_bottom().id], node_temp_diff[elm.nodes[2].find_bottom().id]])
        i += 1

print 'Calculating bottom water salinity'
bwsalt_beg = zeros(num_elm)
# Read full 3D field to start
id = Dataset(directory_beg + oce_file_beg, 'r')
node_salt_beg = id.variables['salt'][0,:]
id.close()
# For each element, find the bottom nodes and calculate average
i = 0
for elm in elements:
    bwsalt_beg[i] = mean([node_salt_beg[elm.nodes[0].find_bottom().id], node_salt_beg[elm.nodes[1].find_bottom().id], node_salt_beg[elm.nodes[2].find_bottom().id]])
    i += 1
bwsalt_diff = zeros([num_expts, num_elm])
for expt in range(num_expts):
    id = Dataset(directories[expt] + oce_file_end, 'r')
    node_salt_end = id.variables['salt'][0,:]
    id.close()
    node_salt_diff = node_salt_end - node_salt_beg
    i = 0
    for elm in elements:
        bwsalt_diff[expt,i] = mean([node_salt_diff[elm.nodes[0].find_bottom().id], node_salt_diff[elm.nodes[1].find_bottom().id], node_salt_diff[elm.nodes[2].find_bottom().id]])
        i += 1

print 'Plotting'
x_min_tmp = -17.5
x_max_tmp = -10.5
y_min_tmp = -11.25
y_max_tmp = -2.25
fig = figure(figsize=(15,10.5))
fig.patch.set_facecolor('white')
# Melt rate
gs_a = GridSpec(1,6)
gs_a.update(left=0.08, right=0.92, bottom=0.62, top=0.89, wspace=0.05)
cbaxes_left = fig.add_axes([0.02, 0.67, 0.02, 0.18])
ticks_left = arange(0, 6+2, 2)
cbaxes_right = fig.add_axes([0.93, 0.67, 0.02, 0.18])
ticks_right = arange(0, 18+6, 6)
plot_melt(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_a, cbaxes_left, ticks_left, cbaxes_right, ticks_right, [1, 2, 4], 'a', 1.15, set_diff_max=18, lon_lines=[-120, -104], lat_lines=[-75])
# Bottom water temperature
x_lonlat_reg, y_lonlat_reg, bathy_reg = interpolate_bathy(-140, -90, -76, -70)
gs_b = GridSpec(1,6)
gs_b.update(left=0.08, right=0.92, bottom=0.32, top=0.59, wspace=0.05)
cbaxes_left = fig.add_axes([0.02, 0.37, 0.02, 0.18])
ticks_left = arange(-1.6, 0.4+0.5, 0.5)
cbaxes_right = fig.add_axes([0.93, 0.37, 0.02, 0.18])
ticks_right = arange(0, 1.8+0.6, 0.6)
plot_bwtemp(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_b, cbaxes_left, ticks_left, cbaxes_right, ticks_right, 'b', bathy_contour=1500)
# Bottom water salinity
gs_c = GridSpec(1,6)
gs_c.update(left=0.08, right=0.92, bottom=0.02, top=0.29, wspace=0.05)
cbaxes_left = fig.add_axes([0.02, 0.07, 0.02, 0.18])
ticks_left = arange(34.1, 34.7+0.2, 0.2)
cbaxes_right = fig.add_axes([0.93, 0.07, 0.02, 0.18])
ticks_right = arange(-0.3, 0.15+0.15, 0.15)
plot_bwsalt(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_c, cbaxes_left, ticks_left, cbaxes_right, ticks_right, 'c', bathy_contour=1500)
suptitle('Amundsen Sea', fontsize=30)
fig.show()
fig.savefig('amundsen.png')


