from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from fesom_grid import *
from fesom_sidegrid import *
from unesco import *
from in_triangle import *

def amundsen_slices_before_after (rcp, model, save=False, fig_name=None):
    
    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/meshB/'
    file_beg = '/short/y99/kaa561/FESOM/highres_spinup/sept.oce.1996.2005.nc'
    file_end = '/short/y99/kaa561/FESOM/rcp'+rcp+'_'+model+'/sept.oce.2091.2100.nc'
    # Spatial bounds and labels
    lon0 = -104
    lon_string = str(int(round(-lon0)))+r'$^{\circ}$W'
    lat_min = -75.1
    lat_max = -70.8
    lat_ticks = arange(-75, -71+1, 1)
    lat_labels = []
    for lat in lat_ticks:
        lat_labels.append(str(int(round(-lat))) + r'$^{\circ}$S')
    depth_min = -1200
    depth_max = 0
    depth_ticks = arange(-1000, 0+250, 250)
    depth_labels = []
    for depth in depth_ticks:
        depth_labels.append(str(int(round(-depth))))
    if model == 'M':
        model_title = 'MMM'
    elif model == 'A':
        model_title = 'ACCESS'
    # Bounds on colourbars
    temp_min = -1.8
    temp_max = 1.1
    temp_ticks = arange(-1.5, 1+0.5, 0.5)
    salt_min = 33.6
    salt_max = 34.6
    salt_ticks = arange(33.6, 34.6+0.2, 0.2)
    # Density contours to plot
    density_contours = [27.45, 27.55]
    # Parameters for regular grid interpolation (needed for density contours)
    num_lat = 500
    num_depth = 250

    print 'Building FESOM mesh'
    elm2D = fesom_grid(mesh_path)
    print 'Reading temperature and salinity data'
    id = Dataset(file_beg, 'r')
    temp_nodes_beg = id.variables['temp'][0,:]
    salt_nodes_beg = id.variables['salt'][0,:]
    id.close()
    id = Dataset(file_end, 'r')
    temp_nodes_end = id.variables['temp'][0,:]
    salt_nodes_end = id.variables['salt'][0,:]

    print 'Calculating density'
    density_nodes_beg = unesco(temp_nodes_beg, salt_nodes_beg, zeros(shape(temp_nodes_beg))) - 1000
    density_nodes_end = unesco(temp_nodes_end, salt_nodes_end, zeros(shape(temp_nodes_end))) - 1000

    print 'Interpolating to ' + str(lon0)
    # Build arrays of SideElements making up zonal slices
    # Start with beginning
    selements_temp_beg = fesom_sidegrid(elm2D, temp_nodes_beg, lon0, lat_max)
    selements_salt_beg = fesom_sidegrid(elm2D, salt_nodes_beg, lon0, lat_max)
    # Build array of quadrilateral patches for the plots, and data values
    # corresponding to each SideElement
    patches = []
    temp_beg = []
    for selm in selements_temp_beg:
        # Make patch
        coord = transpose(vstack((selm.y, selm.z)))
        patches.append(Polygon(coord, True, linewidth=0.))
        # Save data value
        temp_beg.append(selm.var)
    temp_beg = array(temp_beg)
    print 'Temp bounds, beginning: ' + str(amin(temp_beg)) + ' ' + str(amax(temp_beg))
    # Salinity has same patches but different values
    salt_beg = []
    for selm in selements_salt_beg:
        salt_beg.append(selm.var)
    salt_beg = array(salt_beg)
    print 'Salt bounds, beginning: ' + str(amin(salt_beg)) + ' ' + str(amax(salt_beg))
    # Repeat for end
    selements_temp_end = fesom_sidegrid(elm2D, temp_nodes_end, lon0, lat_max)
    selements_salt_end = fesom_sidegrid(elm2D, salt_nodes_end, lon0, lat_max)
    temp_end = []
    for selm in selements_temp_end:
        temp_end.append(selm.var)
    temp_end = array(temp_end)
    print 'Temp bounds, end: ' + str(amin(temp_end)) + ' ' + str(amax(temp_end))
    salt_end = []
    for selm in selements_salt_end:
        salt_end.append(selm.var)
    salt_end = array(salt_end)
    print 'Salt bounds, end: ' + str(amin(salt_end)) + ' ' + str(amax(salt_end))

    print 'Interpolating density to regular grid'
    lat_reg = linspace(lat_min, lat_max, num_lat)
    depth_reg = linspace(-depth_max, -depth_min, num_depth)
    density_reg_beg = zeros([num_depth, num_lat])
    density_reg_end = zeros([num_depth, num_lat])
    density_reg_beg[:,:] = NaN
    density_reg_end[:,:] = NaN
    # For each element, check if a point on the regular grid lies
    # within. If so, do barycentric interpolation to that point, at each
    # depth on the regular grid.
    for elm in elm2D:
        # Check if this element crosses lon0
        if amin(elm.lon) < lon0 and amax(elm.lon) > lon0:
            # Check if we are within the latitude bounds
            if amax(elm.lat) > lat_min and amin(elm.lat) < lat_max:
                # Find largest regular latitude value south of Element
                tmp = nonzero(lat_reg > amin(elm.lat))[0]
                if len(tmp) == 0:
                    # Element crosses the southern boundary
                    jS = 0
                else:
                    jS = tmp[0] - 1
                # Find smallest regular latitude north of Element
                tmp = nonzero(lat_reg > amax(elm.lat))[0]
                if len(tmp) == 0:
                    # Element crosses the northern boundary
                    jN = num_lat
                else:
                    jN = tmp[0]
                for j in range(jS+1,jN):
                    # There is a chance that the regular gridpoint at j
                    # lies within this element
                    lat0 = lat_reg[j]
                    if in_triangle(elm, lon0, lat0):
                        # Yes it does
                        # Get area of entire triangle
                        area = triangle_area(elm.lon, elm.lat)
                        # Get area of each sub-triangle formed by (lon0, lat0)
                        area0 = triangle_area([lon0, elm.lon[1], elm.lon[2]], [lat0, elm.lat[1], elm.lat[2]])
                        area1 = triangle_area([lon0, elm.lon[0], elm.lon[2]], [lat0, elm.lat[0], elm.lat[2]])
                        area2 = triangle_area([lon0, elm.lon[0], elm.lon[1]], [lat0, elm.lat[0], elm.lat[1]])
                        # Find fractional area of each
                        cff = [area0/area, area1/area, area2/area]
                        # Interpolate each depth value
                        for k in range(num_depth):
                            # Linear interpolation in the vertical for the
                            # value at each corner of the triangle
                            node_vals_beg = []
                            node_vals_end = []
                            for n in range(3):
                                id1, id2, coeff1, coeff2 = elm.nodes[n].find_depth(depth_reg[k])
                                if any(isnan(array([id1, id2, coeff1, coeff2]))):
                                    # No ocean data here (seafloor or ice shelf)
                                    node_vals_beg.append(NaN)
                                    node_vals_end.append(NaN)
                                else:
                                    node_vals_beg.append(coeff1*density_nodes_beg[id1] + coeff2*density_nodes_beg[id2])
                                    node_vals_end.append(coeff1*density_nodes_end[id1] + coeff2*density_nodes_end[id2])
                            if any(isnan(node_vals_beg)):
                                pass
                            else:
                                # Barycentric interpolation for the value at
                                # lon0, lat0
                                density_reg_beg[k,j] = sum(array(cff)*array(node_vals_beg))
                                density_reg_end[k,j] = sum(array(cff)*array(node_vals_end))
    density_reg_beg = ma.masked_where(isnan(density_reg_beg), density_reg_beg)
    density_reg_end = ma.masked_where(isnan(density_reg_end), density_reg_end)
    depth_reg = -1*depth_reg

    print 'Plotting'
    fig = figure(figsize=(16,10))    
    # Temperature
    gs_temp = GridSpec(1,2)
    gs_temp.update(left=0.08, right=0.9, bottom=0.5, top=0.88, wspace=0.05, hspace=0.5)
    # Beginning
    ax = subplot(gs_temp[0,0])
    img = PatchCollection(patches, cmap='jet')
    img.set_array(temp_beg)
    img.set_edgecolor('face')
    img.set_clim(vmin=temp_min, vmax=temp_max)
    ax.add_collection(img)
    # Overlay density contours on regular grid
    contour(lat_reg, depth_reg, density_reg_beg, levels=density_contours, colors='black')
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    title(r'Temperature ($^{\circ}$C), 1996-2005', fontsize=24)
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels([])
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels, fontsize=16)
    ylabel('Depth (m)', fontsize=18)
    # End
    ax = subplot(gs_temp[0,1])
    img = PatchCollection(patches, cmap='jet')
    img.set_array(temp_end)
    img.set_edgecolor('face')
    img.set_clim(vmin=temp_min, vmax=temp_max)
    ax.add_collection(img)
    contour(lat_reg, depth_reg, density_reg_end, levels=density_contours, colors='black')
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    title(r'Temperature ($^{\circ}$C), 2091-2100', fontsize=24)
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels([])
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels([])
    # Add a colorbar on the right
    cbaxes = fig.add_axes([0.92, 0.55, 0.02, 0.28])
    cbar = colorbar(img, cax=cbaxes, extend='both', ticks=temp_ticks)
    cbar.ax.tick_params(labelsize=16)
    # Salinity
    gs_salt = GridSpec(1,2)
    gs_salt.update(left=0.08, right=0.9, bottom=0.07, top=0.45, wspace=0.05, hspace=0.5)
    # Beginning
    ax = subplot(gs_salt[0,0])
    img = PatchCollection(patches, cmap='jet')
    img.set_array(salt_beg)
    img.set_edgecolor('face')
    img.set_clim(vmin=salt_min, vmax=salt_max)
    ax.add_collection(img)
    contour(lat_reg, depth_reg, density_reg_beg, levels=density_contours, colors='black')
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    title('Salinity (psu), 1996-2005', fontsize=24)
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    xlabel('Latitude', fontsize=18)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels, fontsize=16)    
    ylabel('Depth (m)', fontsize=18)
    # End
    ax = subplot(gs_salt[0,1])
    img = PatchCollection(patches, cmap='jet')
    img.set_array(salt_end)
    img.set_edgecolor('face')
    img.set_clim(vmin=salt_min, vmax=salt_max)
    ax.add_collection(img)
    contour(lat_reg, depth_reg, density_reg_end, levels=density_contours, colors='black')
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    title('Salinity (psu), 2091-2100', fontsize=24)
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    xlabel('Latitude', fontsize=18)
    ax.set_yticks(lat_ticks)
    ax.set_yticklabels([])
    # Add a colorbar on the right
    cbaxes = fig.add_axes([0.92, 0.12, 0.02, 0.28])
    cbar = colorbar(img, cax=cbaxes, extend='both', ticks=salt_ticks)
    cbar.ax.tick_params(labelsize=16)
    # Main title
    suptitle('RCP ' + rcp[0] + '.' + rcp[1] + ' ' + model_title + ', ' + lon_string + ' (Amundsen Sea), September', fontsize=28)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    key = int(raw_input('RCP 4.5 (4) or 8.5 (8)? '))
    if key == 4:
        rcp = '45'
    elif key == 8:
        rcp = '85'
    key = int(raw_input('Multi-model mean (1) or ACCESS 1.0 (2)? '))
    if key == 1:
        model = 'M'
    elif key == 2:
        model = 'A'
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input('Filename for figure: ')
    elif action == 'd':
        save = False
        fig_name = None
    amundsen_slices_before_after (rcp, model, save, fig_name)
