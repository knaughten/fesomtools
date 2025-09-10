from netCDF4 import Dataset
from numpy import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.cm import *
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import *
from fesom_grid import *

def ross_plots ():

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/meshB/'
    forcing_file_beg = '/short/y99/kaa561/FESOM/highres_spinup/annual_avg.forcing.diag.1996.2005.nc'
    forcing_file_end = '/short/y99/kaa561/FESOM/rcp85_A/annual_avg.forcing.diag.2091.2100.nc'
    forcing_file_2094 = '/short/y99/kaa561/FESOM/rcp85_A/annual_avg.forcing.diag.2094.nc'
    oce_file_beg = '/short/y99/kaa561/FESOM/highres_spinup/annual_avg.oce.mean.1996.2005.nc'
    oce_file_end = '/short/y99/kaa561/FESOM/rcp85_A/annual_avg.oce.mean.2091.2100.nc'
    oce_file_2094 = '/short/y99/kaa561/FESOM/rcp85_A/annual_avg.oce.mean.2094.nc'
    oce2_file_beg = '/short/y99/kaa561/FESOM/highres_spinup/seasonal_climatology_oce_1996_2005.nc'
    oce2_file_end = '/short/y99/kaa561/FESOM/rcp85_A/seasonal_climatology_oce_2091_2100.nc'
    oce2_file_2094 = '/short/y99/kaa561/FESOM/rcp85_A/seasonal_climatology_oce_2094.nc'
    ice_file_beg = '/short/y99/kaa561/FESOM/highres_spinup/seasonal_climatology_ice_1996_2005.nc'
    ice_file_end = '/short/y99/kaa561/FESOM/rcp85_A/seasonal_climatology_ice_2091_2100.nc'
    ice_file_2094 = '/short/y99/kaa561/FESOM/rcp85_A/seasonal_climatology_ice_2094.nc'
    # Bounds on plot (in polar coordinate transformation)
    x_min = -5.5
    x_max = 4
    y_min = -13.8
    y_max = -4.75
    # Plotting parameters
    circumpolar = True
    # Season names for plot titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Seconds per year
    sec_per_year = 365.25*24*3600    

    print('Building mesh')
    elements = fesom_grid(mesh_path, circumpolar)
    # Build one set of plotting patches with all elements, one with
    # ice shelf cavities masked, and one with open ocean masked
    patches_all = []
    patches_ice = []
    patches_ocn = []
    for elm in elements:
        coord = transpose(vstack((elm.x, elm.y)))
        patches_all.append(Polygon(coord, True, linewidth=0.))
        if elm.cavity:
            patches_ice.append(Polygon(coord, True, linewidth=0.))
        else:
            patches_ocn.append(Polygon(coord, True, linewidth=0.))
    num_elm = len(patches_all)
    num_elm_ice = len(patches_ice)
    num_elm_ocn = len(patches_ocn)
    # Build ice shelf front contours
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
    # Set up a grey square to fill the background with land
    x_reg, y_reg = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
    land_square = zeros(shape(x_reg))

    print('Processing ice shelf melt rate')
    # Read annually averaged data, and convert from m/s to m/y
    id = Dataset(forcing_file_beg, 'r')
    wnet_nodes_beg = id.variables['wnet'][0,:]*sec_per_year
    id.close()
    id = Dataset(forcing_file_end, 'r')
    # Get difference from beginning
    wnet_nodes_end_diff = id.variables['wnet'][0,:]*sec_per_year - wnet_nodes_beg
    id.close()
    id = Dataset(forcing_file_2094, 'r')
    wnet_nodes_2094_diff = id.variables['wnet'][0,:]*sec_per_year - wnet_nodes_beg
    id.close()
    # Now average over each cavity element
    ismr_beg = []
    ismr_end_diff = []
    ismr_2094_diff = []
    for elm in elements:
        if elm.cavity:
            ismr_beg.append(mean([wnet_nodes_beg[elm.nodes[0].id], wnet_nodes_beg[elm.nodes[1].id], wnet_nodes_beg[elm.nodes[2].id]]))
            ismr_end_diff.append(mean([wnet_nodes_end_diff[elm.nodes[0].id], wnet_nodes_end_diff[elm.nodes[1].id], wnet_nodes_end_diff[elm.nodes[2].id]]))
            ismr_2094_diff.append(mean([wnet_nodes_2094_diff[elm.nodes[0].id], wnet_nodes_2094_diff[elm.nodes[1].id], wnet_nodes_2094_diff[elm.nodes[2].id]]))
    # Figure out bounds for colour scale
    # Min and max of beginning
    # Initialise with something impossible
    var_min = amax(array(ismr_beg))
    var_max = amin(array(ismr_beg))
    # Modify as needed
    i = 0
    for elm in elements:
        if elm.cavity:
            if any(elm.x >= x_min) and any(elm.x <= x_max) and any(elm.y >= y_min) and any(elm.y <= y_max):
                if ismr_beg[i] < var_min:
                    var_min = ismr_beg[i]
                if ismr_beg[i] > var_max:
                    var_max = ismr_beg[i]
            i += 1
    # Max absolute difference
    diff_max = 0
    i = 0
    for elm in elements:
        if elm.cavity:
            if any(elm.x >= x_min) and any(elm.x <= x_max) and any(elm.y >= y_min) and any(elm.y <= y_max):
                if abs(ismr_end_diff[i]) > diff_max:
                    diff_max = abs(ismr_end_diff[i])
                if abs(ismr_2094_diff[i]) > diff_max:
                    diff_max = abs(ismr_2094_diff[i])
            i += 1
    # Special colour map for absolute melt
    change_points = [0.5, 2, 3.5]
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
    # Plot
    fig = figure(figsize=(22,7))
    # 1996-2005
    ax = fig.add_subplot(1, 3, 1, aspect='equal')
    # Start with land background
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    # Add ice shelf elements
    img = PatchCollection(patches_ice, cmap=mf_cmap)
    img.set_array(array(ismr_beg))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    # Mask out the open ocean in white
    overlay = PatchCollection(patches_ocn, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    title('1996-2005', fontsize=20)
    # Colourbar on the left
    cbaxes = fig.add_axes([0.05, 0.25, 0.02, 0.5])
    cbar = colorbar(img, cax=cbaxes)
    # 2091-2100
    ax = fig.add_subplot(1, 3, 2, aspect='equal')
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches_ice, cmap='RdBu_r')
    img.set_array(array(ismr_end_diff))
    img.set_edgecolor('face')
    img.set_clim(vmin=-diff_max, vmax=diff_max)
    ax.add_collection(img)
    overlay = PatchCollection(patches_ocn, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    title('2091-2100 anomalies', fontsize=20)
    # 2094
    ax = fig.add_subplot(1, 3, 3, aspect='equal')
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches_ice, cmap='RdBu_r')
    img.set_array(array(ismr_2094_diff))
    img.set_edgecolor('face')
    img.set_clim(vmin=-diff_max, vmax=diff_max)
    ax.add_collection(img)
    overlay = PatchCollection(patches_ocn, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    title('2094 anomalies', fontsize=20)
    # Colourbar on the right
    cbaxes = fig.add_axes([0.92, 0.25, 0.02, 0.5])
    cbar = colorbar(img, cax=cbaxes)
    suptitle('Ice shelf melt rate (m/y)', fontsize=24)
    subplots_adjust(wspace=0.02, hspace=0.025)
    fig.show()
    fig.savefig('ross_melt.png')

    print('Processing bottom water temperature')
    # Read annually averaged data
    id = Dataset(oce_file_beg, 'r')
    temp_nodes_beg = id.variables['temp'][0,:]
    id.close()
    id = Dataset(oce_file_end, 'r')
    temp_nodes_end = id.variables['temp'][0,:]
    id.close()
    id = Dataset(oce_file_2094, 'r')
    temp_nodes_2094 = id.variables['temp'][0,:]
    id.close()
    # Now average bottom node temperatures over each element
    bwtemp_beg = []
    bwtemp_end = []
    bwtemp_2094 = []
    for elm in elements:
        bwtemp_beg.append(mean([temp_nodes_beg[elm.nodes[0].find_bottom().id], temp_nodes_beg[elm.nodes[1].find_bottom().id], temp_nodes_beg[elm.nodes[2].find_bottom().id]]))
        bwtemp_end.append(mean([temp_nodes_end[elm.nodes[0].find_bottom().id], temp_nodes_end[elm.nodes[1].find_bottom().id], temp_nodes_end[elm.nodes[2].find_bottom().id]]))
        bwtemp_2094.append(mean([temp_nodes_2094[elm.nodes[0].find_bottom().id], temp_nodes_2094[elm.nodes[1].find_bottom().id], temp_nodes_2094[elm.nodes[2].find_bottom().id]]))
    # Plot
    fig = figure(figsize=(22,7))
    # 1996-2005
    ax = fig.add_subplot(1, 3, 1, aspect='equal')
    # Start with land background
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    # Add all ocean elements
    img = PatchCollection(patches_all, cmap='jet')
    img.set_array(array(bwtemp_beg))
    img.set_edgecolor('face')
    img.set_clim(vmin=-2, vmax=-0.5)
    ax.add_collection(img)
    # Contour ice shelf fronts
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax.add_collection(contours)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    title('1996-2005', fontsize=20)
    # 2091-2100
    ax = fig.add_subplot(1, 3, 2, aspect='equal')
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches_all, cmap='jet')
    img.set_array(array(bwtemp_end))
    img.set_edgecolor('face')
    img.set_clim(vmin=-2, vmax=-0.5)
    ax.add_collection(img)
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax.add_collection(contours)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    title('2091-2100', fontsize=20)
    # 2094
    ax = fig.add_subplot(1, 3, 3, aspect='equal')
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches_all, cmap='jet')
    img.set_array(array(bwtemp_2094))
    img.set_edgecolor('face')
    img.set_clim(vmin=-2, vmax=-0.5)
    ax.add_collection(img)
    contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
    ax.add_collection(contours)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    title('2094', fontsize=20)
    # Horizontal colourbar below
    cbaxes = fig.add_axes([0.35, 0.04, 0.3, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='both')
    suptitle(r'Bottom water temperature ($^{\circ}$C)', fontsize=24)
    subplots_adjust(wspace=0.02, hspace=0.025)
    fig.show()
    fig.savefig('ross_bwtemp.png')

    print('Processing seasonal SSTs')
    # Read seasonally averaged data
    id = Dataset(oce2_file_beg, 'r')
    sst_nodes_beg = id.variables['temp'][:,:]
    id.close()
    id = Dataset(oce2_file_end, 'r')
    sst_nodes_end = id.variables['temp'][:,:]
    id.close()
    id = Dataset(oce2_file_2094, 'r')
    sst_nodes_2094 = id.variables['temp'][:,:]
    id.close()
    # Now average surface nodes over each non-cavity element
    sst_beg = empty([4, num_elm_ocn])
    sst_end = empty([4, num_elm_ocn])
    sst_2094 = empty([4, num_elm_ocn])
    i = 0
    for elm in elements:
        if not elm.cavity:
            sst_beg[:,i] = (sst_nodes_beg[:,elm.nodes[0].id] + sst_nodes_beg[:,elm.nodes[1].id] + sst_nodes_beg[:,elm.nodes[2].id])/3.0
            sst_end[:,i] = (sst_nodes_end[:,elm.nodes[0].id] + sst_nodes_end[:,elm.nodes[1].id] + sst_nodes_end[:,elm.nodes[2].id])/3.0
            sst_2094[:,i] = (sst_nodes_2094[:,elm.nodes[0].id] + sst_nodes_2094[:,elm.nodes[1].id] + sst_nodes_2094[:,elm.nodes[2].id])/3.0
            i += 1
    # Plot
    fig = figure(figsize=(19,11))
    for season in range(4):
        # 1996-2005
        ax = fig.add_subplot(3, 4, season+1, aspect='equal')
        # Start with land background
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        # Add open ocean elements
        img = PatchCollection(patches_ocn, cmap='jet')
        img.set_array(sst_beg[season,:])
        img.set_edgecolor('face')
        img.set_clim(vmin=-1.8, vmax=1.5)
        ax.add_collection(img)
        # Mask out cavities in white
        overlay = PatchCollection(patches_ice, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        ax.set_xticks([])
        ax.set_yticks([])
        title(season_names[season], fontsize=24)
        if season == 0:
            text(x_min-1, 0.5*(y_min+y_max), '1996-2005', fontsize=20, ha='center', rotation=90)
        # 2091-2100
        ax = fig.add_subplot(3, 4, season+5, aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(patches_ocn, cmap='jet')
        img.set_array(sst_end[season,:])
        img.set_edgecolor('face')
        img.set_clim(vmin=-1.8, vmax=1.5)
        ax.add_collection(img)
        overlay = PatchCollection(patches_ice, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        ax.set_xticks([])
        ax.set_yticks([])
        if season == 0:
            text(x_min-1, 0.5*(y_min+y_max), '2091-2100', fontsize=20, ha='center', rotation=90)
        # 2094
        ax = fig.add_subplot(3, 4, season+9, aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(patches_ocn, cmap='jet')
        img.set_array(sst_2094[season,:])
        img.set_edgecolor('face')
        img.set_clim(vmin=-1.8, vmax=1.5)
        ax.add_collection(img)
        overlay = PatchCollection(patches_ice, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        ax.set_xticks([])
        ax.set_yticks([])
        if season == 0:
            text(x_min-1, 0.5*(y_min+y_max), '2094', fontsize=20, ha='center', rotation=90)
        if season == 3:
            # Colourbar below
            cbaxes = fig.add_axes([0.35, 0.04, 0.3, 0.02])
            cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='both')
    suptitle(r'Sea surface temperature ($^{\circ}$C)', fontsize=24)
    subplots_adjust(wspace=0.025, hspace=0.025)
    fig.show()
    fig.savefig('ross_sst.png')

    print('Processing seasonal sea ice concentration')
    # Read seasonally averaged data
    id = Dataset(ice_file_beg, 'r')
    aice_nodes_beg = id.variables['area'][:,:]
    id.close()
    id = Dataset(ice_file_end, 'r')
    aice_nodes_end = id.variables['area'][:,:]
    id.close()
    id = Dataset(ice_file_2094, 'r')
    aice_nodes_2094 = id.variables['area'][:,:]
    id.close()
    # Now average nodes over each non-cavity element
    aice_beg = empty([4, num_elm_ocn])
    aice_end = empty([4, num_elm_ocn])
    aice_2094 = empty([4, num_elm_ocn])
    i = 0
    for elm in elements:
        if not elm.cavity:
            aice_beg[:,i] = (aice_nodes_beg[:,elm.nodes[0].id] + aice_nodes_beg[:,elm.nodes[1].id] + aice_nodes_beg[:,elm.nodes[2].id])/3.0
            aice_end[:,i] = (aice_nodes_end[:,elm.nodes[0].id] + aice_nodes_end[:,elm.nodes[1].id] + aice_nodes_end[:,elm.nodes[2].id])/3.0
            aice_2094[:,i] = (aice_nodes_2094[:,elm.nodes[0].id] + aice_nodes_2094[:,elm.nodes[1].id] + aice_nodes_2094[:,elm.nodes[2].id])/3.0
            i += 1
    # Plot
    fig = figure(figsize=(19,11))
    for season in range(4):
        # 1996-2005
        ax = fig.add_subplot(3, 4, season+1, aspect='equal')
        # Start with land background
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        # Add open ocean elements
        img = PatchCollection(patches_ocn, cmap='jet')
        img.set_array(aice_beg[season,:])
        img.set_edgecolor('face')
        img.set_clim(vmin=0, vmax=1)
        ax.add_collection(img)
        # Mask out cavities in white
        overlay = PatchCollection(patches_ice, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        ax.set_xticks([])
        ax.set_yticks([])
        title(season_names[season], fontsize=24)
        if season == 0:
            text(x_min-1, 0.5*(y_min+y_max), '1996-2005', fontsize=20, ha='left', rotation=90)
        # 2091-2100
        ax = fig.add_subplot(3, 4, season+5, aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(patches_ocn, cmap='jet')
        img.set_array(aice_end[season,:])
        img.set_edgecolor('face')
        img.set_clim(vmin=0, vmax=1)
        ax.add_collection(img)
        overlay = PatchCollection(patches_ice, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        ax.set_xticks([])
        ax.set_yticks([])
        if season == 0:
            text(x_min-1, 0.5*(y_min+y_max), '2091-2100', fontsize=20, ha='left', rotation=90)
        # 2094
        ax = fig.add_subplot(3, 4, season+9, aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(patches_ocn, cmap='jet')
        img.set_array(aice_2094[season,:])
        img.set_edgecolor('face')
        img.set_clim(vmin=0, vmax=1)
        ax.add_collection(img)
        overlay = PatchCollection(patches_ice, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        ax.set_xticks([])
        ax.set_yticks([])
        if season == 0:
            text(x_min-1, 0.5*(y_min+y_max), '2094', fontsize=20, ha='left', rotation=90)
        if season == 3:
            # Colourbar below
            cbaxes = fig.add_axes([0.35, 0.04, 0.3, 0.02])
            cbar = colorbar(img, orientation='horizontal', cax=cbaxes)
    suptitle('Sea ice concentration', fontsize=24)
    subplots_adjust(wspace=0.025, hspace=0.025)
    fig.show()
    fig.savefig('ross_aice.png')


# Command-line interface
if __name__ == "__main__":

    ross_plots()
