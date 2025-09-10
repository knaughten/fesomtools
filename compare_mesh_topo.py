from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.colors import ListedColormap
from matplotlib.patches import Polygon
from fesom_grid import *

def compare_mesh_topo ():

    # File paths
    # I know the copy-paste of 6 meshes is terrible programing practice.
    # But I am in a hurry.
    mesh_path_A = '../FESOM_mesh/low_res_fixsmooth/output/'
    mesh_path_B = '../FESOM_mesh/high_res_fixsmooth/output/'
    mesh_path_C = '../FESOM_mesh/extra_high_res_less_amery/output/'
    mesh_path_D = '../FESOM_mesh/extra_high_res_less_amery_more_acc/output/'
    mesh_path_E = '../FESOM_mesh/extra_high_res/output/'
    mesh_path_F = '../FESOM_mesh/extra_high_res_more_acc/output/'
    rtopo_data_file = '../FESOM_mesh/RTopo105/RTopo105_data.nc'
    rtopo_aux_file = '../FESOM_mesh/RTopo105/RTopo105_aux.nc'
    # Maximum j-index to read from RTopo (can skip most of the world)
    rtopo_max_j = 2000
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Name of each region
    region_names = ['Filchner-Ronne Ice Shelf', 'Eastern Weddell Region', 'Amery Ice Shelf', 'Australian Sector', 'Ross Sea', 'Amundsen Sea', 'Bellingshausen Sea', 'Larsen Ice Shelves']
    num_regions = len(region_names)
    # Beginning of filenames for figures
    fig_heads = ['filchner_ronne', 'eweddell', 'amery', 'australian', 'ross', 'amundsen', 'bellingshausen', 'larsen']
    # Bounds for each region (using polar coordinate transformation as below)
    x_min = [-14, -8, 15.25, 12, -9.5, -15.5, -20.25, -22.5]
    x_max = [-4.5, 13, 20.5, 25.5, 4, -10.5, -15.5, -14.5]
    y_min = [1, 12, 4.75, -20, -13, -11.25, -4.5, 8.3]
    y_max = [10, 21, 8, 4, -4.75, -2.25, 7.6, 13]
    # Size of each plot in the y direction
    ysize = [12, 9, 10, 13, 10, 13, 15, 10]
    # Grey colourmap for RTopo land mask
    grey_cmap = ListedColormap([(0.6, 0.6, 0.6)])

    print('Reading mesh A')
    # Build mesh
    elements_A = fesom_grid(mesh_path_A, circumpolar=True, cross_180=True)
    # Build patches for ocean, open ocean, and ice shelves
    patches_A = []
    shelf_patches_A = []
    ocean_patches_A = []
    for elm in elements_A:
        coord = transpose(vstack((elm.x, elm.y)))
        if elm.cavity:
            shelf_patches_A.append(Polygon(coord, True, linewidth=0.))
        else:
            ocean_patches_A.append(Polygon(coord, True, linewidth=0.))
        patches_A.append(Polygon(coord, True, linewidth=0.))
    # Calculate bathy, draft, wct at each element
    bathy_A = []
    draft_A = []
    wct_A = []
    for elm in elements_A:
        # Bathymetry is depth of bottom layer, averaged over 3 Nodes
        bathy_A.append(mean([elm.nodes[0].find_bottom().depth, elm.nodes[1].find_bottom().depth, elm.nodes[2].find_bottom().depth]))
        if elm.cavity:
            # Ice shelf draft is depth of surface layer
            draft_A.append(mean([elm.nodes[0].depth, elm.nodes[1].depth, elm.nodes[2].depth]))
            # Water column thickness is depth of bottom layer minus depth of
            # surface layer
            wct_A.append(mean([elm.nodes[0].find_bottom().depth - elm.nodes[0].depth, elm.nodes[1].find_bottom().depth - elm.nodes[1].depth, elm.nodes[2].find_bottom().depth - elm.nodes[2].depth]))

    print('Reading mesh B')
    elements_B = fesom_grid(mesh_path_B, circumpolar=True, cross_180=True)
    patches_B = []
    shelf_patches_B = []
    ocean_patches_B = []
    for elm in elements_B:
        coord = transpose(vstack((elm.x, elm.y)))
        if elm.cavity:
            shelf_patches_B.append(Polygon(coord, True, linewidth=0.))
        else:
            ocean_patches_B.append(Polygon(coord, True, linewidth=0.))
        patches_B.append(Polygon(coord, True, linewidth=0.))
    bathy_B = []
    draft_B = []
    wct_B = []
    for elm in elements_B:
        bathy_B.append(mean([elm.nodes[0].find_bottom().depth, elm.nodes[1].find_bottom().depth, elm.nodes[2].find_bottom().depth]))
        if elm.cavity:
            draft_B.append(mean([elm.nodes[0].depth, elm.nodes[1].depth, elm.nodes[2].depth]))
            wct_B.append(mean([elm.nodes[0].find_bottom().depth - elm.nodes[0].depth, elm.nodes[1].find_bottom().depth - elm.nodes[1].depth, elm.nodes[2].find_bottom().depth - elm.nodes[2].depth]))

    print('Reading mesh C')
    elements_C = fesom_grid(mesh_path_C, circumpolar=True, cross_180=True)
    patches_C = []
    shelf_patches_C = []
    ocean_patches_C = []
    for elm in elements_C:
        coord = transpose(vstack((elm.x, elm.y)))
        if elm.cavity:
            shelf_patches_C.append(Polygon(coord, True, linewidth=0.))
        else:
            ocean_patches_C.append(Polygon(coord, True, linewidth=0.))
        patches_C.append(Polygon(coord, True, linewidth=0.))
    bathy_C = []
    draft_C = []
    wct_C = []
    for elm in elements_C:
        bathy_C.append(mean([elm.nodes[0].find_bottom().depth, elm.nodes[1].find_bottom().depth, elm.nodes[2].find_bottom().depth]))
        if elm.cavity:
            draft_C.append(mean([elm.nodes[0].depth, elm.nodes[1].depth, elm.nodes[2].depth]))
            wct_C.append(mean([elm.nodes[0].find_bottom().depth - elm.nodes[0].depth, elm.nodes[1].find_bottom().depth - elm.nodes[1].depth, elm.nodes[2].find_bottom().depth - elm.nodes[2].depth]))

    print('Reading mesh D')
    elements_D = fesom_grid(mesh_path_D, circumpolar=True, cross_180=True)
    patches_D = []
    shelf_patches_D = []
    ocean_patches_D = []
    for elm in elements_D:
        coord = transpose(vstack((elm.x, elm.y)))
        if elm.cavity:
            shelf_patches_D.append(Polygon(coord, True, linewidth=0.))
        else:
            ocean_patches_D.append(Polygon(coord, True, linewidth=0.))
        patches_D.append(Polygon(coord, True, linewidth=0.))
    bathy_D = []
    draft_D = []
    wct_D = []
    for elm in elements_D:
        bathy_D.append(mean([elm.nodes[0].find_bottom().depth, elm.nodes[1].find_bottom().depth, elm.nodes[2].find_bottom().depth]))
        if elm.cavity:
            draft_D.append(mean([elm.nodes[0].depth, elm.nodes[1].depth, elm.nodes[2].depth]))
            wct_D.append(mean([elm.nodes[0].find_bottom().depth - elm.nodes[0].depth, elm.nodes[1].find_bottom().depth - elm.nodes[1].depth, elm.nodes[2].find_bottom().depth - elm.nodes[2].depth]))

    print('Reading mesh E')
    elements_E = fesom_grid(mesh_path_E, circumpolar=True, cross_180=True)
    patches_E = []
    shelf_patches_E = []
    ocean_patches_E = []
    for elm in elements_E:
        coord = transpose(vstack((elm.x, elm.y)))
        if elm.cavity:
            shelf_patches_E.append(Polygon(coord, True, linewidth=0.))
        else:
            ocean_patches_E.append(Polygon(coord, True, linewidth=0.))
        patches_E.append(Polygon(coord, True, linewidth=0.))
    bathy_E = []
    draft_E = []
    wct_E = []
    for elm in elements_E:
        bathy_E.append(mean([elm.nodes[0].find_bottom().depth, elm.nodes[1].find_bottom().depth, elm.nodes[2].find_bottom().depth]))
        if elm.cavity:
            draft_E.append(mean([elm.nodes[0].depth, elm.nodes[1].depth, elm.nodes[2].depth]))
            wct_E.append(mean([elm.nodes[0].find_bottom().depth - elm.nodes[0].depth, elm.nodes[1].find_bottom().depth - elm.nodes[1].depth, elm.nodes[2].find_bottom().depth - elm.nodes[2].depth]))

    print('Reading mesh F')
    elements_F = fesom_grid(mesh_path_F, circumpolar=True, cross_180=True)
    patches_F = []
    shelf_patches_F = []
    ocean_patches_F = []
    for elm in elements_F:
        coord = transpose(vstack((elm.x, elm.y)))
        if elm.cavity:
            shelf_patches_F.append(Polygon(coord, True, linewidth=0.))
        else:
            ocean_patches_F.append(Polygon(coord, True, linewidth=0.))
        patches_F.append(Polygon(coord, True, linewidth=0.))
    bathy_F = []
    draft_F = []
    wct_F = []
    for elm in elements_F:
        bathy_F.append(mean([elm.nodes[0].find_bottom().depth, elm.nodes[1].find_bottom().depth, elm.nodes[2].find_bottom().depth]))
        if elm.cavity:
            draft_F.append(mean([elm.nodes[0].depth, elm.nodes[1].depth, elm.nodes[2].depth]))
            wct_F.append(mean([elm.nodes[0].find_bottom().depth - elm.nodes[0].depth, elm.nodes[1].find_bottom().depth - elm.nodes[1].depth, elm.nodes[2].find_bottom().depth - elm.nodes[2].depth]))

    print('Reading RTopo')
    # Read grid, bathy, draft
    id = Dataset(rtopo_data_file, 'r')
    rtopo_lon = id.variables['lon'][:]
    rtopo_lat = id.variables['lat'][:rtopo_max_j]
    rtopo_bathy = -1*id.variables['bathy'][:rtopo_max_j,:]
    rtopo_draft = -1*id.variables['draft'][:rtopo_max_j,:]
    id.close()
    # Calculate wct
    rtopo_wct = rtopo_bathy - rtopo_draft
    # Read composite mask
    id = Dataset(rtopo_aux_file, 'r')
    rtopo_amask = id.variables['amask'][:rtopo_max_j,:]
    id.close()
    # Create ocean and ice shelf masks
    rtopo_omask = zeros(shape(rtopo_amask))
    rtopo_imask = zeros(shape(rtopo_amask))
    index = rtopo_amask == 0
    rtopo_omask[index] = 1
    index = rtopo_amask == 2
    rtopo_omask[index] = 1    
    rtopo_imask[index] = 1
    # Create land mask for plotting in grey
    rtopo_land = ma.masked_where(rtopo_amask==0, rtopo_amask)
    # Apply masks to each variable
    rtopo_bathy = ma.masked_where(rtopo_omask==0, rtopo_bathy)
    rtopo_draft = ma.masked_where(rtopo_imask==0, rtopo_draft)
    rtopo_wct = ma.masked_where(rtopo_imask==0, rtopo_wct)
    # Calculate polar coordinate transformation
    rtopo_lon2d, rtopo_lat2d = meshgrid(rtopo_lon, rtopo_lat)
    rtopo_x = -(rtopo_lat2d+90)*cos(rtopo_lon2d*deg2rad+pi/2)
    rtopo_y = (rtopo_lat2d+90)*sin(rtopo_lon2d*deg2rad+pi/2)

    # Loop over regions
    for index in range(num_regions):
        print('Processing ' + region_names[index])
        # Set up a grey square to fill the background with land
        x_reg, y_reg = meshgrid(linspace(x_min[index], x_max[index], num=100), linspace(y_min[index], y_max[index], num=100))
        land_square = zeros(shape(x_reg))
        # Find bounds on variables in this region
        # Start with RTopo
        loc = (rtopo_x >= x_min[index])*(rtopo_x <= x_max[index])*(rtopo_y >= y_min[index])*(rtopo_y <= y_max[index])
        bathy_min = amin(rtopo_bathy[loc])
        bathy_max = amax(rtopo_bathy[loc])
        draft_min = amin(rtopo_draft[loc])
        draft_max = amax(rtopo_draft[loc])
        wct_min = amin(rtopo_wct[loc])
        wct_max = amax(rtopo_wct[loc])
        # Modify with each mesh
        i = 0
        j = 0
        for elm in elements_A:
            if any(elm.x >= x_min[index]) and any(elm.x <= x_max[index]) and any(elm.y >= y_min[index]) and any(elm.y <= y_max[index]):
                if bathy_A[i] < bathy_min:
                    bathy_min = bathy_A[i]
                if bathy_A[i] > bathy_max:
                    bathy_max = bathy_A[i]
                if elm.cavity:
                    if draft_A[j] < draft_min:
                        draft_min = draft_A[j]
                    if draft_A[j] > draft_max:
                        draft_max = draft_A[j]
                    if wct_A[j] < wct_min:
                        wct_min = wct_A[j]
                    if wct_A[j] > wct_max:
                        wct_max = wct_A[j]
            i += 1
            if elm.cavity:
                j += 1
        i = 0
        j = 0
        for elm in elements_B:
            if any(elm.x >= x_min[index]) and any(elm.x <= x_max[index]) and any(elm.y >= y_min[index]) and any(elm.y <= y_max[index]):
                if bathy_B[i] < bathy_min:
                    bathy_min = bathy_B[i]
                if bathy_B[i] > bathy_max:
                    bathy_max = bathy_B[i]
                if elm.cavity:
                    if draft_B[j] < draft_min:
                        draft_min = draft_B[j]
                    if draft_B[j] > draft_max:
                        draft_max = draft_B[j]
                    if wct_B[j] < wct_min:
                        wct_min = wct_B[j]
                    if wct_B[j] > wct_max:
                        wct_max = wct_B[j]
            i += 1
            if elm.cavity:
                j += 1
        i = 0
        j = 0
        for elm in elements_C:
            if any(elm.x >= x_min[index]) and any(elm.x <= x_max[index]) and any(elm.y >= y_min[index]) and any(elm.y <= y_max[index]):
                if bathy_C[i] < bathy_min:
                    bathy_min = bathy_C[i]
                if bathy_C[i] > bathy_max:
                    bathy_max = bathy_C[i]
                if elm.cavity:
                    if draft_C[j] < draft_min:
                        draft_min = draft_C[j]
                    if draft_C[j] > draft_max:
                        draft_max = draft_C[j]
                    if wct_C[j] < wct_min:
                        wct_min = wct_C[j]
                    if wct_C[j] > wct_max:
                        wct_max = wct_C[j]
            i += 1
            if elm.cavity:
                j += 1
        i = 0
        j = 0
        for elm in elements_D:
            if any(elm.x >= x_min[index]) and any(elm.x <= x_max[index]) and any(elm.y >= y_min[index]) and any(elm.y <= y_max[index]):
                if bathy_D[i] < bathy_min:
                    bathy_min = bathy_D[i]
                if bathy_D[i] > bathy_max:
                    bathy_max = bathy_D[i]
                if elm.cavity:
                    if draft_D[j] < draft_min:
                        draft_min = draft_D[j]
                    if draft_D[j] > draft_max:
                        draft_max = draft_D[j]
                    if wct_D[j] < wct_min:
                        wct_min = wct_D[j]
                    if wct_D[j] > wct_max:
                        wct_max = wct_D[j]
            i += 1
            if elm.cavity:
                j += 1
        i = 0
        j = 0
        for elm in elements_E:
            if any(elm.x >= x_min[index]) and any(elm.x <= x_max[index]) and any(elm.y >= y_min[index]) and any(elm.y <= y_max[index]):
                if bathy_E[i] < bathy_min:
                    bathy_min = bathy_E[i]
                if bathy_E[i] > bathy_max:
                    bathy_max = bathy_E[i]
                if elm.cavity:
                    if draft_E[j] < draft_min:
                        draft_min = draft_E[j]
                    if draft_E[j] > draft_max:
                        draft_max = draft_E[j]
                    if wct_E[j] < wct_min:
                        wct_min = wct_E[j]
                    if wct_E[j] > wct_max:
                        wct_max = wct_E[j]
            i += 1
            if elm.cavity:
                j += 1
        i = 0
        j = 0
        for elm in elements_F:
            if any(elm.x >= x_min[index]) and any(elm.x <= x_max[index]) and any(elm.y >= y_min[index]) and any(elm.y <= y_max[index]):
                if bathy_F[i] < bathy_min:
                    bathy_min = bathy_F[i]
                if bathy_F[i] > bathy_max:
                    bathy_max = bathy_F[i]
                if elm.cavity:
                    if draft_F[j] < draft_min:
                        draft_min = draft_F[j]
                    if draft_F[j] > draft_max:
                        draft_max = draft_F[j]
                    if wct_F[j] < wct_min:
                        wct_min = wct_F[j]
                    if wct_F[j] > wct_max:
                        wct_max = wct_F[j]
            i += 1
            if elm.cavity:
                j += 1
        # Plot
        fig = figure(figsize=(20, ysize[index]))
        fig.patch.set_facecolor('white')
        gs = GridSpec([3,7])
        gs.update(left=0.05, right=0.9, bottom=0.05, top=0.9, wspace=0.05, hspace=0.1)

        # RTopo
        # Bathymetry
        ax = subplot(gs[0,0], aspect='equal')
        # First plot land in grey
        pcolor(rtopo_x, rtopo_y, rtopo_land, cmap=grey_cmap)
        # Now plot bathymetry
        pcolor(rtopo_x, rtopo_y, rtopo_bathy, vmin=bathy_min, vmax=bathy_max)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        # Ice shelf draft
        ax = subplot(gs[1,0], aspect='equal')
        pcolor(rtopo_x, rtopo_y, rtopo_land, cmap=grey_cmap)
        pcolor(rtopo_x, rtopo_y, rtopo_draft, vmin=draft_min, vmax=draft_max)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        # Water column thickness
        ax = subplot(gs[2,0], aspect='equal')
        pcolor(rtopo_x, rtopo_y, rtopo_land, cmap=grey_cmap)
        pcolor(rtopo_x, rtopo_y, rtopo_wct, vmin=wct_min, vmax=wct_max)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        # Add RTopo label on bottom
        xlabel('RTopo')

        # Mesh A
        # Bathymetry
        ax = subplot(gs[0,1], aspect='equal')
        # Start with land background
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        # Add bathymetry data
        img = PatchCollection(patches_A, cmap='jet')
        img.set_array(array(bathy_A))
        img.set_edgecolor('face')
        img.set_clim(vmin=bathy_min, vmax=bathy_max)
        ax.add_collection(img)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        # Ice shelf draft
        ax = subplot(gs[1,1], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(shelf_patches_A, cmap='jet')
        img.set_array(array(draft_A))
        img.set_edgecolor('face')
        img.set_clim(vmin=draft_min, vmax=draft_max)
        ax.add_collection(img)
        # Mask out the open ocean in white
        overlay = PatchCollection(ocean_patches_A, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        # Water column thickness in cavity
        ax = subplot(gs[2,1], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(shelf_patches_A, cmap='jet')
        img.set_array(array(wct_A))
        img.set_edgecolor('face')
        img.set_clim(vmin=wct_min, vmax=wct_max)
        ax.add_collection(img)
        overlay = PatchCollection(ocean_patches_A, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        # Mesh title
        xlabel('Mesh A')

        # Mesh B
        # Bathymetry
        ax = subplot(gs[0,2], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(patches_B, cmap='jet')
        img.set_array(array(bathy_B))
        img.set_edgecolor('face')
        img.set_clim(vmin=bathy_min, vmax=bathy_max)
        ax.add_collection(img)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        # Ice shelf draft
        ax = subplot(gs[1,2], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(shelf_patches_B, cmap='jet')
        img.set_array(array(draft_B))
        img.set_edgecolor('face')
        img.set_clim(vmin=draft_min, vmax=draft_max)
        ax.add_collection(img)
        overlay = PatchCollection(ocean_patches_B, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        # Water column thickness in cavity
        ax = subplot(gs[2,2], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(shelf_patches_B, cmap='jet')
        img.set_array(array(wct_B))
        img.set_edgecolor('face')
        img.set_clim(vmin=wct_min, vmax=wct_max)
        ax.add_collection(img)
        overlay = PatchCollection(ocean_patches_B, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        xlabel('Mesh B')

        # Mesh C
        # Bathymetry
        ax = subplot(gs[0,3], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(patches_C, cmap='jet')
        img.set_array(array(bathy_C))
        img.set_edgecolor('face')
        img.set_clim(vmin=bathy_min, vmax=bathy_max)
        ax.add_collection(img)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        # Variable title
        title('Bathymetry (m)')
        # Ice shelf draft
        ax = subplot(gs[1,3], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(shelf_patches_C, cmap='jet')
        img.set_array(array(draft_C))
        img.set_edgecolor('face')
        img.set_clim(vmin=draft_min, vmax=draft_max)
        ax.add_collection(img)
        overlay = PatchCollection(ocean_patches_C, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        title('Ice shelf draft (m)')
        # Water column thickness in cavity
        ax = subplot(gs[2,3], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(shelf_patches_C, cmap='jet')
        img.set_array(array(wct_C))
        img.set_edgecolor('face')
        img.set_clim(vmin=wct_min, vmax=wct_max)
        ax.add_collection(img)
        overlay = PatchCollection(ocean_patches_C, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        title('Water column thickness (m)')
        xlabel('Mesh C')

        # Mesh D
        # Bathymetry
        ax = subplot(gs[0,4], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(patches_D, cmap='jet')
        img.set_array(array(bathy_D))
        img.set_edgecolor('face')
        img.set_clim(vmin=bathy_min, vmax=bathy_max)
        ax.add_collection(img)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        # Ice shelf draft
        ax = subplot(gs[1,4], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(shelf_patches_D, cmap='jet')
        img.set_array(array(draft_D))
        img.set_edgecolor('face')
        img.set_clim(vmin=draft_min, vmax=draft_max)
        ax.add_collection(img)
        overlay = PatchCollection(ocean_patches_D, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        # Water column thickness in cavity
        ax = subplot(gs[2,4], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(shelf_patches_D, cmap='jet')
        img.set_array(array(wct_D))
        img.set_edgecolor('face')
        img.set_clim(vmin=wct_min, vmax=wct_max)
        ax.add_collection(img)
        overlay = PatchCollection(ocean_patches_D, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        xlabel('Mesh D')

        # Mesh E
        # Bathymetry
        ax = subplot(gs[0,5], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(patches_E, cmap='jet')
        img.set_array(array(bathy_E))
        img.set_edgecolor('face')
        img.set_clim(vmin=bathy_min, vmax=bathy_max)
        ax.add_collection(img)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        # Ice shelf draft
        ax = subplot(gs[1,5], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(shelf_patches_E, cmap='jet')
        img.set_array(array(draft_E))
        img.set_edgecolor('face')
        img.set_clim(vmin=draft_min, vmax=draft_max)
        ax.add_collection(img)
        overlay = PatchCollection(ocean_patches_E, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        # Water column thickness in cavity
        ax = subplot(gs[2,5], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(shelf_patches_E, cmap='jet')
        img.set_array(array(wct_E))
        img.set_edgecolor('face')
        img.set_clim(vmin=wct_min, vmax=wct_max)
        ax.add_collection(img)
        overlay = PatchCollection(ocean_patches_E, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        xlabel('Mesh E')

        # Mesh F
        # Bathymetry
        ax = subplot(gs[0,6], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(patches_F, cmap='jet')
        img.set_array(array(bathy_F))
        img.set_edgecolor('face')
        img.set_clim(vmin=bathy_min, vmax=bathy_max)
        ax.add_collection(img)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        # Colourbar on right
        cbaxes = fig.add_axes([0.92, 0.7, 0.01, 0.15])
        cbar = colorbar(img, cax=cbaxes)
        # Ice shelf draft
        ax = subplot(gs[1,6], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(shelf_patches_F, cmap='jet')
        img.set_array(array(draft_F))
        img.set_edgecolor('face')
        img.set_clim(vmin=draft_min, vmax=draft_max)
        ax.add_collection(img)
        overlay = PatchCollection(ocean_patches_F, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        cbaxes = fig.add_axes([0.92, 0.4, 0.01, 0.15])
        cbar = colorbar(img, cax=cbaxes)
        # Water column thickness in cavity
        ax = subplot(gs[2,6], aspect='equal')
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(shelf_patches_F, cmap='jet')
        img.set_array(array(wct_F))
        img.set_edgecolor('face')
        img.set_clim(vmin=wct_min, vmax=wct_max)
        ax.add_collection(img)
        overlay = PatchCollection(ocean_patches_F, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        ax.set_xticks([])
        ax.set_yticks([])
        xlabel('Mesh F')
        cbaxes = fig.add_axes([0.92, 0.1, 0.01, 0.15])
        cbar = colorbar(img, cax=cbaxes)

        # Main title
        suptitle(region_names[index])
        fig.show()
        #fig.savefig(fig_heads[index] + '_compare_topo.png')


# Command-line interface
if __name__ == "__main__":

    compare_mesh_topo ()
                

    
    
    
    

    
