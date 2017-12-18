from numpy import *
from netCDF4 import Dataset
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.colors import *
from matplotlib.pyplot import *
from fesom_grid import *
from patches import *

def rcp_map_circles (key=1):

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/meshB/'
    directory_beg = '/short/y99/kaa561/FESOM/highres_spinup/'
    directories = ['/short/y99/kaa561/FESOM/rcp45_M/', '/short/y99/kaa561/FESOM/rcp45_A/', '/short/y99/kaa561/FESOM/rcp85_M/', '/short/y99/kaa561/FESOM/rcp85_A/']
    file_beg = 'annual_avg.forcing.diag.1996.2005.nc'
    file_end = 'annual_avg.forcing.diag.2091.2100.nc'
    # Titles for plotting
    expt_names = ['RCP 4.5 MMM', 'RCP 4.5 ACCESS', 'RCP 8.5 MMM', 'RCP 8.5 ACCESS']
    letters = ['a)', 'b)', 'c)', 'd)']
    num_expts = len(directories)
    # Northern boundary for plot: 60S
    nbdry = -60 + 90
    # Seconds per year
    sec_per_year = 365.25*24*3600
    # Density of ice in kg/m^3
    rho_ice = 916
    # Sectors to split Antarctica into
    sector_names = ['Filchner-Ronne Ice Shelf', 'Eastern Weddell Region', 'Amery Ice Shelf', 'Australian Sector', 'Ross Sea', 'Amundsen Sea', 'Bellingshausen Sea', 'Larsen Ice Shelves']
    # Number of sectors
    num_sectors = len(sector_names)
    # Colours for sectors
    sector_colours = [(0.6, 0.8, 1), (1, 1, 0), (1, 0.76, 0.4), (0, 0.9, 0.9), (0.73, 0.6, 1), (1, 0.4, 0.4), (0.52, 0.88, 0.52), (1, 0.8, 0.9)]
    # x and y coordinates for circles (using polar coordinate transformation)
    if key == 1:
        x_circles = [-12, 10, 25, 22, -2, -17, -23, -21]
        y_circles = [14, 24, 8, -16, -17, -13, 2, 18]
    elif key == 2:
        x_circles = [-12, 10, 25, 22, -2, -19, -24, -22]
        y_circles = [14, 24, 8, -17, -16, -15, 2, 18]
    # Scaling factor for radius of circles
    k = 0.4

    print 'Building mesh'
    # One set of elements for calculating: doesn't cross 180
    elements_calc = fesom_grid(mesh_path, circumpolar=True, cross_180=False)
    # Count the number of ice shelf elements
    num_cavity_elm = 0
    for elm in elements_calc:
        if elm.cavity:
            num_cavity_elm += 1
    # One set of elements and patches for plotting: crosses 180, open ocean
    # masked
    elements_plot, mask_patches = make_patches(mesh_path, circumpolar=True, mask_cavities=True)
    # Unmask ice shelf elements
    patches = iceshelf_mask(elements_plot)
    # Contour ice shelf front
    contour_lines = []
    for elm in elements_plot:
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

    print 'Reading data'
    ismr_elm = zeros([num_expts+1, num_cavity_elm])
    # Loop over experiments
    for expt in range(num_expts+1):
        if expt == 0:
            # 1996-2005
            id = Dataset(directory_beg + file_beg, 'r')            
        else:
            # RCP
            id = Dataset(directories[expt-1] + file_end, 'r')
        ismr_nodes = id.variables['wnet'][0,:]*sec_per_year
        id.close()
        # Loop over elements
        i = 0
        for elm in elements_calc:
            # For each element in an ice shelf cavity, save the mean value
            # for the 3 component Nodes
            if elm.cavity:
                ismr_elm[expt,i] = mean([ismr_nodes[elm.nodes[0].id], ismr_nodes[elm.nodes[1].id], ismr_nodes[elm.nodes[2].id]])
                i += 1

    print 'Integrating mass loss'
    total_massloss = zeros([num_expts+1, num_sectors])
    # Loop over elements
    i = 0
    for elm in elements_calc:
        if elm.cavity:
            # Figure out which sector this ice shelf element falls into
            # First get average lon and lat across 3 Nodes
            lon = mean(elm.lon)
            lat = mean(elm.lat)
            if lon >= -85 and lon < -30 and lat < -74:
                # Filchner-Ronne
                index = 0            
            elif lon >= -30 and lon < 65:
                # Eastern Weddell region
                index = 1            
            elif lon >= 65 and lon < 76:
                # Amery
                index = 2            
            elif lon >= 76 and lon < 165 and lat >= -74:
                # Australian sector
                index = 3            
            elif (lon >= 155 and lon < 165 and lat < -74) or (lon >= 165) or (lon < -140):
                # Ross Sea
                index = 4            
            elif (lon >= -140 and lon < -105) or (lon >= -105 and lon < -98 and lat < -73.1):
                # Amundsen Sea
                index = 5            
            elif (lon >= -104 and lon < -98 and lat >= -73.1) or (lon >= -98 and lon < -66 and lat >= -75):
                # Bellingshausen Sea
                index = 6            
            elif lon >= -66 and lon < -59 and lat >= -74:
                # Larsen Ice Shelves
                index = 7
            else:
                print 'No region found for lon=',str(lon),', lat=',str(lat)
                break #return
            # Integrate total mass loss in this sector
            for expt in range(num_expts+1):
                total_massloss[expt,index] += ismr_elm[expt,i]*elm.area()*rho_ice*1e-12
            i += 1

    # Calculate change in mass loss (absolute and percent) for each sector and
    # each RCP; print to screen
    abs_change = zeros([num_expts, num_sectors])
    percent_change = zeros([num_expts, num_sectors])
    for index in range(num_sectors):
        print sector_names[index]
        # Get mass loss for 1996-2005
        massloss_beg = total_massloss[0,index]
        for expt in range(num_expts):
            # Get mass loss for 2091-2100 in this RCP
            massloss_end = total_massloss[expt+1,index]
            abs_change[expt,index] = massloss_end-massloss_beg
            percent_change[expt,index] = (massloss_end-massloss_beg)/massloss_beg*1e2
            print expt_names[expt] + ': ' + str(abs_change[expt,index]) + ' (' + str(percent_change[expt,index]) + ')'
    # Also calculate total Antarctica
    abs_change_all = zeros(num_expts)
    percent_change_all = zeros(num_expts)
    print 'Total Antarctica'
    massloss_beg = sum(total_massloss[0,:])
    for expt in range(num_expts):
        massloss_end = sum(total_massloss[expt+1,:])
        abs_change_all[expt] = massloss_end-massloss_beg
        percent_change_all[expt] = (massloss_end-massloss_beg)/massloss_beg*1e2
        print expt_names[expt] + ': ' + str(abs_change_all[expt]) + ' (' + str(percent_change_all[expt]) + ')'

    print 'Setting up sectors for plotting'
    patch_sector = zeros(len(patches))
    # Loop over elements (plotting elements this time, which may cross 180)
    i = 0
    for elm in elements_plot:
        if elm.cavity:
            # Figure out which sector this ice shelf element falls into
            # First get average lon and lat across 3 Nodes
            lon = mean(elm.lon)
            lat = mean(elm.lat)
            if lon >= -85 and lon < -30 and lat < -74:
                # Filchner-Ronne
                index = 0            
            elif lon >= -30 and lon < 65:
                # Eastern Weddell region
                index = 1            
            elif lon >= 65 and lon < 76:
                # Amery
                index = 2            
            elif lon >= 76 and lon < 165 and lat >= -74:
                # Australian sector
                index = 3            
            elif (lon >= 155 and lon < 165 and lat < -74) or (lon >= 165) or (lon < -140):
                # Ross Sea
                index = 4            
            elif (lon >= -140 and lon < -105) or (lon >= -105 and lon < -98 and lat < -73.1):
                # Amundsen Sea
                index = 5            
            elif (lon >= -104 and lon < -98 and lat >= -73.1) or (lon >= -98 and lon < -66 and lat >= -75):
                # Bellingshausen Sea
                index = 6            
            elif lon >= -66 and lon < -59 and lat >= -74:
                # Larsen Ice Shelves
                index = 7
            else:
                print 'No region found for lon=',str(lon),', lat=',str(lat)
                break #return
            # Save this value
            patch_sector[i] = index
            i += 1
    # Set up colour map
    sector_cmap = ListedColormap(sector_colours)
    bounds = arange(-0.5, num_sectors)
    norm = BoundaryNorm(bounds, sector_cmap.N)
    # Set up a grey square covering the domain, anything that isn't covered
    # up later is land
    x_reg, y_reg = meshgrid(linspace(-nbdry, nbdry, num=100), linspace(-nbdry, nbdry, num=100))
    land_square = zeros(shape(x_reg))

    print 'Plotting'
    fig = figure(figsize=(12,12))
    fig.patch.set_facecolor('white')
    gs = GridSpec(2,2)
    gs.update(left=0.05, right=0.95, bottom=0, top=0.88, wspace=0.05, hspace=0.05)
    for expt in range(num_expts):
        ax = subplot(gs[expt/2, expt%2], aspect='equal')
        # Start with grey square background for land
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.7', '0.7', '0.7')))
        img = PatchCollection(patches, cmap=sector_cmap)
        img.set_array(array(patch_sector))
        img.set_edgecolor('face')
        img.set_norm(norm)
        ax.add_collection(img)
        # Mask out the open ocean in white
        overlay = PatchCollection(mask_patches, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        # Contour ice shelf fronts
        contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
        ax.add_collection(contours)
        # Add circles
        circles = []
        for index in range(num_sectors):
            if key == 1:
                circles.append(Circle((x_circles[index], y_circles[index]), radius=k*sqrt(abs_change[expt,index])))
            elif key == 2:
                circles.append(Circle((x_circles[index], y_circles[index]), radius=k*sqrt(percent_change[expt,index])))
        circle_patches = PatchCollection(circles, cmap=sector_cmap)
        circle_patches.set_array(arange(num_sectors))
        circle_patches.set_edgecolor('black')
        circle_patches.set_norm(norm)
        ax.add_collection(circle_patches)
        # Add labels on top of circles
        for index in range(num_sectors):
            if key == 1:
                text(x_circles[index], y_circles[index], str(int(round(abs_change[expt,index]))), ha='center', va='center', fontsize=16)
            elif key == 2:
                text(x_circles[index], y_circles[index], str(int(round(percent_change[expt,index]))), ha='center', va='center', fontsize=16)
        # Label for the whole continent
        if key == 1:
            text(3, 2, str(int(round(abs_change_all[expt]))), ha='center', va='center', fontsize=24)
        elif key == 2:
            text(3, 2, str(int(round(percent_change_all[expt]))), ha='center', va='center', fontsize=24)
        xlim([-nbdry, nbdry])
        ylim([-nbdry+8, nbdry])
        axis('off')
        title(letters[expt] + ' ' + expt_names[expt], fontsize=24)
    if key == 1:
        suptitle('Change in ice shelf basal mass loss (Gt/y),\n2091-2100 minus 1996-2005', fontsize=30)
    elif key == 2:
        suptitle('Percent change in ice shelf basal mass loss by sector,\n2091-2100 with respect to 1996-2005', fontsize=30)
    fig.show()
    if key == 1:
        fig.savefig('rcp_map_circles.png')
    elif key == 2:
        fig.savefig('rcp_map_circles_percent.png')


# Command-line interface
if __name__ == "__main__":

    key = int(raw_input('Plot absolute change (1) or percent change (2)? '))
    rcp_map_circles(key)
        



