from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.pyplot import *
from patches import *

def ismr_change ():

    directory_head = '/short/y99/kaa561/FESOM/'
    mesh_path = directory_head + 'mesh/low_res/'
    expt_paths = ['rcp45_M', 'rcp45_A', 'rcp85_M', 'rcp85_A']
    expt_titles = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A']
    output_dir = '/output/'
    filenames = ['avg.forcing.diag.2006.2015.nc', 'avg.forcing.diag.2091.2100.nc']
    lat_max = -63 + 90
    circumpolar = True
    mask_cavities = True
    # Seconds per year
    sec_per_year = 365*24*3600

    # Build FESOM mesh
    # Get separate patches for the open ocean elements so we can mask them out
    elements, mask_patches = make_patches(mesh_path, circumpolar, mask_cavities)
    patches = iceshelf_mask(elements)

    values = empty([4, len(patches)])
    for expt in range(4):
        # Read freshwater flux for 2006-2015 climatology
        file1 = directory_head + expt_paths[expt] + output_dir + filenames[0]
        id = Dataset(file1, 'r')
        # Annually average and convert from m/s to m/y
        data1 = mean(id.variables['wnet'][:,:]*sec_per_year, axis=0)
        id.close()
        # Repeat for 2091-2100 climatology
        file2 = directory_head + expt_paths[expt] + output_dir + filenames[1]
        id = Dataset(file2, 'r')
        # Annually average and convert from m/s to m/y
        data2 = mean(id.variables['wnet'][:,:]*sec_per_year, axis=0)
        id.close()
        # Get difference
        data = data2 - data1
        values_tmp = []
        # Loop over elements
        for elm in elements:
            # For each element in an ice shelf cavity, append the mean value
            # for the 3 component Nodes
            if elm.cavity:
                values_tmp.append(mean([data[elm.nodes[0].id], data[elm.nodes[1].id], data[elm.nodes[2].id]]))
        values[expt,:] = values_tmp

    # Set up a grey square covering the domain, anything that isn't covered
    # up later is land
    x_reg, y_reg = meshgrid(linspace(-lat_max, lat_max, num=100), linspace(-lat_max, lat_max, num=100))
    land_square = zeros(shape(x_reg))

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

    max_val = 3.0

    fig = figure(figsize=(16,12))
    for expt in range(4):
        ax = fig.add_subplot(2, 2, expt+1, aspect='equal')
        # Start with grey square background for land
        contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(patches, cmap='YlOrRd')
        img.set_array(values[expt,:])
        img.set_edgecolor('face')
        img.set_clim(0, max_val)
        ax.add_collection(img)
        overlay = PatchCollection(mask_patches, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        contours = LineCollection(contour_lines, edgecolor='black', linewidth=1)
        ax.add_collection(contours)
        xlim([-lat_max, lat_max])
        ylim([-lat_max, lat_max])
        axis('off')
        title(expt_titles[expt], fontsize=24)
        if expt == 0:
            cbaxes = fig.add_axes([0.07, 0.55, 0.02, 0.3])
            cbar = colorbar(img, cax=cbaxes)
            cbar.ax.tick_params(labelsize=16)
    suptitle('Change in ice shelf melt rate (m/y), 2091-2100 vs 2006-2015', fontsize=30)
    subplots_adjust(wspace=0.025,hspace=0.1)

    #fig.show()
    fig.savefig('ismr_change.png')


if __name__ == "__main__":

    ismr_change()
        
