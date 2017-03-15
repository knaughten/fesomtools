from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.patches import Polygon
from patches import *
from unrotate_vector import *

def rcp_sst ():

    mesh_path = '../FESOM/mesh/low_res/'
    expt_names = ['rcp45_M', 'rcp85_M']
    expt_titles = ['RCP 4.5', 'RCP 8.5']
    file_middle = '/output/avg.oce.mean.'
    start_year = [2006, 2091]
    end_year = [2015, 2100]
    min_abs = -2
    max_abs = 8
    min_anom = -2.5
    max_anom = 2.5
    circumpolar = True
    mask_cavities = True
    lat_max = -60+90

    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)
    mask_patches = iceshelf_mask(elements)
    # Set up a grey square covering the domain, anything that isn't covered
    # up later is land
    x_reg, y_reg = meshgrid(linspace(-lat_max, lat_max, num=100), linspace(-lat_max, lat_max, num=100))
    land_square = zeros(shape(x_reg))

    fig = figure(figsize=(20,12))
    for expt_number in range(len(expt_names)):            
        expt = expt_names[expt_number]
        file_start = '../FESOM/' + expt + file_middle + str(start_year[0]) + '.' + str(end_year[0]) + '.nc'
        id = Dataset(file_start, 'r')
        n3d = id.variables['temp'].shape[1]
        data = ma.empty([3, n3d])
        # DJF: 1/5 of index 67 (1-based), indices 68-73, indices 1-11, and
        # 4/5 of index 12 in 5-day climatology file; 90 days in total
        data[0,:] = (id.variables['temp'][66,:] + sum(id.variables['temp'][67:73,:]*5, axis=0) + sum(id.variables['temp'][0:11,:]*5, axis=0) + id.variables['temp'][11,:]*4)/90.0
        id.close()
        file_end = '../FESOM/' + expt + file_middle + str(start_year[1]) + '.' + str(end_year[1]) + '.nc'
        id = Dataset(file_end, 'r')
        data[1,:] = (id.variables['temp'][66,:] + sum(id.variables['temp'][67:73,:]*5, axis=0) + sum(id.variables['temp'][0:11,:]*5, axis=0) + id.variables['temp'][11,:]*4)/90.0
        id.close()
        data[2,:] = data[1,:] - data[0,:]
        for panel in range(3):
            values = []
            for elm in elements:
                if (mask_cavities and not elm.cavity) or (not mask_cavities):
                    # Surface nodes; this is easy
                    # Average the data value for each of the three component
                    # nodes
                    values.append(mean([data[panel,elm.nodes[0].id], data[panel,elm.nodes[1].id], data[panel,elm.nodes[2].id]]))
            if panel == 2:
                var_min = min_anom
                var_max = max_anom
                colour_map = 'RdBu_r'
            else:
                var_min = min_abs
                var_max = max_abs
                colour_map = 'jet'
            ax = fig.add_subplot(2, 3, expt_number*3+panel+1, aspect='equal')
            contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
            img = PatchCollection(patches, cmap=colour_map)
            img.set_array(array(values))
            img.set_edgecolor('face')
            img.set_clim(vmin=var_min, vmax=var_max)
            ax.add_collection(img)
            if mask_cavities:
                overlay = PatchCollection(mask_patches, facecolor=(0.8, 0.8, 0.8))
                overlay.set_edgecolor('face')
                ax.add_collection(overlay)
            xlim([-lat_max, lat_max])
            ylim([-lat_max, lat_max])
            axis('off')
            if panel == 2:
                title(expt_titles[expt_number] + ' (difference)', fontsize=22)
            else:
                title(expt_titles[expt_number] + ' (' + str(start_year[panel]) + '-' + str(end_year[panel]) + ' average)', fontsize=22)
            if expt_number == 0:
                if panel != 1:
                    tmp = [0, 0.25, 0.025, 0.5]
                    if panel == 0:
                        tmp[0] = 0.05
                    elif panel == 2:
                        tmp[0] = 0.92
                    cbaxes = fig.add_axes(tmp)
                    cbar = colorbar(img, cax=cbaxes, extend='both')
                    cbar.ax.tick_params(labelsize=16)
    suptitle(r'DJF Sea Surface Temperature ($^{\circ}$C)', fontsize=30)
    #fig.show()
    fig.savefig('sst_rcp.png')


if __name__ == "__main__":

    rcp_sst()

        
                
                
                    
                    
        
