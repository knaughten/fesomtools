from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.patches import Polygon
from patches import *
from unrotate_vector import *

def rcp_bw_ts ():

    mesh_path = '../FESOM/mesh/low_res/'
    var_names = ['temp', 'salt']
    var_titles = [r'temperature ($^{\circ}$C)', 'salinity (psu)']
    expt_names = ['rcp45_M', 'rcp85_M']
    expt_titles = ['RCP 4.5', 'RCP 8.5']
    file_middle = '/output/avg.oce.mean.'
    start_year = [2006, 2091]
    end_year = [2015, 2100]
    min_abs = [-2.5, 34]
    max_abs = [1.5, 35]
    min_anom = [-1.5, -0.4]
    max_anom = [1.5, 0.4]
    circumpolar = True
    mask_cavities = False
    lat_max = -60+90

    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)
    # Set up a grey square covering the domain, anything that isn't covered
    # up later is land
    x_reg, y_reg = meshgrid(linspace(-lat_max, lat_max, num=100), linspace(-lat_max, lat_max, num=100))
    land_square = zeros(shape(x_reg))

    for var_number in range(len(var_names)):
        var = var_names[var_number]
        fig = figure(figsize=(20,12))
        for expt_number in range(len(expt_names)):            
            expt = expt_names[expt_number]
            file_start = '../FESOM/' + expt + file_middle + str(start_year[0]) + '.' + str(end_year[0]) + '.nc'
            id = Dataset(file_start, 'r')
            n3d = id.variables[var].shape[1]
            data = ma.empty([3, n3d])
            data[0,:] = mean(id.variables[var][:,:], axis=0)
            id.close()
            file_end = '../FESOM/' + expt + file_middle + str(start_year[1]) + '.' + str(end_year[1]) + '.nc'
            id = Dataset(file_end, 'r')
            data[1,:] = mean(id.variables[var][:,:], axis=0)
            id.close()
            data[2,:] = data[1,:] - data[0,:]
            for panel in range(3):
                values = []
                for elm in elements:
                    values_tmp = []
                    # For each of the three component nodes, find the id of
                    # the bottom node beneath it
                    for node in elm.nodes:
                        id = node.find_bottom().id
                        values_tmp.append(data[panel,id])
                    values.append(mean(values_tmp))
                if panel == 2:
                    var_min = min_anom[var_number]
                    var_max = max_anom[var_number]
                    colour_map = 'RdBu_r'
                else:
                    var_min = min_abs[var_number]
                    var_max = max_abs[var_number]
                    colour_map = 'jet'
                ax = fig.add_subplot(2, 3, expt_number*3+panel+1, aspect='equal')
                contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
                img = PatchCollection(patches, cmap=colour_map)
                img.set_array(array(values))
                img.set_edgecolor('face')
                img.set_clim(vmin=var_min, vmax=var_max)
                ax.add_collection(img)
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
        suptitle('Bottom water ' + var_titles[var_number], fontsize=30)
        #fig.show()
        fig.savefig('bw_' + var_names[var_number] + '_rcp.png')


if __name__ == "__main__":

    rcp_bw_ts()

        
                
                
                    
                    
        
