from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from patches import *

def aice_seasonal_rcp ():

    directory_head = '/short/y99/kaa561/FESOM/'
    mesh_path = directory_head + 'mesh/low_res/'
    expt_paths = ['rcp45_M', 'rcp45_A', 'rcp85_M', 'rcp85_A']
    expt_titles = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A']
    output_dir = '/output/'
    filenames = ['avg.ice.mean.2006.2015.nc', 'avg.ice.mean.2091.2100.nc']

    # FESOM parameters
    circumpolar=True
    mask_cavities=True
    # Season names for plot titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']

    # Build FESOM mesh
    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)

    for expt in range(4):
        # Get seasonal averages of the FESOM output
        # This is hard-coded and ugly
        # First do 2006-2015 climatology
        file1 = directory_head + expt_paths[expt] + output_dir + filenames[0]
        id = Dataset(file1, 'r')
        n2d = id.variables['area'].shape[1]
        aice1 = zeros([4, n2d])
        # DJF: 1/5 of index 67 (1-based) and indices 68-73, plus indices 1-11
         # and 4/5 of index 12; 90 days in total
        aice1[0,:] = id.variables['area'][66,:] + sum(id.variables['area'][67:73,:]*5, axis=0)
        aice1[0,:] += sum(id.variables['area'][0:11,:]*5, axis=0) + id.variables['area'][11,:]*4
        aice1[0,:] /= 90
        # MAM: 1/5 of index 12, indices 13-30, and 1/5 of index 31;
        # 92 days in total
        aice1[1,:] = id.variables['area'][11,:] + sum(id.variables['area'][12:30,:]*5, axis=0) + id.variables['area'][30,:]
        aice1[1,:] /= 92
        # JJA: 4/5 of index 31, indices 32-48, and 3/5 of index 49;
        # 92 days in total
        aice1[2,:] = id.variables['area'][30,:]*4 + sum(id.variables['area'][31:48,:]*5, axis=0) + id.variables['area'][48,:]*3
        aice1[2,:] /= 92
        # SON: 2/5 of index 49, indices 50-66, and 4/5 of index 67;
        # 91 days in total
        aice1[3,:] = id.variables['area'][48,:]*2 + sum(id.variables['area'][49:66,:]*5, axis=0) + id.variables['area'][66,:]*4
        aice1[3,:] /= 91
        id.close()
        # Repeat for 2091-2100 climatology
        file2 = directory_head + expt_paths[expt] + output_dir + filenames[1]
        id = Dataset(file2, 'r')
        n2d = id.variables['area'].shape[1]
        aice2 = zeros([4, n2d])
        # DJF: 1/5 of index 67 (1-based) and indices 68-73, plus indices 1-11
         # and 4/5 of index 12; 90 days in total
        aice2[0,:] = id.variables['area'][66,:] + sum(id.variables['area'][67:73,:]*5, axis=0)
        aice2[0,:] += sum(id.variables['area'][0:11,:]*5, axis=0) + id.variables['area'][11,:]*4
        aice2[0,:] /= 90
        # MAM: 1/5 of index 12, indices 13-30, and 1/5 of index 31;
        # 92 days in total
        aice2[1,:] = id.variables['area'][11,:] + sum(id.variables['area'][12:30,:]*5, axis=0) + id.variables['area'][30,:]
        aice2[1,:] /= 92
        # JJA: 4/5 of index 31, indices 32-48, and 3/5 of index 49;
        # 92 days in total
        aice2[2,:] = id.variables['area'][30,:]*4 + sum(id.variables['area'][31:48,:]*5, axis=0) + id.variables['area'][48,:]*3
        aice2[2,:] /= 92
        # SON: 2/5 of index 49, indices 50-66, and 4/5 of index 67;
        # 91 days in total
        aice2[3,:] = id.variables['area'][48,:]*2 + sum(id.variables['area'][49:66,:]*5, axis=0) + id.variables['area'][66,:]*4
        aice2[3,:] /= 91
        id.close()

        # Plot
        fig = figure(figsize=(20,9))
        # Loop over seasons
        for season in range(4):
            # 2006-2015 climatology
            # Build an array of FESOM data values corresponding to each Element
            values1 = []
            for elm in elements:
                # For each element not in an ice shelf cavity, append the mean
                # value for the 3 component Nodes
                if not elm.cavity:
                    values1.append(mean([aice1[season,elm.nodes[0].id], aice1[season,elm.nodes[1].id], aice1[season,elm.nodes[2].id]]))
            ax = fig.add_subplot(2, 4, season+1, aspect='equal')
            img = PatchCollection(patches, cmap=jet)
            img.set_array(array(values1))
            img.set_clim(vmin=0, vmax=1)
            img.set_edgecolor('face')
            ax.add_collection(img)
            xlim([-35, 35])
            ylim([-33, 37])
            axis('off')
            if season == 0:
                text(-39, 0, '2006-2015', fontsize=21, ha='right')
            title(season_names[season], fontsize=24)
            # 2091-2100 climatology
            values2 = []
            for elm in elements:
                if not elm.cavity:
                    values2.append(mean([aice2[season,elm.nodes[0].id], aice2[season,elm.nodes[1].id], aice2[season,elm.nodes[2].id]]))
            ax = fig.add_subplot(2, 4, season+5, aspect='equal')
            img = PatchCollection(patches, cmap=jet)
            img.set_array(array(values2))
            img.set_clim(vmin=0, vmax=1)
            img.set_edgecolor('face')
            ax.add_collection(img)
            xlim([-35, 35])
            ylim([-33, 37])
            axis('off')
            if season == 0:
                text(-39, 0, '2091-2100', fontsize=21, ha='right')
        # Add a horizontal colorbar at the bottom
        cbaxes = fig.add_axes([0.25, 0.04, 0.5, 0.02])
        cbar = colorbar(img, orientation='horizontal', ticks=arange(0,1+0.25,0.25), cax=cbaxes)
        cbar.ax.tick_params(labelsize=16)
        # Add the main title
        suptitle('Sea ice concentration (' + expt_titles[expt] + ')', fontsize=30)
        # Decrease space between plots
        subplots_adjust(wspace=0.025,hspace=0.025)
        fig.savefig('aice_' + expt_paths[expt] + '.png')
        

# Command-line interface
if __name__ == "__main__":

    aice_seasonal_rcp()

    

    

    
