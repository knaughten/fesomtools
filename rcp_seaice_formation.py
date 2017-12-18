from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.colors import *
from matplotlib.cm import *
from patches import *

def rcp_seaice_formation (rcp, model):

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/meshB/'
    file_beg = '/short/y99/kaa561/FESOM/highres_spinup/seasonal_climatology_ice_diag_1996_2005.nc'
    file_end = '/short/y99/kaa561/FESOM/rcp' + rcp + '_' + model + '/seasonal_climatology_ice_diag_2091_2100.nc'
    # Bound on plot
    nbdry = -64+90
    if model == 'M':
        rcp_title = 'RCP ' + rcp[0] + '.' + rcp[1] + ' MMM'
    elif model == 'A':
        rcp_title = 'RCP ' + rcp[0] + '.' + rcp[1] + ' ACCESS'
    # Mesh parameters
    circumpolar = True
    mask_cavities = True
    # Season names for plot titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    # m/s to m/y conversion factor
    sec_per_year = 365.25*24*60*60
    # Bounds on colour scale
    bound_abs = 20
    bound_anom = 4

    print 'Building mesh'
    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)

    print 'Reading data'
    # Initial thdgr, convert from m/s to m/y
    id = Dataset(file_beg, 'r')
    thdgr_nodes_beg = id.variables['thdgr'][:,:]*sec_per_year
    id.close()
    # Read anomalies in RCP
    id = Dataset(file_end, 'r')
    thdgr_nodes_diff = id.variables['thdgr'][:,:]*sec_per_year - thdgr_nodes_beg
    id.close()

    print 'Calculating element averages'
    thdgr_beg = zeros([4, len(patches)])
    thdgr_diff = zeros([4, len(patches)])
    i = 0
    for elm in elements:
        if not elm.cavity:
            thdgr_beg[:,i] = (thdgr_nodes_beg[:,elm.nodes[0].id] + thdgr_nodes_beg[:,elm.nodes[1].id] + thdgr_nodes_beg[:,elm.nodes[2].id])/3.0
            thdgr_diff[:,i] = (thdgr_nodes_diff[:,elm.nodes[0].id] + thdgr_nodes_diff[:,elm.nodes[1].id] + thdgr_nodes_diff[:,elm.nodes[2].id])/3.0
            i += 1

    # Make a nonlinear colour scale for initial values
    bounds = empty(100)
    bounds[:50] = -1*linspace(bound_abs**(1.0/2), 0, num=50)**2
    bounds[50:] = linspace(0, bound_abs**(1.0/2), num=50)**2
    norm = BoundaryNorm(boundaries=bounds, ncolors=256)

    print 'Plotting'
    fig = figure(figsize=(8,20))
    gs = GridSpec(4,2)
    gs.update(left=0.12, right=0.95, bottom=0.04, top=0.92, wspace=0.025, hspace=0.025)
    for season in range(4):
        # 1996-2005
        ax = subplot(gs[season,0], aspect='equal')
        img = PatchCollection(patches, cmap='PRGn', norm=norm)
        img.set_array(thdgr_beg[season,:])
        img.set_clim(-bound_abs, bound_abs)
        img.set_edgecolor('face')
        ax.add_collection(img)
        xlim([-nbdry, nbdry])
        ylim([-nbdry, nbdry])
        ax.set_xticks([])
        ax.set_yticks([])
        if season == 0:
            title('a) 1996-2005', fontsize=24)
        text(-nbdry-2, 0, season_names[season], fontsize=24, ha='right')
        if season == 3:
            cbaxes1 = fig.add_axes([0.145, 0.02, 0.365, 0.01])
            cbar1 = colorbar(img, orientation='horizontal', cax=cbaxes1, extend='both', ticks=arange(-bound_abs, bound_abs+10, 10))
            cbar1.ax.tick_params(labelsize=16)
        # 2091-2100 anomalies
        ax = subplot(gs[season,1], aspect='equal')
        img = PatchCollection(patches, cmap='RdBu_r')
        img.set_array(thdgr_diff[season,:])
        img.set_clim(-bound_anom, bound_anom)
        img.set_edgecolor('face')
        ax.add_collection(img)
        xlim([-nbdry, nbdry])
        ylim([-nbdry, nbdry])
        ax.set_xticks([])
        ax.set_yticks([])
        if season == 0:
            title('b) ' + rcp_title, fontsize=24)
            text(nbdry+2, -2.7*nbdry, 'anomalies (2091-2100 minus 1996-2005)', fontsize=20, ha='left', va='bottom', rotation=270)
        if season == 3:
            cbaxes2 = fig.add_axes([0.565, 0.02, 0.365, 0.01])
            cbar2 = colorbar(img, orientation='horizontal', cax=cbaxes2, extend='both', ticks=arange(-bound_anom, bound_anom+2, 2))
            cbar2.ax.tick_params(labelsize=16)
    suptitle('Net sea ice formation (+) or melt (-), m/y', fontsize=27)
    fig.show()
    fig.savefig('seaice_formation_rcp'+rcp+'_'+model+'.png')


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
    rcp_seaice_formation(rcp, model)

    
