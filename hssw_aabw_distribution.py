from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from fesom_grid import *

def hssw_aabw_distribution ():

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/high_res/'
    directory_beg = '/short/y99/kaa561/FESOM/highres_spinup/'
    directories = ['/short/y99/kaa561/FESOM/rcp45_M_highres/output/', '/short/y99/kaa561/FESOM/rcp45_A_highres/output/', '/short/y99/kaa561/FESOM/rcp85_M_highres/output/', '/short/y99/kaa561/FESOM/rcp85_A_highres/output/', '/short/y99/kaa561/FESOM/highres_spinup/']
    file_beg = 'annual_avg.oce.mean.1996.2005.nc'
    file_end = 'annual_avg.oce.mean.2091.2100.nc'
    # Titles for plotting
    expt_names = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A', 'CONTROL']
    num_expts = len(directories)
    # Mesh parameters
    circumpolar = False
    cross_180 = False
    # Northern boundary of water masses to consider
    nbdry = -65
    # Number of temperature and salinity bins
    num_bins = 1000
    # Bounds on temperature and salinity bins (pre-computed, change if needed)
    min_salt = 32.3
    max_salt = 35.1
    min_temp = -3.1
    max_temp = 3.8
    # Bounds to plot for HSSW and AABW
    hssw_salt_bounds = [34.3, 35]
    hssw_temp_bounds = [-2.25, -0.75]
    aabw_salt_bounds = [34.55, 34.8]
    aabw_temp_bounds = [-1, 2.5]
    # More readable labels
    hssw_salt_ticks = arange(34.3, 35+0.1, 0.1)
    hssw_salt_labels = ['', '34.4', '', '34.6', '', '34.8', '', '35']
    hssw_temp_ticks = arange(-2.25, -0.75+0.25, 0.25)
    hssw_temp_labels = ['', '-2', '', '-1.5', '', '-1', '']
    aabw_salt_ticks = arange(34.55, 34.8+0.05, 0.05)
    aabw_salt_labels = ['', '34.6', '', '34.7', '', '34.8']
    aabw_temp_ticks = arange(-1, 2.5+0.5, 0.5)
    aabw_temp_labels = ['-1', '', '0', '', '1', '', '2', '']

    print 'Setting up bins'
    # Calculate boundaries of temperature bins
    temp_bins = linspace(min_temp, max_temp, num=num_bins)
    # Calculate centres of temperature bins (for plotting)
    temp_centres = 0.5*(temp_bins[:-1] + temp_bins[1:])
    # Repeat for salinity
    salt_bins = linspace(min_salt, max_salt, num=num_bins)
    salt_centres = 0.5*(salt_bins[:-1] + salt_bins[1:])
    # Set up a 3D array of experiment x temperature bins x salinity bins to
    # increment with volume of water masses
    ts_vals = zeros([num_expts+1, size(temp_centres), size(salt_centres)])
    # Calculate surface freezing point as a function of salinity as seen by
    # sea ice model
    freezing_pt = -0.0575*salt_centres + 1.7105e-3*sqrt(salt_centres**3) - 2.155e-4*salt_centres**2

    print 'Building mesh'
    elements = fesom_grid(mesh_path, circumpolar, cross_180)

    print 'Reading data'
    # 1996-2005
    id = Dataset(directory_beg + file_beg)
    n3d = id.variables['temp'].shape[1]
    temp_nodes = empty([num_expts+1, n3d])
    salt_nodes = empty([num_expts+1, n3d])
    temp_nodes[0,:] = id.variables['temp'][0,:]
    salt_nodes[0,:] = id.variables['salt'][0,:]
    id.close()
    # Loop over RCPs
    for expt in range(num_expts):
        id = Dataset(directories[expt] + file_end)
        temp_nodes[expt+1,:] = id.variables['temp'][0,:]
        salt_nodes[expt+1,:] = id.variables['salt'][0,:]
        id.close()

    print 'Binning elements'
    for elm in elements:
        # See if we're in the region of interest
        if all(elm.lat < nbdry):
            # Get area of 2D triangle
            area = elm.area()
            nodes = [elm.nodes[0], elm.nodes[1], elm.nodes[2]]
            # Loop downward
            while True:
                if nodes[0].below is None or nodes[1].below is None or nodes[2].below is None:
                    # We've reached the bottom
                    break
                # Calculate average temperature and salinity for each
                # experiment, as well as layer thickness, over this 3D
                # triangular prism.
                temp_vals = empty([num_expts+1, 6])
                salt_vals = empty([num_expts+1, 6])
                dz = empty(3)
                for i in range(3):
                    # Loop over experiments
                    for expt in range(num_expts+1):
                        # Average temperature over 6 nodes
                        temp_vals[expt,i] = temp_nodes[expt,nodes[i].id]
                        temp_vals[expt,i+3] = temp_nodes[expt,nodes[i].below.id]
                        salt_vals[expt,i] = salt_nodes[expt,nodes[i].id]
                        salt_vals[expt,i+3] = salt_nodes[expt,nodes[i].below.id]
                    # Average dz over 3 vertical edges
                    dz[i] = abs(nodes[i].depth - nodes[i].below.depth)
                    # Get ready for next repetition of loop
                    nodes[i] = nodes[i].below
                temp_elm = mean(temp_vals, axis=1)
                salt_elm = mean(salt_vals, axis=1)
                # Calculate volume of 3D triangular prism
                volume = area*mean(dz)
                # Loop over experiments again
                for expt in range(num_expts+1):
                    # Figure out which bins this falls into
                    temp_index = nonzero(temp_bins > temp_elm[expt])[0][0] - 1
                    salt_index = nonzero(salt_bins > salt_elm[expt])[0][0] - 1
                    # Increment bins with volume
                    ts_vals[expt, temp_index, salt_index] += volume
    # Mask bins with zero volume
    ts_vals = ma.masked_where(ts_vals==0, ts_vals)

    # Find the volume bounds for plotting
    min_val = log(amin(ts_vals))
    max_val = log(amax(ts_vals))

    print 'Plotting'
    fig = figure(figsize=(20,11))
    # HSSW
    gs_a = GridSpec(1,num_expts+1)
    gs_a.update(left=0.05, right=0.98, bottom=0.62, top=0.86, wspace=0.16)
    for expt in range(num_expts+1):
        ax = subplot(gs_a[0,expt])
        # Log scale is more visible
        img = pcolor(salt_centres, temp_centres, log(ts_vals[expt,:,:]), vmin=min_val, vmax=max_val, cmap='jet')
        plot(salt_centres, freezing_pt, color='black', linestyle='dashed', linewidth=2)
        grid(True)
        xlim(hssw_salt_bounds)
        ylim(hssw_temp_bounds)
        ax.set_xticks(hssw_salt_ticks)
        ax.set_xticklabels(hssw_salt_labels)
        ax.set_yticks(hssw_temp_ticks)
        ax.set_yticklabels(hssw_temp_labels)
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        # Labels and titles
        if expt == 0:
            xlabel('Salinity (psu)', fontsize=16)
            ylabel(r'Temperature ($^{\circ}$C)', fontsize=16)
            title('1996-2005', fontsize=19)
        elif expt == 1:
            title(expt_names[expt-1] + ' (2091-2100)', fontsize=19)
        else:
            title(expt_names[expt-1], fontsize=19)
        # HSSW title
        if expt == 2:
            text(34.83, hssw_temp_bounds[1]+0.3, 'a) HSSW', ha='left', fontsize=30)
    # AABW
    gs_b = GridSpec(1,num_expts+1)
    gs_b.update(left=0.05, right=0.98, bottom=0.12, top=0.5, wspace=0.16)
    for expt in range(num_expts+1):
        ax = subplot(gs_b[0,expt])
        img = pcolor(salt_centres, temp_centres, log(ts_vals[expt,:,:]), vmin=min_val, vmax=max_val, cmap='jet')
        grid(True)
        xlim(aabw_salt_bounds)
        ylim(aabw_temp_bounds)
        ax.set_xticks(aabw_salt_ticks)
        ax.set_xticklabels(aabw_salt_labels)
        ax.set_yticks(aabw_temp_ticks)
        ax.set_yticklabels(aabw_temp_labels)
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        if expt == 0:
            xlabel('Salinity (psu)', fontsize=16)
            ylabel(r'Temperature ($^{\circ}$C)', fontsize=16)
            title('1996-2005', fontsize=19)
        elif expt == 1:
            title(expt_names[expt-1] + ' (2091-2100)', fontsize=19)
        else:
            title(expt_names[expt-1], fontsize=19)
        # AABW title
        if expt == 2:
            text(34.71, aabw_temp_bounds[1]+0.4, 'b) AABW', ha='left', fontsize=30)
        # Horizontal colourbar at the bottom
        if expt == num_expts:
            cbaxes = fig.add_axes([0.35, 0.06, 0.3, 0.02])
            cbar = colorbar(img, cax=cbaxes, orientation='horizontal')
            cbar.ax.tick_params(labelsize=14)
            text(0.5, 0.01, 'log of volume', fontsize=20, transform=fig.transFigure, ha='center')
    # Main title
    suptitle(r'Water masses south of 65$^{\circ}$S', fontsize=30)
    fig.show()
    fig.savefig('hssw_aabw_distribution.png')


# Command-line interface
if __name__ == "__main__":

    hssw_aabw_distribution()
            
    
    
    
    
