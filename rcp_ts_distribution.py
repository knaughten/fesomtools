from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.colors import *
from fesom_grid import *
from unesco import *

def rcp_ts_distribution (key=1):

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/high_res/'
    directory_beg = '/short/y99/kaa561/FESOM/highres_spinup/'
    directories = ['/short/y99/kaa561/FESOM/rcp45_M_highres/output/', '/short/y99/kaa561/FESOM/rcp45_A_highres/output/', '/short/y99/kaa561/FESOM/rcp85_M_highres/output/', '/short/y99/kaa561/FESOM/rcp85_A_highres/output/', '/short/y99/kaa561/FESOM/highres_spinup/']
    file_beg = 'annual_avg.oce.mean.1996.2005.nc'
    file_end = 'annual_avg.oce.mean.2091.2100.nc'
    # Titles for plotting
    expt_names = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A', 'CONTROL']
    num_expts = len(directories)
    # Start and end years for each period
    beg_years = [1996, 2005]
    end_years = [2091, 2100]
    # Northern boundary of water masses to consider
    nbdry = -65
    # Number of temperature and salinity bins
    num_bins = 1000
    # Bounds on temperature and salinity bins (pre-computed, change if needed)
    min_salt = 32.3
    max_salt = 35.1
    min_temp = -3.1
    max_temp = 3.8
    # Bounds to actually plot
    if key==1:
        min_salt_plot = 32.25
        max_salt_plot = 35
        min_temp_plot = -3
        max_temp_plot = 3.25
    elif key==2:
        min_salt_plot = 34
        max_salt_plot = 35
        min_temp_plot = -2.5
        max_temp_plot = -1
    # FESOM grid generation parameters
    circumpolar = False
    cross_180 = False

    print('Setting up bins')
    # Calculate boundaries of temperature bins
    temp_bins = linspace(min_temp, max_temp, num=num_bins)
    # Calculate centres of temperature bins (for plotting)
    temp_centres = 0.5*(temp_bins[:-1] + temp_bins[1:])
    # Repeat for salinity
    salt_bins = linspace(min_salt, max_salt, num=num_bins)
    salt_centres = 0.5*(salt_bins[:-1] + salt_bins[1:])
    # Set up 3D array of experiment x temperature bins x salinity bins to hold
    # average depth of water masses, weighted by volume
    ts_vals = zeros([num_expts+1, size(temp_centres), size(salt_centres)])
    # Also array to integrate volume of each bin
    volume = zeros([num_expts+1, size(temp_centres), size(salt_centres)])
    # Calculate surface freezing point as a function of salinity as seen by
    # sea ice model
    freezing_pt = -0.0575*salt_centres + 1.7105e-3*sqrt(salt_centres**3) - 2.155e-4*salt_centres**2
    # Get 2D versions of the temperature and salinity bins
    salt_2d, temp_2d = meshgrid(salt_centres, temp_centres)
    # Calculate potential density of each combination of temperature and
    # salinity bins
    density = unesco(temp_2d, salt_2d, zeros(shape(temp_centres)))-1000
    # Density contours to plot
    if key == 1:
        density_lev = arange(25.8, 28.4, 0.2)
    elif key == 2:
        density_lev = arange(27.2, 28.4, 0.2)

    print('Building grid')
    elements = fesom_grid(mesh_path, circumpolar, cross_180)

    print('Reading data')
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

    print('Binning elements')
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
                # experiment, as well as depth and layer thickness, over this
                # 3D triangular prism.
                temp_vals = empty([num_expts+1, 6])
                salt_vals = empty([num_expts+1, 6])
                depth_vals = empty(6)
                dz = empty(3)
                for i in range(3):
                    # Loop over experiments
                    for expt in range(num_expts+1):
                        # Average temperature over 6 nodes
                        temp_vals[expt,i] = temp_nodes[expt,nodes[i].id]
                        temp_vals[expt,i+3] = temp_nodes[expt,nodes[i].below.id]
                        salt_vals[expt,i] = salt_nodes[expt,nodes[i].id]
                        salt_vals[expt,i+3] = salt_nodes[expt,nodes[i].below.id]
                    # Average depth over 6 nodes
                    depth_vals[i] = nodes[i].depth
                    depth_vals[i+3] = nodes[i].below.depth
                    # Average dz over 3 vertical edges
                    dz[i] = abs(nodes[i].depth - nodes[i].below.depth)
                    # Get ready for next repetition of loop
                    nodes[i] = nodes[i].below
                temp_elm = mean(temp_vals, axis=1)
                salt_elm = mean(salt_vals, axis=1)
                depth_elm = mean(depth_vals)
                # Calculate volume of 3D triangular prism
                curr_volume = area*mean(dz)
                # Loop over experiments again
                for expt in range(num_expts+1):
                    # Figure out which bins this falls into
                    temp_index = nonzero(temp_bins > temp_elm[expt])[0][0] - 1
                    salt_index = nonzero(salt_bins > salt_elm[expt])[0][0] - 1
                    # Integrate depth*volume in this bin
                    ts_vals[expt, temp_index, salt_index] += depth_elm*curr_volume
                    volume[expt, temp_index, salt_index] += curr_volume
    # Mask bins with zero volume
    ts_vals = ma.masked_where(volume==0, ts_vals)
    volume = ma.masked_where(volume==0, volume)
    # Convert depths from integrals to volume-averages
    ts_vals /= volume

    # Find the maximum depth for plotting
    if key == 1:
        max_depth = amax(ts_vals)
    elif key == 2:
        temp_start = nonzero(temp_bins > min_temp_plot)[0][0]-2
        temp_end = nonzero(temp_bins > max_temp_plot)[0][0]
        salt_start = nonzero(salt_bins > min_salt_plot)[0][0]-2
        salt_end = nonzero(salt_bins > max_salt_plot)[0][0]
        max_depth = amax(ts_vals[:,temp_start:temp_end, salt_start:salt_end])
    # Make a nonlinear colour scale
    bounds = linspace(0, max_depth**(1.0/2.5), num=100)**2.5
    norm = BoundaryNorm(boundaries=bounds, ncolors=256)

    print('Plotting')
    fig = figure(figsize=(24,6))
    gs = GridSpec(1,num_expts+1)
    gs.update(left=0.04, right=0.99, bottom=0.12, top=0.86)
    for expt in range(num_expts+1):
        ax = subplot(gs[0,expt])
        img = pcolor(salt_centres, temp_centres, ts_vals[expt,:,:], norm=norm, vmin=0, vmax=max_depth, cmap='jet')
        plot(salt_centres, freezing_pt, color='black', linestyle='dashed')
        cs = contour(salt_centres, temp_centres, density, density_lev, colors=(0.6,0.6,0.6), linestyles='dotted')
        clabel(cs, inline=1, fontsize=10, color=(0.6,0.6,0.6), fmt='%1.1f')
        xlim([min_salt_plot, max_salt_plot])
        ylim([min_temp_plot, max_temp_plot])
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)
        if expt == 0:
            xlabel('Salinity (psu)', fontsize=14)
            ylabel(r'Temperature ($^{\circ}$C)', fontsize=14)
            title(str(beg_years[0]) + '-' + str(beg_years[1]), fontsize=20)
        elif expt == 1:
            title(expt_names[expt-1] + ' (' + str(end_years[0]) + '-' + str(end_years[1]) + ')', fontsize=20)
        else:
            title(expt_names[expt-1], fontsize=20)
        if expt == num_expts:
            # Add a horizontal colourbar below
            cbaxes = fig.add_axes([0.35, 0.05, 0.3, 0.02])
            if key == 1:
                cbar = colorbar(img, cax=cbaxes, orientation='horizontal', ticks=[0,50,100,200,500,1000,2000,4000])
            elif key == 2:
                cbar = colorbar(img, cax=cbaxes, orientation='horizontal', ticks=[0,50,100,200,500,1000,2000])
            cbar.ax.tick_params(labelsize=14)
    # Add the main title
    if key == 1:
        suptitle(r'Water masses south of 65$^{\circ}$S: depth (m)', fontsize=24)
    elif key == 2:
        suptitle(r'Water masses south of 65$^{\circ}$S, zoomed into HSSW: depth (m)', fontsize=24)
    fig.show()
    if key == 1:
        fig.savefig('ts_distribution_full.png')
    elif key ==2:
        fig.savefig('ts_distribution_hssw.png')


# Command-line interface
if __name__ == "__main__":

    key = int(input('Plot everything south of 65S (1) or just HSSW (2)? '))
    rcp_ts_distribution(key)
    
                    
    
    
    
