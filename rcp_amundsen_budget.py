from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from fesom_grid import *

def rcp_amundsen_budget ():

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/meshB/'
    directory_beg = '/short/y99/kaa561/FESOM/highres_spinup/'
    directories = ['/short/y99/kaa561/FESOM/rcp45_M/', '/short/y99/kaa561/FESOM/rcp45_A/', '/short/y99/kaa561/FESOM/rcp85_M/', '/short/y99/kaa561/FESOM/rcp85_A/']
    o_file_beg = 'avg.forcing.diag.1996.2005.nc'
    o_file_end = 'avg.forcing.diag.2091.2100.nc'
    i_file_beg = 'avg.ice.diag.1996.2005.nc'
    i_file_end = 'avg.ice.diag.2091.2100.nc'
    num_expts = len(directories)
    # Titles for plot
    beg_title = '1996-2005'
    rcp_titles = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A']
    # Colours for plotting
    beg_colour = (0,0,0)
    rcp_colours = [(0, 0.4, 1), (0.6, 0.2, 1), (0, 0.6, 0), (1, 0, 0.4)]
    # Days in each output step
    days_per_output = 5
    # Bounds on Amundsen Sea box
    lon_min = -115
    lon_max = -100
    lat_max = -71
    # FESOM mesh parameters
    circumpolar = True
    cross_180 = False

    print 'Building grid'
    elements = fesom_grid(mesh_path, circumpolar, cross_180)

    print 'Reading data'
    print '...' + beg_title
    id = Dataset(directory_beg + o_file_beg, 'r')
    wnet_nodes_beg = id.variables['wnet'][:,:]
    pminuse_nodes_beg = id.variables['rain'][:,:] + id.variables['snow'][:,:] - id.variables['evap'][:,:]    
    id.close()
    id = Dataset(directory_beg + i_file_beg, 'r')
    ice2ocn_nodes_beg = -1*id.variables['thdgr'][:,:]
    id.close()
    num_time = size(wnet_nodes_beg,0)
    n2d = size(wnet_nodes_beg,1)
    wnet_nodes_end = empty([num_expts, num_time, n2d])
    pminuse_nodes_end = empty([num_expts, num_time, n2d])
    ice2ocn_nodes_end = empty([num_expts, num_time, n2d])
    for expt in range(num_expts):
        print '...' + rcp_titles[expt]
        id = Dataset(directories[expt] + o_file_end, 'r')
        wnet_nodes_end[expt,:,:] = id.variables['wnet'][:,:]
        pminuse_nodes_end[expt,:,:] = id.variables['rain'][:,:] + id.variables['snow'][:,:] - id.variables['evap'][:,:]
        id.close()
        id = Dataset(directories[expt] + i_file_end, 'r')
        ice2ocn_nodes_end[expt,:,:] = -1*id.variables['thdgr'][:,:]
        id.close()

    print 'Averaging over Amundsen Sea'
    wnet_beg = zeros(num_time)
    pminuse_beg = zeros(num_time)
    ice2ocn_beg = zeros(num_time)
    wnet_end = zeros([num_expts, num_time])
    pminuse_end = zeros([num_expts, num_time])
    ice2ocn_end = zeros([num_expts, num_time])
    total_area = 0
    # Loop over elements
    for elm in elements:
        # Ignore ice shelf cavities
        if not elm.cavity:
            # Check if it's within the given lon and lat bounds
            if all(elm.lon >= lon_min) and all(elm.lon <= lon_max) and all(elm.lat <= lat_max):
                # Get area of triangle and integrate total area array
                area = elm.area()
                total_area += area
                # Process 1996-2005
                # Integrate each array with variable*area where variable is
                # averaged over 3 corner Nodes
                wnet_beg[:] += (wnet_nodes_beg[:,elm.nodes[0].id] + wnet_nodes_beg[:,elm.nodes[1].id] + wnet_nodes_beg[:,elm.nodes[2].id])/3.0*area
                pminuse_beg[:] += (pminuse_nodes_beg[:,elm.nodes[0].id] + pminuse_nodes_beg[:,elm.nodes[1].id] + pminuse_nodes_beg[:,elm.nodes[2].id])/3.0*area
                ice2ocn_beg[:] += (ice2ocn_nodes_beg[:,elm.nodes[0].id] + ice2ocn_nodes_beg[:,elm.nodes[1].id] + ice2ocn_nodes_beg[:,elm.nodes[2].id])/3.0*area
                # Process RCPs
                for expt in range(num_expts):
                    wnet_end[expt,:] += (wnet_nodes_end[expt,:,elm.nodes[0].id] + wnet_nodes_end[expt,:,elm.nodes[1].id] + wnet_nodes_end[expt,:,elm.nodes[2].id])/3.0*area
                    pminuse_end[expt,:] += (pminuse_nodes_end[expt,:,elm.nodes[0].id] + pminuse_nodes_end[expt,:,elm.nodes[1].id] + pminuse_nodes_end[expt,:,elm.nodes[2].id])/3.0*area
                    ice2ocn_end[expt,:] += (ice2ocn_nodes_end[expt,:,elm.nodes[0].id] + ice2ocn_nodes_end[expt,:,elm.nodes[1].id] + ice2ocn_nodes_end[expt,:,elm.nodes[2].id])/3.0*area
    # Convert variable arrays from integrals to averages (divide by total area)
    wnet_beg /= total_area
    pminuse_beg /= total_area
    ice2ocn_beg /= total_area
    wnet_end /= total_area
    pminuse_end /= total_area
    ice2ocn_end /= total_area

    # Set up time axis in days
    time = arange(num_time)*days_per_output

    print 'Plotting'
    fig = figure(figsize=(18,8))
    gs = GridSpec(1,3)
    gs.update(left=0.05, right=0.95, bottom=0.15, top=0.83, wspace=0.07)
    # Net freshwater flux, 1996-2005
    ax = subplot(gs[0,0])
    # Start with a horizontal line at zero
    ax.plot(time, zeros(size(time)), color=(0.6, 0.6, 0.6), linewidth=2)
    ax.plot(time, wnet_beg*1e8, color=beg_colour, linewidth=1.5)
    xlim([time[0], time[-1]])
    ylim([-14, 14])
    grid(True)
    xlabel('day of year', fontsize=18)
    ylabel(r'10$^{-8} $m/s', fontsize=18)
    title('a) Net freshwater flux\n1996-2005 climatology', fontsize=20)
    # RCP anomalies in net freshwater flux
    ax = subplot(gs[0,1])
    ax.plot(time, zeros(size(time)), color=(0.6, 0.6, 0.6), linewidth=2)
    for expt in range(num_expts):
        ax.plot(time, (wnet_end[expt,:]-wnet_beg)*1e8, color=rcp_colours[expt], linewidth=1.5)
    xlim([time[0], time[-1]])
    ylim([-14, 14])
    grid(True)
    xlabel('day of year', fontsize=18)
    title('b) Net freshwater flux\n2091-2100 minus 1996-2005', fontsize=20)
    # RCP anomalies in sea ice to ocean freshwater flux
    ax = subplot(gs[0,2])
    ax.plot(time, zeros(size(time)), color=(0.6, 0.6, 0.6), linewidth=2)
    for expt in range(num_expts):
        ax.plot(time, (ice2ocn_end[expt,:]-ice2ocn_beg)*1e8, color=rcp_colours[expt], label=rcp_titles[expt], linewidth=1.5)
    xlim([time[0], time[-1]])
    ylim([-14, 14])
    grid(True)
    xlabel('day of year', fontsize=18)
    title('c) Sea ice to ocean freshwater flux\n2091-2100 minus 1996-2005', fontsize=20)
    # Main title
    suptitle('Amundsen Sea', fontsize=30)
    # Add a legend at the bottom
    ax.legend(bbox_to_anchor=(0.8, -0.1), ncol=4, fontsize=16)
    fig.show()
    fig.savefig('amundsen_budget.png')        

    # One plot with absolute budget
    #fig = figure(figsize=(18,7))
    #gs = GridSpec(1,3)
    #gs.update(left=0.05, right=0.95, bottom=0.2, top=0.88)
    # Total freshwater flux
    #ax = subplot(gs[0,0])
    # Start with a horizontal line at zero
    #ax.plot(time, zeros(size(time)), color=(0.6, 0.6, 0.6), linewidth=2)
    # One line for 1996-2005
    #ax.plot(time, wnet_beg*1e8, color=beg_colour, linewidth=1.5)
    # One line for each RCP
    #for expt in range(num_expts):
        #ax.plot(time, wnet_end[expt,:]*1e8, color=rcp_colours[expt], linewidth=1.5)
    # Configure plot
    #xlim([time[0], time[-1]])
    #ylim([-16, 15])
    #grid(True)
    #xlabel('day of year', fontsize=14)
    #ylabel(r'10$^{-8} $m/s', fontsize=14)
    #title('Net freshwater flux', fontsize=18)
    # Precipitation minus evaporation
    #ax = subplot(gs[0,1])
    #ax.plot(time, zeros(size(time)), color=(0.6, 0.6, 0.6), linewidth=2)
    #ax.plot(time, pminuse_beg*1e8, color=beg_colour, linewidth=1.5)
    #for expt in range(num_expts):
        #ax.plot(time, pminuse_end[expt,:]*1e8, color=rcp_colours[expt], linewidth=1.5)
    #xlim([time[0], time[-1]])
    #ylim([-16, 15])
    #grid(True)
    #xlabel('day of year', fontsize=14)
    #ylabel(r'10$^{-8} $m/s', fontsize=14)
    #title('Precipitation minus evaporation', fontsize=18)
    # Sea ice to ocean freshwater flux
    #ax = subplot(gs[0,2])
    #ax.plot(time, zeros(size(time)), color=(0.6, 0.6, 0.6), linewidth=2)
    #ax.plot(time, ice2ocn_beg*1e8, color=beg_colour, label=beg_title, linewidth=1.5)
    #for expt in range(num_expts):
        #ax.plot(time, ice2ocn_end[expt,:]*1e8, color=rcp_colours[expt], label=rcp_titles[expt], linewidth=1.5)
    #xlim([time[0], time[-1]])
    #ylim([-16, 15])
    #grid(True)
    #xlabel('day of year', fontsize=14)
    #ylabel(r'10$^{-8} $m/s', fontsize=14)
    #title('Sea ice to ocean freshwater flux', fontsize=18)
    # Main title
    #suptitle('Amundsen Sea, 2091-2100 climatologies', fontsize=24)
    # Add a legend at the bottom
    #ax.legend(bbox_to_anchor=(-0.1, -0.1), ncol=3, fontsize=14)
    #fig.show()
    #fig.savefig('amundsen_budget.png')

    # One plot with anomalies in budget
    #fig = figure(figsize=(18,7))
    #gs = GridSpec(1,3)
    #gs.update(left=0.05, right=0.95, bottom=0.2, top=0.88)
    # Total freshwater flux
    #ax = subplot(gs[0,0])
    # Start with a horizontal line at zero
    #ax.plot(time, zeros(size(time)), color=(0.6, 0.6, 0.6), linewidth=2)
    # One line for each RCP
    #for expt in range(num_expts):
        #ax.plot(time, (wnet_end[expt,:]-wnet_beg)*1e8, color=rcp_colours[expt], linewidth=1.5)
    # Configure plot
    #xlim([time[0], time[-1]])
    #ylim([-10, 14])
    #grid(True)
    #xlabel('day of year', fontsize=14)
    #ylabel(r'10$^{-8} $m/s', fontsize=14)
    #title('Net freshwater flux', fontsize=18)
    # Precipitation minus evaporation
    #ax = subplot(gs[0,1])
    #ax.plot(time, zeros(size(time)), color=(0.6, 0.6, 0.6), linewidth=2)
    #for expt in range(num_expts):
        #ax.plot(time, (pminuse_end[expt,:]-pminuse_beg)*1e8, color=rcp_colours[expt], linewidth=1.5)
    #xlim([time[0], time[-1]])
    #ylim([-10, 14])
    #grid(True)
    #xlabel('day of year', fontsize=14)
    #ylabel(r'10$^{-8} $m/s', fontsize=14)
    #title('Precipitation minus evaporation', fontsize=18)
    # Sea ice to ocean freshwater flux
    #ax = subplot(gs[0,2])
    #ax.plot(time, zeros(size(time)), color=(0.6, 0.6, 0.6), linewidth=2)
    #for expt in range(num_expts):
        #ax.plot(time, (ice2ocn_end[expt,:]-ice2ocn_beg)*1e8, color=rcp_colours[expt], label=rcp_titles[expt], linewidth=1.5)
    #xlim([time[0], time[-1]])
    #ylim([-10, 14])
    #grid(True)
    #xlabel('day of year', fontsize=14)
    #ylabel(r'10$^{-8} $m/s', fontsize=14)
    #title('Sea ice to ocean freshwater flux', fontsize=18)
    # Main title
    #suptitle('Amundsen Sea, 2091-2100 minus 1996-2005', fontsize=24)
    # Add a legend at the bottom
    #ax.legend(bbox_to_anchor=(-0.3, -0.1), ncol=2, fontsize=14)
    #fig.show()
    #fig.savefig('amundsen_budget_anom.png')


# Command-line interface
if __name__ == "__main__":

    rcp_amundsen_budget()
    
    
    
