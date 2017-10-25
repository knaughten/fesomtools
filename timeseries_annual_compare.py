from numpy import *
from matplotlib.pyplot import *

# Compare all high-res RCP experiments and control by plotting timeseries on
# the same axes: annual averages of Drake Passage transport and ice shelf mass
# loss, annual max/mins of sea ice area, extent, and volume.
def timeseries_annual_compare ():

    # Paths to RCP experiment directories
    directory_head = '/short/y99/kaa561/FESOM/'
    rcp_expt = ['rcp45_M', 'rcp45_A', 'rcp85_M', 'rcp85_A']
    # Titles for plot
    rcp_titles = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A']
    # Colours for plotting
    rcp_colours = ['blue', 'cyan', 'green', 'magenta']
    # Paths to control experiment directories
    control_expt = ['highres_spinup']
    # Titles for plot
    control_titles = ['CONTROL']
    # Colours for plotting
    control_colours = ['black']
    # Years to plot
    year_start = 1992
    year_end = 2100
    # Spinup years to discard (first 2 repetitions of 1992-2005)
    control_skipyears = 28
    # Output steps per year
    peryear = 365/5

    # Names of ice shelves for plotting
    names = ['All Ice Shelves', 'Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    # Filenames for ice shelf plots
    fig_names = ['total_massloss.png', 'larsen_d.png', 'larsen_c.png', 'wilkins_georgevi_stange.png', 'ronne_filchner.png', 'abbot.png', 'pig.png', 'thwaites.png', 'dotson.png', 'getz.png', 'nickerson.png', 'sulzberger.png', 'mertz.png', 'totten_moscowuni.png', 'shackleton.png', 'west.png', 'amery.png', 'princeharald.png', 'baudouin_borchgrevink.png', 'lazarev.png', 'nivl.png', 'fimbul_jelbart_ekstrom.png', 'brunt_riiserlarsen.png', 'ross.png']

    # Build time axis: 1 point per year
    num_years = year_end - year_start + 1
    time = arange(year_start, year_end+1)

    # Drake Passage transport, annual averages
    dpt = empty([len(rcp_expt)+len(control_expt), num_years])
    # Loop over RCP experiments
    for expt in range(len(rcp_expt)):
        # Read logfile
        dpt_tmp = []
        f = open(directory_head + rcp_expt[expt] + '/dpt.log')
        f.readline()
        for line in f:
            dpt_tmp.append(float(line))
        f.close()
        # Calculate annual averages
        for year in range(num_years):
            dpt[expt,year] = mean(array(dpt_tmp[peryear*year:peryear*(year+1)]))
    # Loop over control experiments
    for expt in range(len(control_expt)):
        # Read logfile
        dpt_tmp = []
        f = open(directory_head + control_expt[expt] + '/dpt.log')
        f.readline()
        for line in f:
            dpt_tmp.append(float(line))
        f.close()
        # Throw away first 2 repetitions
        dpt_tmp = dpt_tmp[control_skipyears*peryear:]
        # Calculate annual averages
        for year in range(num_years):
            dpt[len(rcp_expt)+expt,year] = mean(array(dpt_tmp[peryear*year:peryear*(year+1)]))

    # Plot
    fig, ax = subplots(figsize=(10,6))
    # One line for each RCP
    for expt in range(len(rcp_expt)):
        ax.plot(time, dpt[expt,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
    # One line for each control experiment
    for expt in range(len(control_expt)):
        ax.plot(time, dpt[len(rcp_expt)+expt,:], label=control_titles[expt], color=control_colours[expt], linewidth=2)
    # Configure plot
    title('Drake Passage Transport (annually averaged)', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel('Sv', fontsize=14)
    xlim([year_start, amax(time)])
    grid(True)
    # Move plot over to make room for legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    # Make legend
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('drakepstrans.png')

    # Sea ice area and volume, annual min/max
    area_min = empty([len(rcp_expt)+len(control_expt), num_years])
    area_max = empty([len(rcp_expt)+len(control_expt), num_years])
    volume_min = empty([len(rcp_expt)+len(control_expt), num_years])
    volume_max = empty([len(rcp_expt)+len(control_expt), num_years])
    # Loop over RCP experiments
    for expt in range(len(rcp_expt)):
        # Read logfile
        area_tmp = []
        volume_tmp = []        
        f = open(directory_head + rcp_expt[expt] + '/seaice.log')
        f.readline()
        for line in f:
            try:
                area_tmp.append(float(line))
            except(ValueError):
                break
        for line in f:
            volume_tmp.append(float(line))
        f.close()
        # Get annual min/max of area and volume
        for year in range(num_years):
            area_min[expt,year] = amin(array(area_tmp[peryear*year:peryear*(year+1)]))
            area_max[expt,year] = amax(array(area_tmp[peryear*year:peryear*(year+1)]))
            volume_min[expt,year] = amin(array(volume_tmp[peryear*year:peryear*(year+1)]))
            volume_max[expt,year] = amax(array(volume_tmp[peryear*year:peryear*(year+1)]))
    # Loop over control experiments
    for expt in range(len(control_expt)):
        # Read logfile
        area_tmp = []
        volume_tmp = []
        f = open(directory_head + control_expt[expt] + '/seaice.log')
        f.readline()
        for line in f:
            try:
                area_tmp.append(float(line))
            except(ValueError):
                break
        for line in f:
            volume_tmp.append(float(line))
        f.close()
        # Throw away spinup years
        area_tmp = area_tmp[control_skipyears*peryear:]
        volume_tmp = volume_tmp[control_skipyears*peryear:]
        # Get annual min/max of area and volume
        for year in range(num_years):
            area_min[len(rcp_expt)+expt,year] = amin(array(area_tmp[peryear*year:peryear*(year+1)]))
            area_max[len(rcp_expt)+expt,year] = amax(array(area_tmp[peryear*year:peryear*(year+1)]))
            volume_min[len(rcp_expt)+expt,year] = amin(array(volume_tmp[peryear*year:peryear*(year+1)]))
            volume_max[len(rcp_expt)+expt,year] = amax(array(volume_tmp[peryear*year:peryear*(year+1)]))

    # Plot area min             
    fig, ax = subplots(figsize=(10,6))
    # One line for each RCP
    for expt in range(len(rcp_expt)):
        ax.plot(time, area_min[expt,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
    # One line for each control experiment
    for expt in range(len(control_expt)):
        ax.plot(time, area_min[len(rcp_expt)+expt,:], label=control_titles[expt], color=control_colours[expt], linewidth=2)
    # Configure plot
    title('Annual minimum sea ice area', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel(r'million km$^2$', fontsize=14)
    xlim([year_start, amax(time)])
    grid(True)
    # Move the plot over to make room for legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    # Make legend
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('seaice_area_min.png')

    # Repeat for area max
    fig, ax = subplots(figsize=(10,6))
    for expt in range(len(rcp_expt)):
        ax.plot(time, area_max[expt,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
    for expt in range(len(control_expt)):
        ax.plot(time, area_max[len(rcp_expt)+expt,:], label=control_titles[expt], color=control_colours[expt], linewidth=2)
    title('Annual maximum sea ice area', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel(r'million km$^2$', fontsize=14)
    xlim([year_start, amax(time)])
    grid(True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('seaice_area_max.png')

    # Repeat for volume min
    fig, ax = subplots(figsize=(10,6))
    for expt in range(len(rcp_expt)):
        ax.plot(time, volume_min[expt,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
    for expt in range(len(control_expt)):
        ax.plot(time, volume_min[len(rcp_expt)+expt,:], label=control_titles[expt], color=control_colours[expt], linewidth=2)
    title('Annual minimum sea ice volume', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel(r'million km$^3$', fontsize=14)
    xlim([year_start, amax(time)])
    grid(True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('seaice_volume_min.png')

    # Repeat for volume max
    fig, ax = subplots(figsize=(10,6))
    for expt in range(len(rcp_expt)):
        ax.plot(time, volume_max[expt,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
    for expt in range(len(control_expt)):
        ax.plot(time, volume_max[len(rcp_expt)+expt,:], label=control_titles[expt], color=control_colours[expt], linewidth=2)
    title('Annual maximum sea ice volume', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel(r'million km$^3$', fontsize=14)
    xlim([year_start, amax(time)])
    grid(True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('seaice_volume_max.png')

    '''# Sea ice extent, annual min/max
    extent_min = empty([len(rcp_expt)+len(control_expt), num_years])
    extent_max = empty([len(rcp_expt)+len(control_expt), num_years])
    # Loop over RCP experiments
    for expt in range(len(rcp_expt)):
        # Read logfile
        extent_tmp = []
        f = open(directory_head + rcp_expt[expt] + '/seaice_extent.log')
        f.readline()
        for line in f:
            extent_tmp.append(float(line))
        f.close()
        # Get annual min/max of extent
        for year in range(num_years):
            extent_min[expt,year] = amin(array(extent_tmp[peryear*year:peryear*(year+1)]))
            extent_max[expt,year] = amax(array(extent_tmp[peryear*year:peryear*(year+1)]))
    # Loop over control experiments
    for expt in range(len(control_expt)):
        # Read logfile
        extent_tmp = []
        f = open(directory_head + control_expt[expt] + '/seaice_extent.log')
        f.readline()
        for line in f:
            extent_tmp.append(float(line))
        f.close()
        # Throw away spinup years
        extent_tmp = extent_tmp[control_skipyears*peryear:]
        # Get annual min/max of extent
        for year in range(num_years):
            extent_min[len(rcp_expt)+expt,year] = amin(array(extent_tmp[peryear*year:peryear*(year+1)]))
            extent_max[len(rcp_expt)+expt,year] = amax(array(extent_tmp[peryear*year:peryear*(year+1)]))

    # Plot extent min             
    fig, ax = subplots(figsize=(10,6))
    # One line for each RCP
    for expt in range(len(rcp_expt)):
        ax.plot(time, extent_min[expt,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
    # One line for each control experiment
    for expt in range(len(control_expt)):
        ax.plot(time, extent_min[len(rcp_expt)+expt,:], label=control_titles[expt], color=control_colours[expt], linewidth=2)
    # Configure plot
    title('Annual minimum sea ice extent', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel(r'million km$^2$', fontsize=14)
    xlim([year_start, amax(time)])
    grid(True)
    # Move the plot over to make room for legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    # Make legend
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('seaice_extent_min.png')

    # Repeat for extent max
    fig, ax = subplots(figsize=(10,6))
    for expt in range(len(rcp_expt)):
        ax.plot(time, extent_max[expt,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
    for expt in range(len(control_expt)):
        ax.plot(time, extent_max[len(rcp_expt)+expt,:], label=control_titles[expt], color=control_colours[expt], linewidth=2)
    title('Annual maximum sea ice extent', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel(r'million km$^2$', fontsize=14)
    xlim([year_start, amax(time)])
    grid(True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('seaice_extent_max.png')'''

    # Ice shelf basal mass loss
    massloss = empty([len(rcp_expt)+len(control_expt), len(names), num_years])
    # Loop over RCP experiments
    for expt in range(len(rcp_expt)):
        # Read logfile
        f = open(directory_head + rcp_expt[expt] + '/massloss.log')
        f.readline()
        # Loop over ice shelves
        for index in range(len(names)):
            massloss_tmp = []
            for line in f:
                try:
                    massloss_tmp.append(float(line))
                except(ValueError):
                    break
            # Calculate annual averages
            for year in range(num_years):
                massloss[expt,index,year] = mean(array(massloss_tmp[peryear*year:peryear*(year+1)]))
        f.close()
    # Loop over control experiments
    for expt in range(len(control_expt)):
        # Read logfile
        f = open(directory_head + control_expt[expt] + '/massloss.log')
        f.readline()
        # Loop over ice shelves
        for index in range(len(names)):
            massloss_tmp = []
            for line in f:
                try:
                    massloss_tmp.append(float(line))
                except(ValueError):
                    break
            # Throw away spinup years
            massloss_tmp = massloss_tmp[control_skipyears*peryear:]
            # Calculate annual averages
            for year in range(num_years):
                massloss[len(rcp_expt)+expt,index,year] = mean(array(massloss_tmp[peryear*year:peryear*(year+1)]))
        f.close()

    # One plot for each ice shelf
    for index in range(len(names)):
        fig, ax = subplots(figsize=(10, 6))
        # One line for each RCP
        for expt in range(len(rcp_expt)):
            ax.plot(time, massloss[expt,index,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
        # One line for each control experiment
        for expt in range(len(control_expt)):
            ax.plot(time, massloss[len(rcp_expt)+expt,index,:], label=control_titles[expt], color=control_colours[expt], linewidth=2)
        # Configure plot
        xlim([year_start, amax(time)])
        xlabel('Year', fontsize=14)
        ylabel('Gt/y', fontsize=14)
        grid(True)
        title(names[index] + '\nBasal Mass Loss (annually averaged)', fontsize=18)
        # Move plot over to make room for legend
        box = ax.get_position()
        # Make legend
        ax.set_position([box.x0, box.y0, box.width*0.75, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1.1,0.5))
        fig.savefig(fig_names[index])

    
# Command-line interface
if __name__ == "__main__":

    timeseries_annual_compare()
    
        

    
