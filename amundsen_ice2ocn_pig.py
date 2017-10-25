from numpy import *
from matplotlib.pyplot import *

def amundsen_ice2ocn_pig (smoothing_years, fig_name):

    # Paths to RCP experiment directories
    directory_head = '/short/y99/kaa561/FESOM/'
    rcp_expt = ['rcp45_M/', 'rcp45_A/', 'rcp85_M/', 'rcp85_A/']
    num_rcps = len(rcp_expt)
    # Titles for plot
    rcp_titles = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A']
    # Colours for plotting
    rcp_colours = ['blue', 'cyan', 'green', 'magenta']
    # Path to control experiment directory
    control_expt = 'highres_spinup/'
    # Title for plot
    control_title = 'CONTROL'
    # Colour for plotting
    control_colour = 'black'
    # Logfile names
    amundsen_log = 'amundsen.log'
    massloss_log = 'massloss.log'
    # Index of PIG in timeseries_massloss.py
    pig_index = 6
    # Years to plot
    year_start = 1992
    year_end = 2100
    # Spinup years to discard (first 2 repetitions of 1992-2005)
    control_skipyears = 28
    # Output steps per year
    peryear = 365/5

    # Figure out window size on smoothing radius; treat even and odd #bins
    # differently
    if smoothing_years % 2 == 0:
        even = True
        window = smoothing_years/2
    else:
        even = False
        window = (smoothing_years-1)/2

    # Build time axis: 1 point per year
    num_years = year_end - year_start + 1
    time_full = arange(year_start, year_end+1)

    # Amundsen Sea ice-to-ocean freshwater flux
    ice2ocn_full = empty([num_rcps+1, num_years])
    # Loop over RCP experiments
    for expt in range(num_rcps):
        # Read logfile
        ice2ocn_tmp = []
        f = open(directory_head + rcp_expt[expt] + amundsen_log, 'r')
        f.readline()
        for line in f:
            ice2ocn_tmp.append(float(line))
        f.close()
        # Calculate annual averages
        for year in range(num_years):
            ice2ocn_full[expt,year] = mean(array(ice2ocn_tmp[peryear*year:peryear*(year+1)]))
    # Control experiment
    ice2ocn_tmp = []
    f = open(directory_head + control_expt + amundsen_log, 'r')
    f.readline()
    for line in f:
        ice2ocn_tmp.append(float(line))
    f.close()
    # Throw away first 2 repetitions
    ice2ocn_tmp = ice2ocn_tmp[control_skipyears*peryear:]
    for year in range(num_years):
        ice2ocn_full[-1,year] = mean(array(ice2ocn_tmp[peryear*year:peryear*(year+1)]))    

    # PIG basal mass loss
    pig_massloss_full = empty([num_rcps+1, num_years])
    for expt in range(num_rcps):
        f = open(directory_head + rcp_expt[expt] + massloss_log, 'r')
        f.readline()
        # Read the ice shelves before PIG
        for index in range(pig_index):
            for line in f:
                try:
                    tmp = float(line)
                except(ValueError):
                    break
        # Now read PIG
        pig_massloss_tmp = []
        for line in f:
            try:
                pig_massloss_tmp.append(float(line))
            except(ValueError):
                break
        f.close()
        for year in range(num_years):
            pig_massloss_full[expt,year] = mean(array(pig_massloss_tmp[peryear*year:peryear*(year+1)]))
    f = open(directory_head + control_expt + massloss_log, 'r')
    f.readline()
    for index in range(pig_index):
        for line in f:
            try:
                tmp = float(line)
            except(ValueError):
                break
    pig_massloss_tmp = []
    for line in f:
        try:
            pig_massloss_tmp.append(float(line))
        except(ValueError):
            break
    f.close()
    pig_massloss_tmp = pig_massloss_tmp[control_skipyears*peryear:]
    for year in range(num_years):
        pig_massloss_full[-1,year] = mean(array(pig_massloss_tmp[peryear*year:peryear*(year+1)]))

    # Running mean
    if even:
        num_years_central = num_years - 2*window + 1
        time = empty([num_years_central])
        ice2ocn = empty([num_rcps+1,num_years_central])
        pig_massloss = empty([num_rcps+1,num_years_central])
        t = 0
        # Moving average between the existing timesteps, considering
        # <window> points on either side
        for year in range(window-1, num_years-window):
            time[t] = time_full[year]+0.5
            ice2ocn[:,t] = mean(ice2ocn_full[:,year-window+1:year+window+1], axis=1)
            pig_massloss[:,t] = mean(pig_massloss_full[:,year-window+1:year+window+1], axis=1)
            t += 1
    else:
        num_years_central = num_years - 2*window
        time = empty([num_years_central])
        ice2ocn = empty([num_rcps+1,num_years_central])
        pig_massloss = empty([num_rcps+1,num_years_central])
        t = 0
        # Moving average at the existing timesteps, considering central
        # point as well as <window> points on either side
        for year in range(window, num_years-window):
            time[t] = time_full[year]
            ice2ocn[:,t] = mean(ice2ocn_full[:,year-window:year+window+1], axis=1)
            pig_massloss[:,t] = mean(pig_massloss_full[:,year-window:year+window+1], axis=1)
            t += 1

    # Plot
    fig = figure(figsize=(12,8))
    gs = GridSpec(2,1)
    gs.update(left=0.1, right=0.95, bottom=0.15, top=0.95, hspace=0.02)
    # Ice to ocean freshwater flux
    ax = subplot(gs[0,0])
    # Horizontal line at "threshold"
    #ax.axhline(-4, color=(0.6,0.6,0.6), linewidth=2, linestyle='dashed')
    # One line for each RCP
    for expt in range(num_rcps):
        ax.plot(time, ice2ocn[expt,:], color=rcp_colours[expt], linewidth=1.5)
    # One line for control experiment
    ax.plot(time, ice2ocn[-1,:], color=control_colour, linewidth=1.5)
    xlim([amin(time), amax(time)])
    ax.set_xticks(arange(1995, amax(time), 5))
    ax.set_xticklabels([])
    grid(True)    
    ylabel(r'10$^{-8}$ m/s', fontsize=18)
    text(0.02, 0.9, 'Amundsen Sea ice-to-ocean freshwater flux ('+str(smoothing_years)+'-year running mean)', horizontalalignment='left', transform=ax.transAxes, fontsize=18)
    # PIG mass loss
    ax = subplot(gs[1,0])
    # Horizontal line at "threshold"
    #ax.axhline(15, color=(0.6,0.6,0.6), linewidth=2, linestyle='dashed')
    for expt in range(num_rcps):
        ax.plot(time, pig_massloss[expt,:], color=rcp_colours[expt], label=rcp_titles[expt], linewidth=1.5)
    ax.plot(time, pig_massloss[-1,:], color=control_colour, label=control_title, linewidth=1.5)
    xlim([amin(time), amax(time)])
    ax.set_xticks(arange(1995, amax(time), 5))
    labels = []
    for year in arange(1995, amax(time), 5):
        if year % 10 == 0:
            labels.append(str(int(year)))
        else:
            labels.append('')
    ax.set_xticklabels(labels)
    grid(True)
    xlabel('Year',fontsize=18)
    ylabel('Gt/y', fontsize=18)
    text(0.02, 0.9, 'Pine Island Ice Shelf basal mass loss ('+str(smoothing_years)+'-year running mean)', horizontalalignment='left', transform=ax.transAxes, fontsize=18)
    # Legend at bottom
    ax.legend(bbox_to_anchor=(0.95, -0.15), ncol=5, fontsize=14)
    fig.show()
    fig.savefig(fig_name)


# Command-line interface
if __name__ == "__main__":

    smoothing_years = int(raw_input("Number of years for running mean: "))
    fig_name = raw_input("Filename for figure: ")
    amundsen_ice2ocn_pig(smoothing_years, fig_name)
    
                

    
