from numpy import *
from matplotlib.pyplot import *

def rcp_9pt_timeseries ():

    # Paths to RCP experiment directories
    directory_head = '/short/y99/kaa561/FESOM/'
    rcp_expt = ['rcp45_M_highres', 'rcp45_A_highres', 'rcp85_M_highres', 'rcp85_A_highres']
    num_rcps = len(rcp_expt)
    # Titles for plot
    rcp_titles = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A']
    # Colours for plotting
    rcp_colours = [(0, 0.4, 1), (0.6, 0.2, 1), (0, 0.6, 0), (1, 0, 0.4)]
    # Paths to control experiment directory
    control_expt = 'highres_spinup'
    # Titles for plot
    control_title = 'CONTROL'
    # Colours for plotting
    control_colour = (0,0,0)
    # Years to plot
    year_start = 1992
    year_end = 2100
    # Spinup years to discard (first 2 repetitions of 1992-2005)
    control_skipyears = 28
    # Output steps per year
    peryear = 365//5
    # Number of ice shelves in timeseries_massloss.py
    num_shelves_all = 24
    # Indices of ice shelves to plot from timeseries_massloss.py
    shelf_indices = [4, 21, 16, 13, 23, 6, 3, 2, 0]
    # Titles for each ice shelf
    shelf_names = ['a) Filchner-Ronne', 'b) Fimbul & Jelbart & Ekstrom', 'c) Amery', 'd) Totten & Moscow University', 'e) Ross', 'f) Pine Island Glacier', 'g) Wilkins & George VI & Stange', 'h) Larsen C', 'i) Total Antarctica']
    num_shelves_plot = len(shelf_names)
    # Bounds on each plot
    min_bounds = [85, 30, 70, 4, 95, 4, 35, 35, 750]
    max_bounds = [315, 135, 200, 25, 265, 44, 115, 160, 1910]

    # Build time axis: 1 point per year
    num_years = year_end - year_start + 1
    time = arange(year_start, year_end+1)

    # Read timeseries for each experiment
    massloss = empty([num_rcps+1, num_shelves_all, num_years])
    # Loop over RCP experiments
    for expt in range(num_rcps):
        # Read logfile
        f = open(directory_head + rcp_expt[expt] + '/massloss.log')
        f.readline()
        # Loop over ice shelves
        for index in range(num_shelves_all):
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
    # Process control experiment
    # Read logfile
    f = open(directory_head + control_expt + '/massloss.log')
    f.readline()
    for index in range(num_shelves_all):
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
            massloss[-1,index,year] = mean(array(massloss_tmp[peryear*year:peryear*(year+1)]))
    f.close()

    # Set up plot
    fig = figure(figsize=(14,10))
    for index in range(num_shelves_plot):
        ax = fig.add_subplot(3,3,index+1)
        # One line for each RCP
        for expt in range(num_rcps):
            ax.plot(time, massloss[expt,shelf_indices[index],:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=1.3)
        # One line for control experiment
        ax.plot(time, massloss[-1,shelf_indices[index],:], label=control_title, color=control_colour, linewidth=1.3)
        # Configure plot
        xlim([amin(time), amax(time)])
        ylim([min_bounds[index], max_bounds[index]])
        grid(True)
        if index < num_shelves_plot-3:
            ax.set_xticklabels([])
        title(shelf_names[index], fontsize=17)
    suptitle('Ice shelf basal mass loss (Gt/y), annual means', fontsize=24)
    subplots_adjust(wspace=0.15, hspace=0.2)
    # Add legend at bottom
    ax.legend(bbox_to_anchor=(0.2,-0.1), ncol=3, fontsize=14)
    fig.show()
    fig.savefig('9pt_timeseries.png')


# Command-line interface
if __name__ == "__main__":

    rcp_9pt_timeseries()

    
    
