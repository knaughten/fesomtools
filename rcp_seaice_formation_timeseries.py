from numpy import *
from matplotlib.pyplot import *

def rcp_seaice_formation_timeseries ():

    # Paths to RCP experiment directories
    directory_beg = '/short/y99/kaa561/FESOM/highres_spinup/'
    directories = ['/short/y99/kaa561/FESOM/rcp45_M/', '/short/y99/kaa561/FESOM/rcp45_A/', '/short/y99/kaa561/FESOM/rcp85_M/', '/short/y99/kaa561/FESOM/rcp85_A/']
    num_expts = len(directories)
    # Name of logfile from timeseries_seaice_formation.py
    log_name = 'seaice_formation.log'
    # Titles for plot
    beg_title = 'CONTROL'
    rcp_titles = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A']
    # Colours for plotting
    beg_colour = (0,0,0)
    rcp_colours = [(0, 0.4, 1), (0.6, 0.2, 1), (0, 0.6, 0), (1, 0, 0.4)]
    # Years to plot
    year_start = 1992
    rcp_year_start = 2006
    year_end = 2100

    num_years = year_end - year_start + 1
    num_years_rep3 = rcp_year_start  - year_start
    time = arange(year_start, year_end+1)

    # Set up arrays    
    formation = empty([num_expts+1, num_years])
    # Read control experiment (correct years already selected)
    f = open(directory_beg + log_name, 'r')
    f.readline()
    t = 0
    for line in f:
        formation[-1,t] = float(line)
        t += 1
    f.close()
    # Loop over RCPs
    for expt in range(num_expts):
        # Copy the beginning of the control
        formation[expt, 0:num_years_rep3] = formation[-1, 0:num_years_rep3]
        # Now read logfile for the rest of it
        f = open(directories[expt] + log_name, 'r')
        f.readline()
        t = num_years_rep3
        for line in f:
            formation[expt,t] = float(line)
            t += 1
        f.close()

    beg_mean = mean(formation[-1,1996-year_start:rcp_year_start-year_start])
    print('Mean over 1996-2005: ' + str(beg_mean) + ' thousand km^3/y')
    print('Mean over 2091-2100:')
    for expt in range(num_expts):
        end_mean = mean(formation[expt,-10:])
        print('...' + rcp_titles[expt] + ': ' + str(end_mean) + ', change of ' + str((end_mean-beg_mean)/beg_mean*100))

    # Plot        
    fig, ax = subplots()
    # One line for each RCP
    for expt in range(num_expts):
        ax.plot(time, formation[expt,:], color=rcp_colours[expt], label=rcp_titles[expt], linewidth=1.8)
    # One line for control experiment
    ax.plot(time, formation[-1,:], color=beg_colour, label=beg_title, linewidth=1.8)
    title('Net sea ice formation on continental shelf', fontsize=20)
    xlabel('Year', fontsize=14)
    ylabel(r'10$^3$ km$^3$/y', fontsize=14)
    xlim([year_start, year_end])
    grid(True)
    ax.legend(loc='lower left', fontsize=14)
    fig.show()
    fig.savefig('rcp_seaice_formation.png')


# Command-line interface
if __name__ == "__main__":

    rcp_seaice_formation_timeseries()
    
    
        
    

    
