from numpy import *
from matplotlib.pyplot import *

# Plot Drake Passage transport and ice shelf melt rates/mass loss for the third
# repetition of the spinup forcing (1992-2005), for both the low-res and
# high-res control simulations.
def timeseries_rep3_compare ():

    # Paths to experiment directories
    directory_head = '/short/y99/kaa561/FESOM/'
    expt_dir = ['lowres_spinup/', 'highres_spinup/']
    # Titles for plotting
    expt_titles = ['low res', 'high res']
    # Colours for plotting
    expt_colours = ['blue', 'green']
    year_start = 1992
    year_end = 2005
    # Skip the first 2 repetitions
    skipyears = 28
    # Number of records per year (assumes 5-day averages)
    peryear = 365//5

    # Bounds of observations for Drake Passage transport
    dpt_low = 134
    dpt_high = 173.3

    # Titles for each ice shelf
    names = ['All Ice Shelves', 'Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    # Figure names for each ice shelf
    fig_names = ['total_massloss.png', 'larsen_d.png', 'larsen_c.png', 'wilkins_georgevi_stange.png', 'ronne_filchner.png', 'abbot.png', 'pig.png', 'thwaites.png', 'dotson.png', 'getz.png', 'nickerson.png', 'sulzberger.png', 'mertz.png', 'totten_moscowuni.png', 'shackleton.png', 'west.png', 'amery.png', 'princeharald.png', 'baudouin_borchgrevink.png', 'lazarev.png', 'nivl.png', 'fimbul_jelbart_ekstrom.png', 'brunt_riiserlarsen.png', 'ross.png']
    # Observed mass loss (Rignot 2013) and uncertainty for each ice shelf, in
    # Gt/y
    obs_massloss = [1325, 1.4, 20.7, 135.4, 155.4, 51.8, 101.2, 97.5, 45.2, 144.9, 4.2, 18.2, 7.9, 90.6, 72.6, 27.2, 35.5, -2, 21.6, 6.3, 3.9, 26.8, 9.7, 47.7]
    obs_massloss_error = [235, 14, 67, 40, 45, 19, 8, 7, 4, 14, 2, 3, 3, 8, 15, 10, 23, 3, 18, 2, 2, 14, 16, 34]

    # Make time axis
    time = arange(year_start, year_end+1, 1.0/peryear)    
    num_time = size(time)
    num_years = year_end - year_start + 1

    # Drake Passage transport
    dpt = empty([2, num_time])
    # Loop over experiments
    for expt in range(2):
        # Read logfile
        dpt_tmp = []
        f = open(directory_head + expt_dir[expt] + 'dpt.log')
        f.readline()
        for line in f:
            dpt_tmp.append(float(line))
        f.close()
        # Select third repetition
        dpt_tmp = dpt_tmp[skipyears*peryear:(skipyears+num_years)*peryear]
        dpt[expt,:] = array(dpt_tmp)

    # Plot
    fig, ax = subplots(figsize=(10,6))
    for expt in range(2):
        ax.plot(time, dpt[expt,:], label=expt_titles[expt], color=expt_colours[expt], linewidth=2)
    # Add lines for range of observations
    ax.axhline(dpt_low, color='red', linestyle='dashed', linewidth=2, label='observations')
    ax.axhline(dpt_high, color='red', linewidth=2, linestyle='dashed')
    title('Drake Passage Transport', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel('Sv', fontsize=14)
    xlim([year_start, amax(time)])
    grid(True)
    # Move plot over to make room for legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    # Make legend
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('drakepsgtrans.png')

    # Ice shelf mass loss
    massloss = empty([2, len(names), num_time])
    # Loop over experiments
    for expt in range(2):
        # Read logfile
        f = open(directory_head + expt_dir[expt] + 'massloss.log')
        f.readline()
        # Loop over ice shelves
        for index in range(len(names)):
            massloss_tmp = []
            for line in f:
                try:
                    massloss_tmp.append(float(line))
                except(ValueError):
                    break
            # Select third repetition
            massloss_tmp = massloss_tmp[skipyears*peryear:(skipyears+num_years)*peryear]
            massloss[expt,index,:] = array(massloss_tmp)
        f.close()

    # One plot for each ice shelf
    for index in range(len(names)):
        # Calculate range of observations
        massloss_low = obs_massloss[index] - obs_massloss_error[index]
        massloss_high = obs_massloss[index] + obs_massloss_error[index]
        fig, ax = subplots(figsize=(10,6))
        for expt in range(2):
            ax.plot(time, massloss[expt,index,:], label=expt_titles[expt], color=expt_colours[expt], linewidth=2)
        # Add lines for range of observations
        ax.axhline(massloss_low, color='red', linestyle='dashed', linewidth=2, label='observations')
        ax.axhline(massloss_high, color='red', linewidth=2, linestyle='dashed')
        title(names[index] + '\nBasal Mass Loss', fontsize=18)
        xlabel('Year', fontsize=14)
        ylabel('Gt/y', fontsize=14)
        xlim([year_start, amax(time)])
        grid(True)
        # Move plot over to make room for legend
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        # Make legend
        ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
        fig.savefig(fig_names[index])


# Command-line interface
if __name__ == "__main__":

    timeseries_rep3_compare()
        

    
    
