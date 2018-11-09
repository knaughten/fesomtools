from numpy import *
from matplotlib.pyplot import *

def watermass_mp_9pt (rcp, model, fig_name):

    # Paths to logfiles
    directory_head = '/short/y99/kaa561/FESOM/'
    beg_log = directory_head + 'highres_spinup/rep3/water_masses_mp.log'
    rcp_log = directory_head + 'rcp' + rcp + '_' + model + '/water_masses_mp.log'
    # Years to consider
    start_year_control = 1992
    end_year_control = 2005
    start_year_rcp = 2006
    end_year_rcp = 2100
    # Sectors
    sector_names = ['Filchner-Ronne Ice Shelf', 'Eastern Weddell Region', 'Amery Ice Shelf', 'Australian Sector', 'Ross Sea', 'Amundsen Sea', 'Bellingshausen Sea', 'Larsen Ice Shelves', 'All Ice Shelves']
    num_sectors = len(sector_names)
    # Water masses
    wm_names = ['ISW', 'HSSW', 'LSSW', 'AASW', 'MCDW', 'CDW']
    num_watermasses = len(wm_names)
    wm_colours = ['cyan', 'black', 'blue', 'green', 'magenta', 'red']
    wm_linestyles = ['-', '-', '--', '-', '--', '-']
    # Title for model
    if model == 'M':
        model_title = 'MMM'
    elif model == 'A':
        model_title = 'ACCESS'

    num_years_control = end_year_control - start_year_control + 1
    num_years_rcp = end_year_rcp - start_year_rcp + 1

    # Build time axis: 1 point per year
    time = arange(start_year_control, end_year_rcp+1)
    # Set up array to hold initial timeseries
    mp_watermass = empty([num_watermasses, num_sectors, num_years_control+num_years_rcp])
    # Read present-day timeseries (rep3, 1992-2005)
    f = open(beg_log, 'r')
    f.readline()
    wm_key = 0
    while wm_key < num_watermasses:
        sector = 0
        while sector < num_sectors:
            year = 0
            for line in f:
                try:
                    mp_watermass[wm_key, sector, year] = float(line)
                    year += 1
                except(ValueError):
                    break
            sector += 1
        wm_key += 1
    f.close()
    # Read RCP timeseries (2006-2100)
    f = open(rcp_log, 'r')
    f.readline()
    wm_key = 0
    while wm_key < num_watermasses:
        sector = 0
        while sector < num_sectors:
            year = 0
            for line in f:
                try:
                    mp_watermass[wm_key, sector, num_years_control+year] = float(line)
                    year += 1
                except(ValueError):
                    break
            sector += 1
        wm_key += 1
    f.close()
    
    # Calculate total melt potential of each region, 1992-2005 average
    mp_watermass_beg = mean(mp_watermass[:,:,0:num_years_control], axis=2)
    mp_total_beg = sum(mp_watermass_beg, axis=0)
    # Scale timeseries as percent of this initial total
    mp_watermass_percent = zeros(shape(mp_watermass))
    for sector in range(num_sectors):
        mp_watermass_percent[:,sector,:] = mp_watermass[:,sector,:]/mp_total_beg[sector]*100

    # Print change in total melt potential, 2091-2100 versus 1992-2005
    mp_watermass_end = mean(mp_watermass[:,:,-10:], axis=2)
    mp_total_end = sum(mp_watermass_end, axis=0)
    for sector in range(num_sectors):
        print 'Change in total melt potential for ' + sector_names[sector] + ': ' + str((mp_total_end[sector]-mp_total_beg[sector])/mp_total_beg[sector]*100)

    # Set up plot
    fig = figure(figsize=(12,10))
    for sector in range(num_sectors):
        ax = fig.add_subplot(3,3,sector+1)
        for wm_key in range(num_watermasses):
            plot(time, mp_watermass_percent[wm_key, sector, :], color=wm_colours[wm_key], label=wm_names[wm_key], linewidth=2, linestyle=wm_linestyles[wm_key])
        xlim([start_year_control, end_year_rcp])
        ylim([0, 165])
        grid(True)
        axvspan(start_year_control, start_year_rcp, facecolor='b', alpha=0.1)
        if sector < num_sectors-3:
            ax.set_xticklabels([])
        else:
            xlabel('Year', fontsize=14)
        if sector % 3 != 0:
            ax.set_yticklabels([])
        if sector == 3:
            ylabel('% of initial heat content in cavity', fontsize=20)
        title(sector_names[sector], fontsize=17)
    suptitle('Heat content (above freezing point)\nof water masses in ice shelf cavities: RCP ' + rcp[0] + '.' + rcp[1] + ' ' + model_title, fontsize=24)
    subplots_adjust(wspace=0.1, hspace=0.2, top=0.87)
    # Add legend at bottom
    ax.legend(bbox_to_anchor=(0.85,-0.2), ncol=6, fontsize=14)
    fig.show()
    fig.savefig(fig_name)


# Command-line interface
if __name__ == "__main__":

    control = False
    key = int(raw_input('RCP 4.5 (4) or 8.5 (8)? '))
    if key == 4:
        rcp = '45'
    elif key == 8:
        rcp = '85'
    model = raw_input('Multi-model mean (M) or ACCESS 1.0 (A)? ')
    fig_name = raw_input('Filename for figure: ')
    watermass_mp_9pt(rcp, model, fig_name)
        
    
