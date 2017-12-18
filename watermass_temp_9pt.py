from numpy import *
from matplotlib.pyplot import *

def watermass_temp_9pt (rcp, model, fig_name):

    # Paths to logfiles
    directory_head = '/short/y99/kaa561/FESOM/'
    beg_log = directory_head + 'highres_spinup/rep3/water_masses_temp_salt.log'
    rcp_log = directory_head + 'rcp' + rcp + '_' + model + '/water_masses_temp_salt.log'
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
    # Title for model
    if model == 'M':
        model_title = 'MMM'
    elif model == 'A':
        model_title = 'ACCESS'

    num_years_control = end_year_control - start_year_control + 1
    num_years_rcp = end_year_rcp - start_year_rcp + 1

    # Build time axis: 1 point per year
    time = arange(start_year_control, end_year_rcp+1)
    # Set up arrays to hold timeseries
    temp_watermass = empty([num_watermasses, num_sectors, num_years_control + num_years_rcp])
    salt_watermass = empty([num_watermasses, num_sectors, num_years_control + num_years_rcp])
    # Read present-day timeseries (rep3, 1992-2005)
    f = open(beg_log, 'r')
    f.readline()
    # Temperature
    wm_key = 0
    while wm_key < num_watermasses:
        sector = 0
        while sector < num_sectors:
            year = 0
            for line in f:
                try:
                    temp_watermass[wm_key, sector, year] = float(line)
                    year += 1
                except(ValueError):
                    break
            sector += 1
        wm_key += 1
    # Salinity
    wm_key = 0
    while wm_key < num_watermasses:
        sector = 0
        while sector < num_sectors:
            year = 0
            for line in f:
                try:
                    salt_watermass[wm_key, sector, year] = float(line)
                    year += 1
                except(ValueError):
                    break
            sector += 1
        wm_key += 1
    f.close()
    # Read RCP timeseries (2006-2100)
    f = open(rcp_log, 'r')
    f.readline()
    # Temperature
    wm_key = 0
    while wm_key < num_watermasses:
        sector = 0
        while sector < num_sectors:
            year = 0
            for line in f:
                try:
                    temp_watermass[wm_key, sector, num_years_control+year] = float(line)
                    year += 1
                except(ValueError):
                    break
            sector += 1
        wm_key += 1
    # Salinity
    wm_key = 0
    while wm_key < num_watermasses:
        sector = 0
        while sector < num_sectors:
            year = 0
            for line in f:
                try:
                    salt_watermass[wm_key, sector, num_years_control+year] = float(line)
                    year += 1
                except(ValueError):
                    break
            sector += 1
        wm_key += 1
    f.close()

    # Set up plot
    fig = figure(figsize=(12,10))
    for sector in range(num_sectors):
        ax = fig.add_subplot(3,3,sector+1)
        for wm_key in range(num_watermasses):
            plot(time, temp_watermass[wm_key, sector, :], color=wm_colours[wm_key], label=wm_names[wm_key], linewidth=2)        
        xlim([start_year_control, end_year_rcp])
        ylim([-2.1, 0.23])
        grid(True)
        axvspan(start_year_control, start_year_rcp, facecolor='b', alpha=0.1)
        if sector < num_sectors-3:
            ax.set_xticklabels([])
        else:
            xlabel('Year', fontsize=14)
        if sector % 3 != 0:
            ax.set_yticklabels([])
        else:
            ylabel(r'$^{\circ}$C', fontsize=14)
        title(sector_names[sector], fontsize=17)
    suptitle('Average water mass temperatures in ice shelf cavities: RCP ' + rcp[0] + '.' + rcp[1] + ' ' + model_title, fontsize=24)
    subplots_adjust(wspace=0.1, hspace=0.2)
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
    watermass_temp_9pt(rcp, model, fig_name)
        
    
