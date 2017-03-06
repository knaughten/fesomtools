from numpy import *
from matplotlib.pyplot import *

def timeseries_rep3_compare_annual ():

    directory_head = '/short/y99/kaa561/FESOM/'
    expt_dir = ['lowres_spinup/', 'highres_spinup/']
    expt_titles = ['low res', 'high res']
    expt_colours = ['blue', 'green']
    year_start = 1992
    year_end = 2005
    skipyears = 28
    peryear = 365/5

    dpt_low = 134
    dpt_high = 173.3

    names = ['All Ice Shelves', 'Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    fig_names = ['total_massloss.png', 'larsen_d.png', 'larsen_c.png', 'wilkins_georgevi_stange.png', 'ronne_filchner.png', 'abbot.png', 'pig.png', 'thwaites.png', 'dotson.png', 'getz.png', 'nickerson.png', 'sulzberger.png', 'mertz.png', 'totten_moscowuni.png', 'shackleton.png', 'west.png', 'amery.png', 'princeharald.png', 'baudouin_borchgrevink.png', 'lazarev.png', 'nivl.png', 'fimbul_jelbart_ekstrom.png', 'brunt_riiserlarsen.png', 'ross.png']
    obs_massloss = [1325, 1.4, 20.7, 135.4, 155.4, 51.8, 101.2, 97.5, 45.2, 144.9, 4.2, 18.2, 7.9, 90.6, 72.6, 27.2, 35.5, -2, 21.6, 6.3, 3.9, 26.8, 9.7, 47.7]
    obs_massloss_error = [235, 14, 67, 40, 45, 19, 8, 7, 4, 14, 2, 3, 3, 8, 15, 10, 23, 3, 18, 2, 2, 14, 16, 34]

    time = arange(year_start, year_end+1)
    num_years = year_end - year_start + 1

    dpt = empty([2, num_years])
    for expt in range(2):
        dpt_tmp = []
        f = open(directory_head + expt_dir[expt] + 'dpt.log')
        f.readline()
        for line in f:
            dpt_tmp.append(float(line))
        f.close()
        dpt_tmp = dpt_tmp[skipyears*peryear:(skipyears+num_years)*peryear]
        for year in range(num_years):
            dpt[expt,year] = mean(array(dpt_tmp[peryear*year:peryear*(year+1)]))

    fig, ax = subplots(figsize=(10,6))
    for expt in range(2):
        ax.plot(time, dpt[expt,:], label=expt_titles[expt], color=expt_colours[expt], linewidth=2)
    ax.axhline(dpt_low, color='red', linestyle='dashed', linewidth=2, label='observations')
    ax.axhline(dpt_high, color='red', linewidth=2, linestyle='dashed')
    title('Drake Passage Transport (annually averaged)', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel('Sv', fontsize=14)
    xlim([year_start, year_end])
    grid(True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('drakepsgtrans.png')

    area_min = empty([2, num_years])
    area_max = empty([2, num_years])
    volume_min = empty([2, num_years])
    volume_max = empty([2, num_years])
    for expt in range(2):
        area_tmp = []
        volume_tmp = []
        f = open(directory_head + expt_dir[expt] + 'seaice.log')
        f.readline()
        for line in f:
            try:
                area_tmp.append(float(line))
            except(ValueError):
                break
        for line in f:
            volume_tmp.append(float(line))
        f.close()
        area_tmp = area_tmp[skipyears*peryear:(skipyears+num_years)*peryear]
        volume_tmp = volume_tmp[skipyears*peryear:(skipyears+num_years)*peryear]
        for year in range(num_years):
            area_min[expt,year] = amin(array(area_tmp[peryear*year:peryear*(year+1)]))
            area_max[expt,year] = amax(array(area_tmp[peryear*year:peryear*(year+1)]))
            volume_min[expt,year] = amin(array(volume_tmp[peryear*year:peryear*(year+1)]))
            volume_max[expt,year] = amax(array(volume_tmp[peryear*year:peryear*(year+1)]))

    fig, ax = subplots(figsize=(10,6))
    for expt in range(2):
        ax.plot(time, area_min[expt,:], label=expt_titles[expt], color=expt_colours[expt], linewidth=2)
    title('Annual minimum sea ice area', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel(r'million km$^2$', fontsize=14)
    xlim([year_start, year_end])
    grid(True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('seaice_area_min.png')

    fig, ax = subplots(figsize=(10,6))
    for expt in range(2):
        ax.plot(time, area_max[expt,:], label=expt_titles[expt], color=expt_colours[expt], linewidth=2)
    title('Annual maximum sea ice area', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel(r'million km$^2$', fontsize=14)
    xlim([year_start, year_end])
    grid(True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('seaice_area_max.png')

    fig, ax = subplots(figsize=(10,6))
    for expt in range(2):
        ax.plot(time, volume_min[expt,:], label=expt_titles[expt], color=expt_colours[expt], linewidth=2)
    title('Annual minimum sea ice volume', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel(r'million km$^3$', fontsize=14)
    xlim([year_start, year_end])
    grid(True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('seaice_volume_min.png')

    fig, ax = subplots(figsize=(10,6))
    for expt in range(2):
        ax.plot(time, volume_max[expt,:], label=expt_titles[expt], color=expt_colours[expt], linewidth=2)
    title('Annual maximum sea ice volume', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel(r'million km$^3$', fontsize=14)
    xlim([year_start, year_end])
    grid(True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('seaice_volume_max.png')        

    massloss = empty([2, len(names), num_years])
    for expt in range(2):
        f = open(directory_head + expt_dir[expt] + 'massloss.log')
        f.readline()
        for index in range(len(names)):
            massloss_tmp = []
            for line in f:
                try:
                    massloss_tmp.append(float(line))
                except(ValueError):
                    break
            massloss_tmp = massloss_tmp[skipyears*peryear:(skipyears+num_years)*peryear]
            for year in range(num_years):
                massloss[expt,index,year] = mean(array(massloss_tmp[peryear*year:peryear*(year+1)]))
        f.close()

    for index in range(len(names)):
        massloss_low = obs_massloss[index] - obs_massloss_error[index]
        massloss_high = obs_massloss[index] + obs_massloss_error[index]
        fig, ax = subplots(figsize=(10,6))
        for expt in range(2):
            ax.plot(time, massloss[expt,index,:], label=expt_titles[expt], color=expt_colours[expt], linewidth=2)
        ax.axhline(massloss_low, color='red', linestyle='dashed', linewidth=2, label='observations')
        ax.axhline(massloss_high, color='red', linewidth=2, linestyle='dashed')
        title(names[index] + '\nBasal Mass Loss (annually averaged)', fontsize=18)
        xlabel('Year', fontsize=14)
        ylabel('Gt/y', fontsize=14)
        xlim([year_start, year_end])
        grid(True)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
        fig.savefig(fig_names[index])


if __name__ == "__main__":

    timeseries_rep3_compare_annual()
        

    
    
