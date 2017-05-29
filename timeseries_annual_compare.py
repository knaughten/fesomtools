from numpy import *
from matplotlib.pyplot import *

def timeseries_annual_compare ():

    directory_head = '/short/y99/kaa561/FESOM/'
    rcp_expt = ['rcp45_M', 'rcp85_M', 'rcp45_A', 'rcp85_A', 'rcp85_M_highres']
    rcp_titles = ['RCP_45_M', 'RCP_85_M', 'RCP_45_A', 'RCP_85_A', 'RCP_85_M_HR']
    rcp_colours = ['blue', 'red', 'green', 'magenta', 'cyan']
    control_expt = ['lowres_spinup', 'highres_spinup']
    control_titles = ['CONTROL_LR', 'CONTROL_HR']
    control_colours = ['black', (0.5, 0.5, 0.5)]
    year_start = 1992
    year_end = 2100
    control_skipyears = 28
    peryear = 365/5

    names = ['All Ice Shelves', 'Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    fig_names = ['total_massloss.png', 'larsen_d.png', 'larsen_c.png', 'wilkins_georgevi_stange.png', 'ronne_filchner.png', 'abbot.png', 'pig.png', 'thwaites.png', 'dotson.png', 'getz.png', 'nickerson.png', 'sulzberger.png', 'mertz.png', 'totten_moscowuni.png', 'shackleton.png', 'west.png', 'amery.png', 'princeharald.png', 'baudouin_borchgrevink.png', 'lazarev.png', 'nivl.png', 'fimbul_jelbart_ekstrom.png', 'brunt_riiserlarsen.png', 'ross.png']

    num_years = year_end - year_start + 1
    time = arange(year_start, year_end+1)

    dpt = empty([len(rcp_expt)+len(control_expt), num_years]) 
    for expt in range(len(rcp_expt)):
        dpt_tmp = []
        f = open(directory_head + rcp_expt[expt] + '/dpt.log')
        f.readline()
        for line in f:
            dpt_tmp.append(float(line))
        f.close()
        for year in range(num_years):
            dpt[expt,year] = mean(array(dpt_tmp[peryear*year:peryear*(year+1)]))
    for expt in range(len(control_expt)):
        dpt_tmp = []
        f = open(directory_head + control_expt[expt] + '/dpt.log')
        f.readline()
        for line in f:
            dpt_tmp.append(float(line))
        f.close()
        dpt_tmp = dpt_tmp[control_skipyears*peryear:]
        for year in range(num_years):
            dpt[len(rcp_expt)+expt,year] = mean(array(dpt_tmp[peryear*year:peryear*(year+1)]))

    fig, ax = subplots(figsize=(10,6))
    for expt in range(len(rcp_expt)):
        ax.plot(time, dpt[expt,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
    for expt in range(len(control_expt)):
        ax.plot(time, dpt[len(rcp_expt)+expt,:], label=control_titles[expt], color=control_colours[expt], linewidth=2)
    title('Drake Passage Transport (annually averaged)', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel('Sv', fontsize=14)
    xlim([year_start, year_end])
    grid(True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('drakepstrans.png')

    area_min = empty([len(rcp_expt)+len(control_expt), num_years])
    area_max = empty([len(rcp_expt)+len(control_expt), num_years])
    volume_min = empty([len(rcp_expt)+len(control_expt), num_years])
    volume_max = empty([len(rcp_expt)+len(control_expt), num_years])
    for expt in range(len(rcp_expt)):
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
        for year in range(num_years):
            area_min[expt,year] = amin(array(area_tmp[peryear*year:peryear*(year+1)]))
            area_max[expt,year] = amax(array(area_tmp[peryear*year:peryear*(year+1)]))
            volume_min[expt,year] = amin(array(volume_tmp[peryear*year:peryear*(year+1)]))
            volume_max[expt,year] = amax(array(volume_tmp[peryear*year:peryear*(year+1)]))
    for expt in range(len(control_expt)):
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
        area_tmp = area_tmp[control_skipyears*peryear:]
        volume_tmp = volume_tmp[control_skipyears*peryear:]
        for year in range(num_years):
            area_min[len(rcp_expt)+expt,year] = amin(array(area_tmp[peryear*year:peryear*(year+1)]))
            area_max[len(rcp_expt)+expt,year] = amax(array(area_tmp[peryear*year:peryear*(year+1)]))
            volume_min[len(rcp_expt)+expt,year] = amin(array(volume_tmp[peryear*year:peryear*(year+1)]))
            volume_max[len(rcp_expt)+expt,year] = amax(array(volume_tmp[peryear*year:peryear*(year+1)]))

    fig, ax = subplots(figsize=(10,6))
    for expt in range(len(rcp_expt)):
        ax.plot(time, area_min[expt,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
    for expt in range(len(control_expt)):
        ax.plot(time, area_min[len(rcp_expt)+expt,:], label=control_titles[expt], color=control_colours[expt], linewidth=2)
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
    for expt in range(len(rcp_expt)):
        ax.plot(time, area_max[expt,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
    for expt in range(len(control_expt)):
        ax.plot(time, area_max[len(rcp_expt)+expt,:], label=control_titles[expt], color=control_colours[expt], linewidth=2)
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
    for expt in range(len(rcp_expt)):
        ax.plot(time, volume_min[expt,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
    for expt in range(len(control_expt)):
        ax.plot(time, volume_min[len(rcp_expt)+expt,:], label=control_titles[expt], color=control_colours[expt], linewidth=2)
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
    for expt in range(len(rcp_expt)):
        ax.plot(time, volume_max[expt,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
    for expt in range(len(control_expt)):
        ax.plot(time, volume_max[len(rcp_expt)+expt,:], label=control_titles[expt], color=control_colours[expt], linewidth=2)
    title('Annual maximum sea ice volume', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel(r'million km$^3$', fontsize=14)
    xlim([year_start, year_end])
    grid(True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('seaice_volume_max.png')

    massloss = empty([len(rcp_expt)+len(control_expt), len(names), num_years])
    for expt in range(len(rcp_expt)):
        f = open(directory_head + rcp_expt[expt] + '/massloss.log')
        f.readline()
        for index in range(len(names)):
            massloss_tmp = []
            for line in f:
                try:
                    massloss_tmp.append(float(line))
                except(ValueError):
                    break
            for year in range(num_years):
                massloss[expt,index,year] = mean(array(massloss_tmp[peryear*year:peryear*(year+1)]))
        f.close()
    for expt in range(len(control_expt)):
        f = open(directory_head + control_expt[expt] + '/massloss.log')
        f.readline()
        for index in range(len(names)):
            massloss_tmp = []
            for line in f:
                try:
                    massloss_tmp.append(float(line))
                except(ValueError):
                    break
            massloss_tmp = massloss_tmp[control_skipyears*peryear:]
            for year in range(num_years):
                massloss[len(rcp_expt)+expt,index,year] = mean(array(massloss_tmp[peryear*year:peryear*(year+1)]))
        f.close()

    for index in range(len(names)):
        fig, ax = subplots(figsize=(10, 6))
        for expt in range(len(rcp_expt)):
            ax.plot(time, massloss[expt,index,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
        for expt in range(len(control_expt)):
            ax.plot(time, massloss[len(rcp_expt)+expt,index,:], label=control_titles[expt], color=control_colours[expt], linewidth=2)
        xlim([year_start, year_end])
        xlabel('Year', fontsize=14)
        ylabel('Gt/y', fontsize=14)
        grid(True)
        title(names[index] + '\nBasal Mass Loss (annually averaged)', fontsize=18)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.75, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1.1,0.5))
        fig.savefig(fig_names[index])

    


if __name__ == "__main__":

    timeseries_annual_compare()
    
        

    
