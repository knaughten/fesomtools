from numpy import *
from matplotlib.pyplot import *

def timeseries_annual_compare ():

    directory_head = '/short/y99/kaa561/FESOM/'
    rcp_expt = ['rcp45_M', 'rcp85_M'] #'rcp45_A', 'rcp85_M', 'rcp85_A']
    rcp_titles = ['RCP 4.5', 'RCP 8.5'] #['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A']
    rcp_colours = ['black', 'red'] #['red', 'magenta', 'blue', 'green']
    control_expt = 'lowres_spinup'
    control_title = 'CONTROL'
    control_colour = 'blue' #'black'
    year_start = 1992
    year_end = 2100
    control_skipyears = 28
    peryear = 365/5

    names = ['All Ice Shelves', 'Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    fig_names = ['total_massloss.png', 'larsen_d.png', 'larsen_c.png', 'wilkins_georgevi_stange.png', 'ronne_filchner.png', 'abbot.png', 'pig.png', 'thwaites.png', 'dotson.png', 'getz.png', 'nickerson.png', 'sulzberger.png', 'mertz.png', 'totten_moscowuni.png', 'shackleton.png', 'west.png', 'amery.png', 'princeharald.png', 'baudouin_borchgrevink.png', 'lazarev.png', 'nivl.png', 'fimbul_jelbart_ekstrom.png', 'brunt_riiserlarsen.png', 'ross.png']
    area = [1.41327580791e12, 9410595961.42, 48147893361.6, 46287951910.9, 429798928470.0, 27030080949.7, 3839594948.81, 2499220358.3, 3908582947.28, 29823059449.5, 4268520899.01, 11108310834.3, 3102730054.84, 4632897701.36, 26030138936.9, 11651566872.8, 64322690314.0, 2957848286.81, 40563562257.6, 6778604330.34, 5671169444.22, 52720412012.5, 72401508276.4, 475666675975.0]
    rho_ice = 916

    num_years = year_end - year_start + 1
    time = arange(year_start, year_end+1)

    dpt = empty([3, num_years]) #empty([5, num_years])
    for expt in range(2): #4):
        dpt_tmp = []
        f = open(directory_head + rcp_expt[expt] + '/dpt.log')
        f.readline()
        for line in f:
            dpt_tmp.append(float(line))
        f.close()
        for year in range(num_years):
            dpt[expt,year] = mean(array(dpt_tmp[peryear*year:peryear*(year+1)]))
    dpt_tmp = []
    f = open(directory_head + control_expt + '/dpt.log')
    f.readline()
    for line in f:
        dpt_tmp.append(float(line))
    f.close()
    dpt_tmp = dpt_tmp[control_skipyears*peryear:]
    for year in range(num_years):
        dpt[-1,year] = mean(array(dpt_tmp[peryear*year:peryear*(year+1)]))

    fig, ax = subplots(figsize=(10,6))
    for expt in range(2): #4):
        ax.plot(time, dpt[expt,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
    ax.plot(time, dpt[-1,:], label=control_title, color=control_colour, linewidth=2)
    title('Drake Passage Transport (annually averaged)', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel('Sv', fontsize=14)
    xlim([year_start, year_end])
    grid(True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('drakepstrans.png')

    area_min = empty([3, num_years]) #5, num_years])
    area_max = empty([3, num_years]) #5, num_years])
    volume_min = empty([3, num_years]) #5, num_years])
    volume_max = empty([3, num_years]) #5, num_years])
    for expt in range(2): #4):
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
    area_tmp = []
    volume_tmp = []
    f = open(directory_head + control_expt + '/seaice.log')
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
        area_min[-1,year] = amin(array(area_tmp[peryear*year:peryear*(year+1)]))
        area_max[-1,year] = amax(array(area_tmp[peryear*year:peryear*(year+1)]))
        volume_min[-1,year] = amin(array(volume_tmp[peryear*year:peryear*(year+1)]))
        volume_max[-1,year] = amax(array(volume_tmp[peryear*year:peryear*(year+1)]))

    fig, ax = subplots(figsize=(10,6))
    for expt in range(2): #4):
        ax.plot(time, area_min[expt,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
    ax.plot(time, area_min[-1,:], label=control_title, color=control_colour, linewidth=2)
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
    for expt in range(2): #4):
        ax.plot(time, area_max[expt,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
    ax.plot(time, area_max[-1,:], label=control_title, color=control_colour, linewidth=2)
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
    for expt in range(2): #4):
        ax.plot(time, volume_min[expt,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
    ax.plot(time, volume_min[-1,:], label=control_title, color=control_colour, linewidth=2)
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
    for expt in range(2): #4):
        ax.plot(time, volume_max[expt,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
    ax.plot(time, volume_max[-1,:], label=control_title, color=control_colour, linewidth=2)
    title('Annual maximum sea ice volume', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel(r'million km$^3$', fontsize=14)
    xlim([year_start, year_end])
    grid(True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('seaice_volume_max.png')

    massloss = empty([3, len(names), num_years]) #5, len(names), num_years])
    for expt in range(2): #4):
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
    f = open(directory_head + control_expt + '/massloss.log')
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
            massloss[-1,index,year] = mean(array(massloss_tmp[peryear*year:peryear*(year+1)]))
    f.close()

    for index in range(len(names)):
        factor = 1e12/(rho_ice*area[index])
        fig, ax1 = subplots(figsize=(10, 6)) #.5,6))
        for expt in range(2): #4):
            ax1.plot(time, massloss[expt,index,:], label=rcp_titles[expt], color=rcp_colours[expt], linewidth=2)
        ax1.plot(time, massloss[-1,index,:], label=control_title, color=control_colour, linewidth=2)
        ax1.set_xlim([year_start, year_end])
        ylimits = ax1.get_ylim()
        ax1.set_xlabel('Year', fontsize=14)
        ax1.set_ylabel('Gt/y', fontsize=14) #Basal Mass Loss (Gt/y)', fontsize=14)
        ax1.grid(True)
        #ax2 = ax1.twinx()
        #ax2.set_ylim([ylimits[0]*factor, ylimits[1]*factor])
        #ax2.set_ylabel('Area-Averaged Ice Shelf Melt Rate (m/y)', fontsize=14)
        title(names[index] + '\nBasal Mass Loss (annually averaged)', fontsize=18)
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width*0.75, box.height])
        #ax2.set_position([box.x0, box.y0, box.width*0.75, box.height])
        ax1.legend(loc='center left', bbox_to_anchor=(1.1,0.5))
        fig.savefig(fig_names[index])

    


if __name__ == "__main__":

    timeseries_annual_compare()
    
        

    
