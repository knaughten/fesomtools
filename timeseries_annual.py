from numpy import *
from matplotlib.pyplot import *

# Given the logfiles for 5-day timeseries created using timeseries_dpt.py,
# timeseries_seaice.py, and timeseries_massloss.py, plot annual averages
# (for Drake Passage transport and ice shelf melt rates/mass loss) or annual
# max/mins (for sea ice area and volume).
# Input:
# dpt_log = path to logfile from timeseries_dpt.py
# seaice_log = path to logfile from timeseries_seaice.py
# massloss_log = path to logfile from timeseries_massloss.py
# res_flag = integer flag indicating low resolution mesh (1) or high (2)
def timeseries_annual (dpt_log, seaice_log, massloss_log, res_flag):

    # Number of records per year (assumes 5-day averages)
    peryear = 365/5
    # Titles for each ice shelf
    names = ['All Ice Shelves', 'Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    # Beginning of figure names for each ice shelf
    fig_heads = ['total_massloss', 'larsen_d', 'larsen_c', 'wilkins_georgevi_stange', 'ronne_filchner', 'abbot', 'pig', 'thwaites', 'dotson', 'getz', 'nickerson', 'sulzberger', 'mertz', 'totten_moscowuni', 'shackleton', 'west', 'amery', 'princeharald', 'baudouin_borchgrevink', 'lazarev', 'nivl', 'fimbul_jelbart_ekstrom', 'brunt_riiserlarsen', 'ross']
    # Area of each ice shelf in m^2 (printed to screen during
    # timeseries_massloss.py, update if the mesh changes)
    if res_flag == 1:
        # Low resolution mesh
        area = [1.41327580791e12, 9410595961.42, 48147893361.6, 46287951910.9, 429798928470.0, 27030080949.7, 3839594948.81, 2499220358.3, 3908582947.28, 29823059449.5, 4268520899.01, 11108310834.3, 3102730054.84, 4632897701.36, 26030138936.9, 11651566872.8, 64322690314.0, 2957848286.81, 40563562257.6, 6778604330.34, 5671169444.22, 52720412012.5, 72401508276.4, 475666675975.0]
    elif res_flag == 2:
        # High resolution mesh
        area = [1.43919497627e12, 10456222697.1, 50141041705.5, 48253708618.4, 429898475473.0, 28942129634.4, 4388809435.81, 3172043475.79, 4290109356.52, 31663268041.5, 5985509656.36, 12669899186.6, 3911331361.75, 4974780745.66, 27070168363.3, 12236727597.8, 64795820721.7, 3070583821.88, 40914157792.4, 6896940796.67, 5942569502.6, 53559454524.6, 72896644960.9, 476899018047.0]
    # Observed mass loss (Rignot 2013) and uncertainty for each ice shelf, in
    # Gt/y
    obs_massloss = [1325, 1.4, 20.7, 135.4, 155.4, 51.8, 101.2, 97.5, 45.2, 144.9, 4.2, 18.2, 7.9, 90.6, 72.6, 27.2, 35.5, -2, 21.6, 6.3, 3.9, 26.8, 9.7, 47.7]
    obs_massloss_error = [235, 14, 67, 40, 45, 19, 8, 7, 4, 14, 2, 3, 3, 8, 15, 10, 23, 3, 18, 2, 2, 14, 16, 34]
    # Observed melt rate (Rignot 2013) and uncertainty for each ice shelf, in
    # m/y
    obs_ismr = [0.85, 0.1, 0.4, 3.1, 0.3, 1.7, 16.2, 17.7, 7.8, 4.3, 0.6, 1.5, 1.4, 7.7, 2.8, 1.7, 0.6, -0.4, 0.4, 0.7, 0.5, 0.5, 0.1, 0.1]
    obs_ismr_error = [0.1, 0.6, 1, 0.8, 0.1, 0.6, 1, 1, 0.6, 0.4, 0.3, 0.3, 0.6, 0.7, 0.6, 0.7, 0.4, 0.6, 0.4, 0.2, 0.2, 0.2, 0.2, 0.1]
    # Density of ice in kg/m^3
    rho_ice = 916

    # Drake Passage transport
    dpt = []
    # Read log file
    f = open(dpt_log, 'r')
    # Skip header
    f.readline()
    for line in f:
        dpt.append(float(line))
    f.close()
    # Calculate how many years of output there are
    num_output = len(dpt)
    num_years = num_output/peryear
    print str(num_years) + ' years of output'
    # Calculate annual averages
    dpt_avg = []
    for year in range(num_years):
        dpt_avg.append(mean(array(dpt[peryear*year:peryear*(year+1)])))
    # Make time array
    time = range(num_years)
    # Plot
    clf()
    plot(time, dpt_avg)
    xlabel('Years')
    ylabel('Drake Passage Transport (Sv)')
    grid(True)
    savefig('drakepsgtrans_avg.png')

    # Sea ice area and volume
    seaice_area = []
    seaice_volume = []
    # Read log file
    f = open(seaice_log, 'r')
    # Skip header for sea ice area
    f.readline()
    for line in f:
        try:
            seaice_area.append(float(line))
        except(ValueError):
            # Header for sea ice volume
            break
    for line in f:
        seaice_volume.append(float(line))
    f.close()
    # Calculate annual max/mins
    area_min = []
    area_max = []
    volume_min = []
    volume_max = []
    for year in range(num_years):
        area_min.append(amin(array(seaice_area[peryear*year:peryear*(year+1)])))
        area_max.append(amax(array(seaice_area[peryear*year:peryear*(year+1)])))
        volume_min.append(amin(array(seaice_volume[peryear*year:peryear*(year+1)])))
        volume_max.append(amax(array(seaice_volume[peryear*year:peryear*(year+1)])))
    # Plot
    clf()
    plot(time, area_min)
    plot(time, area_max)
    xlabel('Years')
    ylabel(r'Annual min and max sea ice area (million km$^2$)')
    grid(True)
    savefig('seaice_area_minmax.png')
    clf()
    plot(time, volume_min)
    plot(time, volume_max)
    xlabel('Years')
    ylabel(r'Annual min and max sea ice volume (million km$^3$)')
    grid(True)
    savefig('seaice_volume_minmax.png')

    # Mass loss
    massloss = empty([len(names), num_output])
    # Read logfile
    f = open(massloss_log, 'r')
    # Skip header
    f.readline()
    index = 0
    # Loop over ice shelves
    while index < len(names):
        t = 0
        for line in f:
            try:
                massloss[index,t] = float(line)
                t += 1
            except(ValueError):
                # Reached the header for the next ice shelf
                break
        index += 1
    f.close()
    # Calculate annual averages
    massloss_avg = empty([len(names), num_years])
    for index in range(len(names)):
        for year in range(num_years):
            massloss_avg[index,year] = mean(massloss[index,peryear*year:peryear*(year+1)])
    # Plot each ice shelf
    for index in range(len(names)):
        # Calculate conversion factor from mass loss to area-averaged melt rate
        factor = 1e12/(rho_ice*area[index])
        # Calculate the bounds on observed mass loss and melt rates
        massloss_low = obs_massloss[index] - obs_massloss_error[index]
        massloss_high = obs_massloss[index] + obs_massloss_error[index]
        ismr_low = obs_ismr[index] - obs_ismr_error[index]
        ismr_high = obs_ismr[index] + obs_ismr_error[index]
        # Set up plot: mass loss and melt rate are directly proportional
        # (with a different constant of proportionality for each ice shelf
        # depending on its area) so plot one line with two y-axes
        fig, ax1 = subplots()
        ax1.plot(time, massloss_avg[index,:], color='black')
        # In blue, add dashed lines for observed mass loss
        ax1.axhline(massloss_low, color='b', linestyle='dashed')
        ax1.axhline(massloss_high, color='b', linestyle='dashed')
        # Make sure y-limits won't cut off observed melt rate
        ymin = amin([ismr_low/factor, massloss_low, amin(massloss_avg[index,:])])
        ymax = amax([ismr_high/factor, massloss_high, amax(massloss_avg[index,:])])
        # Adjust y-limits to line up with ticks
        ticks = ax1.get_yticks()
        min_tick = ticks[0]
        max_tick = ticks[-1]
        dtick = ticks[1]-ticks[0]
        while min_tick >= ymin:
            min_tick -= dtick
        while max_tick <= ymax:
            max_tick += dtick
        ax1.set_ylim([min_tick, max_tick])
        # Title and ticks in blue for this side of the plot
        ax1.set_ylabel('Basal Mass Loss (Gt/y)', color='b')
        for t1 in ax1.get_yticklabels():
            t1.set_color('b')
        ax1.set_xlabel('Years')
        ax1.grid(True)
        # Twin axis for melt rates
        ax2 = ax1.twinx()
        # Make sure the scales line up
        limits = ax1.get_ylim()
        ax2.set_ylim([limits[0]*factor, limits[1]*factor])
        # In red, add dashed lines for observed ice shelf melt rates
        ax2.axhline(ismr_low, color='r', linestyle='dashed')
        ax2.axhline(ismr_high, color='r', linestyle='dashed')
        # Title and ticks in red for this side of the plot
        ax2.set_ylabel('Area-Averaged Ice Shelf Melt Rate (m/y)', color='r')
        for t2 in ax2.get_yticklabels():
            t2.set_color('r')
        # Name of the ice shelf for the main title
        title(names[index])
        fig.savefig(fig_heads[index] + '_avg.png')


# Command-line interface
if __name__ == '__main__':

    dpt_log = raw_input("Path to logfile for timeseries_dpt.py: ")
    seaice_log = raw_input("Path to logfile for timeseries_seaice.py: ")
    massloss_log = raw_input("Path to logfile for timeseries_massloss.py: ")
    res_flag = int(raw_input("Low resolution (1) or high (2)? "))
    timeseries_annual(dpt_log, seaice_log, massloss_log, res_flag)

    
        
        
