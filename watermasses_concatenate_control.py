from numpy import *

def watermasses_concatenate_control ():

    # File paths
    directory_head = '/short/y99/kaa561/FESOM/highres_spinup/'
    log_name = 'water_masses.log'
    # Years to consider
    start_year = 1992
    end_year = 2100
    # Spinup forcing repetitions to consider
    start_rep = 3
    end_rep = 10
    # Sectors
    sector_names = ['Filchner-Ronne Ice Shelf', 'Eastern Weddell Region', 'Amery Ice Shelf', 'Australian Sector', 'Ross Sea', 'Amundsen Sea', 'Bellingshausen Sea', 'Larsen Ice Shelf', 'All Ice Shelves']
    num_sectors = len(sector_names)
    # Water masses
    wm_names = ['ISW', 'AASW', 'CDW', 'MCDW', 'WW', 'HSSW']
    num_watermasses = len(wm_names)

    num_years = end_year-start_year+1

    # Set up array to hold timeseries
    percent_watermass = empty([num_watermasses, num_sectors, num_years])
    # Read all but the last repetition
    year = 0
    for rep in range(start_rep, end_rep):
        f = open(directory_head + 'rep' + str(rep) + '/' + log_name, 'r')
        f.readline()
        wm_key = 0
        while wm_key < num_watermasses:
            sector = 0
            while sector < num_sectors:
                t = year
                for line in f:
                    try:
                        percent_watermass[wm_key, sector, t] = float(line)
                        t += 1
                    except(ValueError):
                        break
                sector += 1
            wm_key += 1
        year = t
    # Read what we need from the last repetition
    f = open(directory_head + 'rep' + str(end_rep) + '/' + log_name, 'r')
    f.readline()
    wm_key = 0
    while wm_key < num_watermasses:
        sector = 0
        while sector < num_sectors:
            t = year
            for line in f:
                try:
                    tmp = float(line)
                except(ValueError):
                    break                
                if t < num_years:
                    percent_watermass[wm_key, sector, t] = tmp
                    t += 1
            sector += 1
        wm_key += 1
    f.close()

    # Write output
    f = open(directory_head + log_name, 'w')
    for wm_key in range(num_watermasses):
        for sector in range(num_sectors):
            f.write(wm_names[wm_key] + 'in ' + sector_names[sector] + '(%)\n')
            for t in range(num_years):
                f.write(str(percent_watermass[wm_key, sector, t]) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    watermasses_concatenate_control()
