from numpy import *

def watermass_combine_mcdw (log_file):

    # Sectors
    sector_names = ['Filchner-Ronne Ice Shelf Cavity', 'Eastern Weddell Region Cavities', 'Amery Ice Shelf Cavity', 'Australian Sector Cavities', 'Ross Sea Cavities', 'Amundsen Sea Cavities', 'Bellingshausen Sea Cavities', 'Larsen Ice Shelf Cavities', 'All Ice Shelf Cavities']
    num_sectors = len(sector_names)
    # Water masses
    wm_names = ['ISW', 'AASW', 'CDW', 'MCDW', 'WW', 'HSSW']
    num_watermasses = len(wm_names)
    # Indices of CDW and MCDW
    cdw_key = wm_names.index('CDW')
    mcdw_key = wm_names.index('MCDW')

    # Read the first variable just to figure out how many timesteps there are
    num_years = 0
    f = open(log_file, 'r')
    f.readline()
    for line in f:
        try:
            tmp = float(line)
            num_years += 1
        except(ValueError):
            break
    f.close()
    # Set up array to hold timeseries
    percent_watermass = empty([num_watermasses, num_sectors, num_years])
    # Read timeseries
    f = open(log_file, 'r')
    f.readline()
    wm_key = 0
    while wm_key < num_watermasses:
        sector = 0
        while sector < num_sectors:
            year = 0
            for line in f:
                try:
                    percent_watermass[wm_key, sector, year] = float(line)
                    year += 1
                except(ValueError):
                    break
            sector += 1
        wm_key += 1
    f.close()

    # Add CDW to MCDW
    percent_watermass_new = empty([num_watermasses-1, num_sectors, num_years])
    for wm_key in range(cdw_key):
        percent_watermass_new[wm_key,:,:] = percent_watermass[wm_key,:,:]
    percent_watermass_new[cdw_key,:,:] = percent_watermass[cdw_key,:,:] + percent_watermass[mcdw_key,:,:]
    for wm_key in range(mcdw_key+1, num_watermasses):
        percent_watermass_new[wm_key-1,:,:] = percent_watermass[wm_key,:,:]

    # Update list of water mass names
    wm_names.remove('CDW')
    num_watermasses = len(wm_names)

    # Write to file
    f = open(log_file, 'w')
    for wm_key in range(num_watermasses):
        for sector in range(num_sectors):
            f.write(wm_names[wm_key] + 'in ' + sector_names[sector] + '(%)\n')
            for t in range(num_years):
                f.write(str(percent_watermass_new[wm_key, sector, t]) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    log_file = input("Path to logfile from timeseries_watermass_sectors.py: ")
    watermass_combine_mcdw(log_file)
                
        
    
    

    
