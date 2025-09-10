from numpy import *

def bugs_calc_aice_min (extent_log_old, extent_log_new):

    start_year = 1992
    end_year = 2016
    num_years = end_year - start_year + 1

    extent_5day_old = []
    f = open(extent_log_old, 'r')
    f.readline()
    for line in f:
        extent_5day_old.append(float(line))
    f.close()
    feb_extent_old = zeros(num_years)
    for year in range(start_year, end_year+1):
        # First timestep of year in 5-day logfile
        t0 = (year-start_year)*73
        # Feburary: 4/5 of index 7, indices 8-11, and 4/5 of index 12
        feb_extent_old[year-start_year] = (extent_5day_old[t0+6]*4 + sum(extent_5day_old[t0+7:t0+11]*5) + extent_5day_old[t0+11]*4)/28.0
    mean_old = mean(feb_extent_old)
    print('Mean February sea ice extent in old simulation: ' + str(mean_old))

    extent_5day_new = []
    f = open(extent_log_new, 'r')
    f.readline()
    for line in f:
        extent_5day_new.append(float(line))
    f.close()
    feb_extent_new = zeros(num_years)
    for year in range(start_year, end_year+1):
        t0 = (year-start_year)*73
        feb_extent_new[year-start_year] = (extent_5day_new[t0+6]*4 + sum(extent_5day_new[t0+7:t0+11]*5) + extent_5day_new[t0+11]*4)/28.0
    mean_new = mean(feb_extent_new)
    print('Mean February sea ice extent in new simulation: ' + str(mean_new))
    print('Increase of ' + str((mean_new-mean_old)/mean_old*100) + '%')


# Command-line interface
if __name__ == "__main__":

    extent_log_old = input("Path to sea ice extent timeseries from old simulation: ")
    extent_log_new = input("Path to sea ice extent timeseries from new simulation: ")
    bugs_calc_aice_min(extent_log_old, extent_log_new)
    

    
