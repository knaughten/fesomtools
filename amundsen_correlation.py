from numpy import *
from scipy.stats import linregress

def amundsen_correlation ():

    # Paths to RCP experiment directories
    directory_head = '/short/y99/kaa561/FESOM/'
    rcp_expt = ['rcp45_M_highres/', 'rcp45_A_highres/', 'rcp85_M_highres/', 'rcp85_A_highres/']
    num_rcps = len(rcp_expt)
    # Titles 
    rcp_titles = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A']
    # Path to control experiment directory
    control_expt = 'highres_spinup/'
    # Title 
    control_title = 'CONTROL'
    # Logfile names
    amundsen_log = 'amundsen.log'
    massloss_log = 'massloss.log'
    # Index of PIG in timeseries_massloss.py
    pig_index = 6
    # Years to consider
    year_start = 1992
    year_end = 2100
    # Spinup years to discard (first 2 repetitions of 1992-2005)
    control_skipyears = 28
    # Output steps per year
    peryear = 365/5

    num_years = year_end - year_start + 1

    # Amundsen Sea ice-to-ocean freshwater flux
    ice2ocn = empty([num_rcps+1, num_years])
    # Loop over RCP experiments
    for expt in range(num_rcps):
        # Read logfile
        ice2ocn_tmp = []
        f = open(directory_head + rcp_expt[expt] + amundsen_log, 'r')
        f.readline()
        for line in f:
            ice2ocn_tmp.append(float(line))
        f.close()
        # Calculate annual averages
        for year in range(num_years):
            ice2ocn[expt,year] = mean(array(ice2ocn_tmp[peryear*year:peryear*(year+1)]))
    # Control experiment
    ice2ocn_tmp = []
    f = open(directory_head + control_expt + amundsen_log, 'r')
    f.readline()
    for line in f:
        ice2ocn_tmp.append(float(line))
    f.close()
    # Throw away first 2 repetitions
    ice2ocn_tmp = ice2ocn_tmp[control_skipyears*peryear:]
    for year in range(num_years):
        ice2ocn[-1,year] = mean(array(ice2ocn_tmp[peryear*year:peryear*(year+1)]))

    # PIG basal mass loss
    pig_massloss = empty([num_rcps+1, num_years])
    for expt in range(num_rcps):
        f = open(directory_head + rcp_expt[expt] + massloss_log, 'r')
        f.readline()
        # Read the ice shelves before PIG
        for index in range(pig_index):
            for line in f:
                try:
                    tmp = float(line)
                except(ValueError):
                    break
        # Now read PIG
        pig_massloss_tmp = []
        for line in f:
            try:
                pig_massloss_tmp.append(float(line))
            except(ValueError):
                break
        f.close()
        for year in range(num_years):
            pig_massloss[expt,year] = mean(array(pig_massloss_tmp[peryear*year:peryear*(year+1)]))
    f = open(directory_head + control_expt + massloss_log, 'r')
    f.readline()
    for index in range(pig_index):
        for line in f:
            try:
                tmp = float(line)
            except(ValueError):
                break
    pig_massloss_tmp = []
    for line in f:
        try:
            pig_massloss_tmp.append(float(line))
        except(ValueError):
            break
    f.close()
    pig_massloss_tmp = pig_massloss_tmp[control_skipyears*peryear:]
    for year in range(num_years):
        pig_massloss[-1,year] = mean(array(pig_massloss_tmp[peryear*year:peryear*(year+1)]))

    for expt in range(num_rcps):
        slope, intercept, r_value, p_value, std_err = linregress(ice2ocn[expt,:], pig_massloss[expt,:])
        print rcp_titles[expt] + ': r^2 = ' + str(r_value**2)
    slope, intercept, r_value, p_value, std_err = linregress(ice2ocn[-1,:], pig_massloss[-1,:])
    print control_title + ': r^2 = ' + str(r_value**2)


# Command-line interface
if __name__ == "__main__":

    amundsen_correlation()
        
    
