from numpy import *

def dpt_change ():

    # Paths to RCP experiment directories
    directory_head = '/short/y99/kaa561/FESOM/'
    rcp_expt = ['rcp45_M/', 'rcp45_A/', 'rcp85_M/', 'rcp85_A/']
    num_rcps = len(rcp_expt)
    # Titles 
    rcp_titles = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A']
    # Path to control experiment directory
    control_expt = 'highres_spinup/'
    # Years to consider
    year_start = 1992
    present_end = 2005
    rcp_end = 2100
    num_years_present = present_end - year_start + 1
    num_years_rcp = rcp_end - year_start + 1
    # Spinup years to discard (first 2 repetitions of 1992-2005)
    control_skipyears = 28
    # Output steps per year
    peryear = 365/5

    f = open(directory_head + control_expt + 'dpt.log', 'r')
    f.readline()
    dpt_tmp = []
    for line in f:
        dpt_tmp.append(float(line))
    f.close()
    dpt_tmp = dpt_tmp[control_skipyears*peryear:]
    dpt_baseline = zeros(num_years_present)
    for year in range(num_years_present):
        dpt_baseline[year] = mean(array(dpt_tmp[peryear*year:peryear*(year+1)]))
    dpt_beg = mean(dpt_baseline[-10:])
    print '1996-2005: ' + str(dpt_beg)
    dpt_control = zeros(num_years_rcp)
    for year in range(num_years_rcp):
        dpt_control[year] = mean(array(dpt_tmp[peryear*year:peryear*(year+1)]))
    dpt_drift = mean(dpt_control[-10:])
    print '2091-2100, CONTROL : ' + str(dpt_drift) + ', change of ' + str((dpt_drift-dpt_beg)/dpt_beg*100) + '%'

    dpt = empty([num_rcps, num_years_rcp])
    for expt in range(num_rcps):
        f = open(directory_head + rcp_expt[expt] + 'dpt.log', 'r')
        f.readline()
        dpt_tmp = []
        for line in f:
            dpt_tmp.append(float(line))
        f.close()
        for year in range(num_years_rcp):
            dpt[expt, year] = mean(array(dpt_tmp[peryear*year:peryear*(year+1)]))
        dpt_end = mean(dpt[expt,-10:])
        print '2091-2100, ' + str(rcp_titles[expt]) + ': ' + str(dpt_end) + ', change of ' + str((dpt_end-dpt_beg)/dpt_beg*100) + '%'


# Command-line interface
if __name__ == "__main__":

    dpt_change()
