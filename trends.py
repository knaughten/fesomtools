from numpy import *
from scipy.stats import linregress
from matplotlib.pyplot import *

def trends ():

    # Paths to RCP experiment directories
    directory_head = '/short/y99/kaa561/FESOM/'
    rcp_expt = ['rcp45_M/', 'rcp45_A/', 'rcp85_M/', 'rcp85_A/']
    num_rcps = len(rcp_expt)
    # Titles 
    rcp_titles = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A']
    # Path to control experiment directory
    control_expt = 'highres_spinup/'
    # Title
    control_title = 'CONTROL'
    # Years to consider
    calc_start = 2006
    calc_end = 2100
    # Logfile names
    massloss_log = 'massloss_sectors.log'
    dpt_log = 'dpt.log'
    subpolar_log = 'subpolar_gyres.log'
    # Present-day years to discard (1992-2005)
    beg_skipyears = 14
    # Spinup years to discard (first 2 repetitions of 1992-2005)
    control_skipyears = 28
    # Output steps per year
    peryear = 365/5
    # Names of sectors
    sector_names = ['Filchner-Ronne Ice Shelf', 'Eastern Weddell Region', 'Amery Ice Shelf', 'Australian Sector', 'Ross Sea', 'Amundsen Sea', 'Bellingshausen Sea', 'Larsen Ice Shelves', 'Total Antarctica']
    num_sectors = len(sector_names)

    time = arange(calc_start, calc_end+1)
    num_years = calc_end - calc_start + 1

    # Read all the timeseries for the correct years

    # Mass loss timeseries for each sector
    massloss = empty([num_rcps+1, num_sectors, num_years])
    # 1992-2005 average for each sector
    massloss_baseline = empty(num_sectors)
    # Loop over RCP experiments
    for expt in range(num_rcps):
        # Read logfile
        f = open(directory_head + rcp_expt[expt] + massloss_log, 'r')
        f.readline()
        # Loop over sectors
        for sector in range(num_sectors):
            massloss_tmp = []
            for line in f:
                try:
                    massloss_tmp.append(float(line))
                except(ValueError):
                    break
            # Throw away 1992-2005
            massloss_tmp = massloss_tmp[beg_skipyears*peryear:]
            # Calculate annual averages
            for year in range(num_years):
                massloss[expt, sector, year] = mean(array(massloss_tmp[peryear*year:peryear*(year+1)]))
        f.close()
    # Control experiment
    f = open(directory_head + control_expt + massloss_log, 'r')
    f.readline()
    for sector in range(num_sectors):
        massloss_tmp = []
        for line in f:
            try:
                massloss_tmp.append(float(line))
            except(ValueError):
                break
        # Throw away first two repetitions
        massloss_tmp = massloss_tmp[control_skipyears*peryear:]
        # Calculate 1992-2005 average
        massloss_baseline[sector] = mean(array(massloss_tmp[:beg_skipyears*peryear]))
        # Throw away 1992-2005
        massloss_tmp = massloss_tmp[beg_skipyears*peryear:]
        # Calculate annual averages
        for year in range(num_years):
            massloss[-1, sector, year] = mean(array(massloss_tmp[peryear*year:peryear*(year+1)]))
    f.close()
    # Now scale as a percentage of 1992-2005 average for each sector
    for sector in range(num_sectors):
        massloss[:,sector,:] = massloss[:,sector,:]/massloss_baseline[sector]*100

    # Drake Passage transport
    dpt = empty([num_rcps+1, num_years])
    for expt in range(num_rcps):
        f = open(directory_head + rcp_expt[expt] + dpt_log, 'r')
        f.readline()
        dpt_tmp = []
        for line in f:
            dpt_tmp.append(float(line))
        f.close()
        dpt_tmp = dpt_tmp[beg_skipyears*peryear:]
        for year in range(num_years):
            dpt[expt, year] = mean(array(dpt_tmp[peryear*year:peryear*(year+1)]))
    f = open(directory_head + control_expt + dpt_log, 'r')
    f.readline()
    dpt_tmp = []
    for line in f:
        dpt_tmp.append(float(line))
    f.close()
    dpt_tmp = dpt_tmp[control_skipyears*peryear:]
    dpt_baseline = mean(array(dpt_tmp[:beg_skipyears*peryear]))
    dpt_tmp = dpt_tmp[beg_skipyears*peryear:]
    for year in range(num_years):
        dpt[-1, year] = mean(array(dpt_tmp[peryear*year:peryear*(year+1)]))
    dpt = dpt/dpt_baseline*100

    # Subpolar gyre transport
    # These are already annually averaged: 2006-2100 for the RCPs, 1992-2100
    # for the control experiment
    ws_trans = empty([num_rcps+1, num_years])
    rs_trans = empty([num_rcps+1, num_years])
    for expt in range(num_rcps):
        f = open(directory_head + rcp_expt[expt] + subpolar_log, 'r')
        f.readline()
        year = 0
        for line in f:
            try:
                ws_trans[expt, year] = float(line)
                year += 1
            except(ValueError):
                break
        year = 0
        for line in f:
            rs_trans[expt, year] = float(line)
            year += 1
        f.close()
    f = open(directory_head + control_expt + subpolar_log, 'r')
    f.readline()
    ws_trans_tmp = []
    for line in f:
        try:
            ws_trans_tmp.append(float(line))
        except(ValueError):
            break
    rs_trans_tmp = []
    for line in f:
        rs_trans_tmp.append(float(line))
    # Calculate 1992-2005 averages
    ws_trans_baseline = mean(array(ws_trans_tmp[:beg_skipyears]))
    rs_trans_baseline = mean(array(rs_trans_tmp[:beg_skipyears]))
    # Throw away 1992-2005
    ws_trans[-1,:] = ws_trans_tmp[beg_skipyears:]
    rs_trans[-1,:] = rs_trans_tmp[beg_skipyears:]
    # Scale against 1992-2005 average
    ws_trans = ws_trans/ws_trans_baseline*100
    rs_trans = rs_trans/rs_trans_baseline*100

    '''# Calculate and print trends, including p-values
    print 'Trends in basal mass loss:\n'
    for sector in range(num_sectors):
        print '  ' + sector_names[sector] + ':'
        for expt in range(num_rcps):
            slope, intercept, r_value, p_value, std_err = linregress(time, massloss[expt,sector,:])
            print '    ' + rcp_titles[expt] + ': ' + str(slope*10) + ' %/decade, p=' + str(p_value)
        slope, intercept, r_value, p_value, std_err = linregress(time, massloss[-1,sector,:])
        print '    ' + control_title + ': ' + str(slope*10) + ' %/decade, p=' + str(p_value)
        print '\n'
    print 'Trends in Drake Passage transport:'
    for expt in range(num_rcps):
        slope, intercept, r_value, p_value, std_err = linregress(time, dpt[expt,:])
        print '    ' + rcp_titles[expt] + ': ' + str(slope*10) + ' %/decade, p=' + str(p_value)
    slope, intercept, r_value, p_value, std_err = linregress(time, dpt[-1,:])
    print '    ' + control_title + ': ' + str(slope*10) + ' %/decade, p=' + str(p_value)
    print '\n'
    print 'Trends in Weddell Sea Gyre transport:'
    for expt in range(num_rcps):
        slope, intercept, r_value, p_value, std_err = linregress(time, ws_trans[expt,:])
        print '    ' + rcp_titles[expt] + ': ' + str(slope*10) + ' %/decade, p=' + str(p_value)
    slope, intercept, r_value, p_value, std_err = linregress(time, ws_trans[-1,:])
    print '    ' + control_title + ': ' + str(slope*10) + ' %/decade, p=' + str(p_value)
    print '\n'
    print 'Trends in Ross Sea Gyre transport:'
    for expt in range(num_rcps):
        slope, intercept, r_value, p_value, std_err = linregress(time, rs_trans[expt,:])
        print '    ' + rcp_titles[expt] + ': ' + str(slope*10) + ' %/decade, p=' + str(p_value)
    slope, intercept, r_value, p_value, std_err = linregress(time, rs_trans[-1,:])
    print '    ' + control_title + ': ' + str(slope*10) + ' %/decade, p=' + str(p_value)'''

    print 'r-squared values for mass loss: '
    r2_min = 1
    r2_max = 0
    r2_sum = 0
    r2_num = 0
    for sector in range(num_sectors-1):
        print '  ' + sector_names[sector] + ':'
        for expt in range(num_rcps):
            slope, intercept, r_value, p_value, std_err = linregress(time, massloss[expt,sector,:])
            r2_linear = r_value**2
            if p_value < 0.05:
                r2_min = min(r2_min, r2_linear)
                r2_max = max(r2_max, r2_linear)
                r2_sum += r2_linear
                r2_num += 1
            '''# Now quadratic fit (a bit more work)
            y = massloss[expt,sector,:]
            coeffs = polyfit(time, y, 2)
            p = poly1d(coeffs)
            yhat = p(time)
            ybar = mean(y)
            ssreg = sum((yhat-ybar)**2)
            sstot = sum((y-ybar)**2)
            r2_quadratic = ssreg/sstot'''
    print 'r-squared values from ' + str(r2_min) + ' to ' + str(r2_max) + ', mean of ' + str(r2_sum/r2_num)
    


# Command-line interface
if __name__ == "__main__":

    trends()
    

    
            
    
            
