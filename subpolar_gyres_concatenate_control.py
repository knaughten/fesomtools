from numpy import *

def subpolar_gyres_concatenate_control ():

    # File paths
    directory_head = '/short/y99/kaa561/FESOM/highres_spinup/'
    log_name = 'subpolar_gyres.log'
    # Years to consider
    start_year = 1992
    end_year = 2100
    # Spinup forcing repetitions to consider
    start_rep = 3
    end_rep = 10

    num_years = end_year-start_year+1
    # Set up arrays to hold timeseries
    ws_trans = empty(num_years)
    rs_trans = empty(num_years)
    # Read all but the last repetition
    year = 0
    for rep in range(start_rep, end_rep):
        f = open(directory_head + 'rep' + str(rep) + '/' + log_name, 'r')
        f.readline()
        t = year
        for line in f:
            try:
                ws_trans[t] = float(line)
                t += 1
            except(ValueError):
                break
        t = year
        for line in f:
            try:
                rs_trans[t] = float(line)
                t += 1
            except(ValueError):
                break
        year = t
    # Read what we need from the last repetition
    f = open(directory_head + 'rep' + str(end_rep) + '/' + log_name, 'r')
    f.readline()
    t = year
    for line in f:
        try:
            tmp = float(line)
        except(ValueError):
            break
        if t < num_years:
            ws_trans[t] = tmp
            t += 1
    t = year
    for line in f:
        try:
            tmp = float(line)
        except(ValueError):
            break
        if t < num_years:
            rs_trans[t] = tmp
            t += 1
    f.close()

    # Write output
    f = open(directory_head + log_name, 'w')
    f.write('Weddell Sea Gyre transport (Sv)\n')
    for t in range(num_years):
        f.write(str(ws_trans[t]) + '\n')
    f.write('Ross Sea Gyre transport (Sv)\n')
    for t in range(num_years):
        f.write(str(rs_trans[t]) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    subpolar_gyres_concatenate_control()
