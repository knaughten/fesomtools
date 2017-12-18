from numpy import *

def seaice_formation_concatenate_control ():

    # File paths
    directory_head = '/short/y99/kaa561/FESOM/highres_spinup/'
    log_name = 'seaice_formation.log'
    # Years to consider
    start_year = 1992
    end_year = 2100
    # Spinup forcing repetitions to consider
    start_rep = 3
    end_rep = 10

    num_years = end_year-start_year+1
    # Set up array to hold timeseries
    formation = empty(num_years)
    # Read all but the last repetition
    year = 0
    for rep in range(start_rep, end_rep):
        f = open(directory_head + 'rep' + str(rep) + '/' + log_name, 'r')
        f.readline()
        t = year
        for line in f:
            formation[t] = float(line)
            t += 1
        year = t
        f.close()
    # Read what we need from the last repetition
    f = open(directory_head + 'rep' + str(end_rep) + '/' + log_name, 'r')
    f.readline()
    t = year
    for line in f:
        if t < num_years:
            formation[t] = float(line)
            t += 1
    f.close()

    # Write output
    f = open(directory_head + log_name, 'w')
    f.write('Net sea ice formation on continental shelf (thousand km^3/y):\n')
    for t in range(num_years):
        f.write(str(formation[t]) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    seaice_formation_concatenate_control()
