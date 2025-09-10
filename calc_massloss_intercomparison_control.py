from numpy import *

def calc_massloss_intercomparison_control (intercomparison_logfile, control_logfile):

    # Year simulations start
    year_start = 1992
    # Years to average over
    calc_start = 1992
    calc_end = 2005
    # Number of output steps per year in FESOM
    peryear = 365//5
    # Name of each ice shelf
    names = ['Total Mass Loss', 'Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange & Bach Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    num_shelves = len(names)

    # Read timeseries
    # Intercomparison (high-res)
    tmp = []
    f = open(intercomparison_logfile, 'r')
    # Skip the first line (header)
    f.readline()
    # Read total mass loss
    num_time = 0
    for line in f:
        try:
            tmp.append(float(line))
            num_time += 1
        except(ValueError):
            # Reached the header for the next variable
            break
    # Set up array for mass loss values at each ice shelf
    massloss_ts_intercomparison = empty([num_shelves, num_time])
    # Save the total values in the first index
    massloss_ts_intercomparison[0,:] = array(tmp)
    # Loop over ice shelves
    index = 1
    while index < num_shelves:
        t = 0
        for line in f:
            try:
                massloss_ts_intercomparison[index,t] = float(line)
                t += 1
            except(ValueError):
                # Reached the header for the next ice shelf
                break
        index += 1
    f.close()
    # Average between given years
    massloss_intercomparison = mean(massloss_ts_intercomparison[:,peryear*(calc_start-year_start):peryear*(calc_end+1-year_start)], axis=1)
    # Control simulation from RCP paper
    tmp = []
    f = open(control_logfile, 'r')
    f.readline()
    num_time = 0
    for line in f:
        try:
            tmp.append(float(line))
            num_time += 1
        except(ValueError):
            break
    massloss_ts_control = empty([num_shelves, num_time])
    massloss_ts_control[0,:] = array(tmp)
    index = 1
    while index < num_shelves:
        t = 0
        for line in f:
            try:
                massloss_ts_control[index,t] = float(line)
                t += 1
            except(ValueError):
                break
        index += 1
    f.close()
    massloss_control = mean(massloss_ts_control[:,peryear*(calc_start-year_start):peryear*(calc_end+1-year_start)], axis=1)

    # Print results
    for index in range(num_shelves):
        print(names[index])
        print('Intercomparison: ' + str(massloss_intercomparison[index]))
        print('Control: ' + str(massloss_control[index]))
        print(str((massloss_control[index]-massloss_intercomparison[index])/massloss_intercomparison[index]*100) + '% change')
        print('\n')


# Command-line interface
if __name__ == "__main__":

    intercomparison_logfile = input("Path to massloss logfile for high-res intercomparison experiment: ")
    control_logfile = input("Path to massloss logfile for RCP control experiment: ")
    calc_massloss_intercomparison_control(intercomparison_logfile, control_logfile)
