from numpy import *

def timeseries_concat (spinup_path, rcp_path):

    names = ['All Ice Shelves', 'Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']

    year_start = 1992
    year_end = 2005
    skipyears = 28
    peryear = 365/5
    numyears = year_end - year_start + 1

    dpt = []
    f = open(spinup_path + 'dpt.log', 'r')
    f.readline()
    for line in f:
        dpt.append(float(line))
    f.close()
    dpt = dpt[skipyears*peryear:(skipyears+numyears)*peryear]
    num_time_spinup = len(dpt)
    f = open(rcp_path + 'dpt.log', 'r')
    f.readline()
    for line in f:
        dpt.append(float(line))
    f.close()
    num_time = len(dpt)
    f = open(rcp_path + 'dpt.log', 'w')
    f.write('Drake Passage Transport (Sv):\n')
    for elm in dpt:
        f.write(str(elm) + '\n')
    f.close()

    seaice_area = []
    seaice_volume = []
    f = open(spinup_path + 'seaice.log', 'r')
    f.readline()
    for line in f:
        try:
            seaice_area.append(float(line))
        except(ValueError):
            break
    for line in f:
        seaice_volume.append(float(line))
    f.close()
    seaice_area = seaice_area[skipyears*peryear:(skipyears+numyears)*peryear]
    seaice_volume = seaice_volume[skipyears*peryear:(skipyears+numyears)*peryear]
    f = open(rcp_path + 'seaice.log', 'r')
    f.readline()
    for line in f:
        try:
            seaice_area.append(float(line))
        except(ValueError):
            break
    for line in f:
        seaice_volume.append(float(line))
    f.close()
    f = open(rcp_path + 'seaice.log', 'w')
    f.write('Total Sea Ice Area (million km^2):\n')
    for elm in seaice_area:
        f.write(str(elm) + '\n')
    f.write('Total Sea Ice Volume (million km^3):\n')
    for elm in seaice_volume:
        f.write(str(elm) + '\n')
    f.close()

    massloss = empty([len(names), num_time])
    f = open(spinup_path + 'massloss.log', 'r')
    f.readline()
    for index in range(len(names)):
        massloss_tmp = []
        for line in f:
            try:
                massloss_tmp.append(float(line))
            except(ValueError):
                break
        massloss_tmp = massloss_tmp[skipyears*peryear:(skipyears+numyears)*peryear]
        massloss[index, 0:num_time_spinup] = massloss_tmp
    f.close()
    f = open(rcp_path + 'massloss.log', 'r')
    f.readline()
    for index in range(len(names)):
        massloss_tmp = []
        for line in f:
            try:
                massloss_tmp.append(float(line))
            except(ValueError):
                break
        massloss[index, num_time_spinup:] = massloss_tmp
    f.close()
    f = open(rcp_path + 'massloss.log', 'w')
    for index in range(len(names)):
        f.write(names[index] + ' Basal Mass Loss\n')
        for t in range(num_time):
            f.write(str(massloss[index, t]) + '\n')
    f.close()

    

if __name__ == "__main__":

    spinup_path = raw_input("Path to control simulation directory: ")
    rcp_path = raw_input("Path to RCP simulation directory: ")
    timeseries_concat(spinup_path, rcp_path)
    
