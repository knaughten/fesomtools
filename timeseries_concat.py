from numpy import *

# For the logfiles from timeseries_dpt.py, timeseries_seaice.py,
# timeseries_seaice_extent.py, timeseries_massloss.py,
# timeseries_massloss_sectors.py, and timeseries_amundsen.py, append the third
# repetition of the forcing in the control simulation (1992-2005) to the
# beginning of the RCP (2006-2100).
# Input:
# spinup_path, rcp_path = paths to experiment directories for the control
#                         simulation and RCP, each containing the following
#                         logfiles: dpt.log, seaice.log, seaice_extent.log,
#                         massloss.log, massloss_sectors.log, amundsen.log.
def timeseries_concat (spinup_path, rcp_path):

    # Name of each ice shelf cavity
    names = ['All Ice Shelves', 'Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    # Name of each ice shelf sector
    sector_names = ['Filchner-Ronne Ice Shelf', 'Eastern Weddell Region', 'Amery Ice Shelf', 'Australian Sector', 'Ross Sea', 'Amundsen Sea', 'Bellingshausen Sea', 'Larsen Ice Shelves', 'Total Antarctica']
    # Forcing years for spinup
    year_start = 1992
    year_end = 2005
    # Skip the first 2 repetitions of the control simulation
    skipyears = 28
    peryear = 365/5
    numyears = year_end - year_start + 1

    # Drake Passage transport
    # Read logfile for control simulation
    dpt = []    
    f = open(spinup_path + 'dpt.log', 'r')
    f.readline()
    for line in f:
        dpt.append(float(line))
    f.close()
    # Select the third repetition
    dpt = dpt[skipyears*peryear:(skipyears+numyears)*peryear]
    num_time_spinup = len(dpt)
    # Read logfile for RCP simulation and append these values to existing array
    f = open(rcp_path + 'dpt.log', 'r')
    f.readline()
    for line in f:
        dpt.append(float(line))
    f.close()
    num_time = len(dpt)
    # Write a new logfile for the RCP simulation
    f = open(rcp_path + 'dpt.log', 'w')
    f.write('Drake Passage Transport (Sv):\n')
    for elm in dpt:
        f.write(str(elm) + '\n')
    f.close()

    # Repeat for sea ice area and volume
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

    # Repeat for sea ice extent
    seaice_extent = []
    f = open(spinup_path + 'seaice_extent.log', 'r')
    f.readline()
    for line in f:
        seaice_extent.append(float(line))
    f.close()
    seaice_extent = seaice_extent[skipyears*peryear:(skipyears+numyears)*peryear]
    f = open(rcp_path + 'seaice_extent.log', 'r')
    f.readline()
    for line in f:
        seaice_extent.append(float(line))
    f.close()
    f = open(rcp_path + 'seaice_extent.log', 'w')
    f.write('Sea Ice Extent (million km^2):\n')
    for elm in seaice_extent:
        f.write(str(elm) + '\n')
    f.close()

    # Repeat for ice shelf mass loss for each ice shelf
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

    # Repeat for ice shelf mass loss for each sector
    massloss_sectors = empty([len(sector_names), num_time])
    f = open(spinup_path + 'massloss_sectors.log', 'r')
    f.readline()
    for index in range(len(sector_names)):
        massloss_tmp = []
        for line in f:
            try:
                massloss_tmp.append(float(line))
            except(ValueError):
                break
        massloss_tmp = massloss_tmp[skipyears*peryear:(skipyears+numyears)*peryear]
        massloss_sectors[index, 0:num_time_spinup] = massloss_tmp
    f.close()
    f = open(rcp_path + 'massloss_sectors.log', 'r')
    f.readline()
    for index in range(len(sector_names)):
        massloss_tmp = []
        for line in f:
            try:
                massloss_tmp.append(float(line))
            except(ValueError):
                break
        massloss_sectors[index, num_time_spinup:] = massloss_tmp
    f.close()
    f = open(rcp_path + 'massloss_sectors.log', 'w')
    for index in range(len(sector_names)):
        f.write(names[index] + ' Basal Mass Loss (Gt/y)\n')
        for t in range(num_time):
            f.write(str(massloss_sectors[index, t]) + '\n')
    f.close()

    # Repeat for Amundsen Sea ice-to-ocean freshwater flux
    ice2ocn = []    
    f = open(spinup_path + 'amundsen.log', 'r')
    f.readline()
    for line in f:
        ice2ocn.append(float(line))
    f.close()
    # Select the third repetition
    ice2ocn = ice2ocn[skipyears*peryear:(skipyears+numyears)*peryear]
    f = open(rcp_path + 'amundsen.log', 'r')
    f.readline()
    for line in f:
        ice2ocn.append(float(line))
    f.close()
    f = open(rcp_path + 'amundsen.log', 'w')
    f.write('Average ice-to-ocean freshwater flux (1e-8 m/s):\n')
    for elm in ice2ocn:
        f.write(str(elm) + '\n')
    f.close()

    
# Command-line interface
if __name__ == "__main__":

    spinup_path = raw_input("Path to control simulation directory: ")
    rcp_path = raw_input("Path to RCP simulation directory: ")
    timeseries_concat(spinup_path, rcp_path)
    
