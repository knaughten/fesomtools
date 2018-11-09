from netCDF4 import Dataset
from unrotate_grid import *

def format_output_nick (model_dir, start_year, end_year, output_head):

    mesh_path = '/short/y99/kaa561/FESOM/mesh/meshB/'
    expt_name = 'MK44005.'
    melt_tail = '.forcing.diag.nc'
    temp_tail = '.oce.mean.nc'

    print 'Preparing mesh'
    # Read rotated lon and lat
    rlon = []
    rlat = []
    f = open(mesh_path + 'nod2d.out', 'r')
    n2d = int(f.readline())
    for line in f:
        tmp = line.split()
        lon_tmp = float(tmp[1])
        lat_tmp = float(tmp[2])
        # Make sure longitude is in the range [-180, 180]
        if lon_tmp < -180:
            lon_tmp += 360
        elif lon_tmp > 180:
            lon_tmp -= 360
        rlon.append(lon_tmp)
        rlat.append(lat_tmp)
    f.close()
    rlon = array(rlon)
    rlat = array(rlat)
    # Unrotate
    lon, lat = unrotate_grid(rlon, rlat)
    # Read cavity flag
    f = open(mesh_path + 'cavity_flag_nod2d.out', 'r')
    cavity = []
    for line in f:
        tmp = int(line)
        if tmp == 1:
            cavity.append(True)
        elif tmp == 0:
            cavity.append(False)
        else:
            print 'Problem'
    f.close()

    # Loop over years
    for year in range(start_year, end_year+1):
        print 'Processing ' + str(year)
        # Read melt rate (annually averaged) and convert to m/y
        id = Dataset(model_dir + expt_name + str(year) + melt_tail, 'r')
        melt = mean(id.variables['wnet'][:,:], axis=0)*24*60*60*365.25
        id.close()
        # Read first n2d nodes of temperature (annually averaged)
        id = Dataset(model_dir + expt_name + str(year) + temp_tail, 'r')
        temp = mean(id.variables['temp'][:,:n2d], axis=0)
        id.close()
        # Set up output file
        f = open(output_head + str(year), 'w')
        # Loop over nodes
        for n in range(n2d):
            # Only care about cavity nodes
            if cavity[n]:
                # Longitude
                f.write('{0:16}'.format(str(round(lon[n],8))))
                # Latitude
                f.write('{0:16}'.format(str(round(lat[n],8))))
                # Melt rate
                f.write('{0:16}'.format(str(round(melt[n],8))))
                # Temperature
                f.write('{0:16}'.format(str(round(temp[n],8))))
                # Finished
                f.write('\n')
        f.close()


# Command-line interface
if __name__ == "__main__":

    model_dir = raw_input("Path to FESOM output directory: ")
    start_year = int(raw_input("First year to process: "))
    end_year = int(raw_input("Last year to process: "))
    output_head = raw_input("Beginning of output text file name: ")
    format_output_nick(model_dir, start_year, end_year, output_head)
        
    
