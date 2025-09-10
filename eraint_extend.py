from netCDF4 import Dataset
from numpy import *

# Given ERA-Interim files for 1992-2016, each containing atmospheric forcing
# data for multiple variables but one year at a time, unravel into FESOM
# forcing files, each containing one variable and one time of the day but
# all years.
def eraint_extend ():

    # Path to directory containing ERA-Interim subdaily files
    input_dir = '/short/m68/kaa561/metroms_iceshelf/data/originals/ERA_Interim/'
    # Path to directory for output files
    output_dir = '/short/y99/kaa561/FESOM/ERA_Interim_1992_2016/'
    # Variable names
    var_names_6h = ['sp', 't2m', 'd2m', 'u10', 'v10']
    # Variable units
    var_units_6h = ['Pa', 'K', 'K', 'm/s', 'm/s']
    # Beginning of output filename for each variable
    file_heads_6h = ['pair', 'tair', 'tdew', 'uwind', 'vwind']
    # End of output filename for each time of day
    file_tails_6h = ['_00.nc', '_06.nc', '_12.nc', '_18.nc']
    # Same for 12-hourly variables
    var_names_12h = ['tp', 'sf', 'e', 'ssrd', 'strd']
    var_units_12h = ['m/12h', 'm/12h', 'm/12h', 'J/m^2/12h', 'J/m^2/12h']
    file_heads_12h = ['precip', 'snow', 'evap', 'dswrf', 'dlwrf']
    file_tails_12h = ['_00_12.nc', '_12_12.nc']
    # Beginning of input filenames
    infile_head_6h = 'AN_'
    infile_head_12h_1 = 'FC_'
    infile_head_12h_2 = 'ER_'
    # End of input filenames
    infile_tail = '_subdaily_orig.nc'
    # Years to process
    year_start = 1992
    year_end = 2016

    # Read grid
    id = Dataset(input_dir + infile_head_6h + str(year_start) + infile_tail, 'r')
    lon = id.variables['longitude'][:]
    lat = id.variables['latitude'][:]
    id.close()

    # Loop over 6-hourly variables
    for i in range(len(var_names_6h)):
        var = var_names_6h[i]
        print('Processing variable ' + var)
        # Loop over time of day
        for j in range(4):
            # Construct output filename
            file_name = output_dir + file_heads_6h[i] + file_tails_6h[j]
            print('Setting up ' + file_name)
            o_id = Dataset(file_name, 'w')
            o_id.createDimension('longitude', size(lon))
            o_id.createDimension('latitude', size(lat))
            o_id.createDimension('time', None)
            o_id.createVariable('longitude', 'f8', ('longitude'))
            o_id.variables['longitude'].units = 'degrees'
            o_id.variables['longitude'][:] = lon
            o_id.createVariable('latitude', 'f8', ('latitude'))
            o_id.variables['latitude'].units = 'degrees'
            o_id.variables['latitude'][:] = lat
            o_id.createVariable('time', 'f8', ('time'))
            o_id.variables['time'].units = 'hours since 1900-1-1 00:00:00'
            o_id.variables['time'].calendar = 'standard'
            o_id.createVariable(var, 'f8', ('time', 'latitude', 'longitude'))
            o_id.variables[var].units = var_units_6h[i]
            t_posn = 0  # Day of simulation
            # Loop over years
            for year in range(year_start, year_end+1):
                print('Year ' + str(year))
                # Open input file for this year
                infile_name = input_dir + infile_head_6h + str(year) + infile_tail
                id = Dataset(infile_name, 'r')
                # Figure out number of days in this year
                num_days = 365
                if year % 4 == 0:
                    num_days = 366
                # Loop over days
                for t in range(num_days):
                    # Select correct index and save to output file
                    time = id.variables['time'][4*t+j]
                    o_id.variables['time'][t_posn] = time
                    data = id.variables[var][4*t+j,:,:]
                    o_id.variables[var][t_posn,:,:] = data
                    t_posn += 1
            o_id.close()

    # Loop over 12-hourly variables
    for i in range(len(var_names_12h)):
        var = var_names_12h[i]
        print('Processing variable ' + var)
        # Loop over time of day
        for j in range(2):
            # Construct output filename
            file_name = output_dir + file_heads_12h[i] + file_tails_12h[j]
            print('Setting up ' + file_name)
            o_id = Dataset(file_name, 'w')
            o_id.createDimension('longitude', size(lon))
            o_id.createDimension('latitude', size(lat))
            o_id.createDimension('time', None)
            o_id.createVariable('longitude', 'f8', ('longitude'))
            o_id.variables['longitude'].units = 'degrees'
            o_id.variables['longitude'][:] = lon
            o_id.createVariable('latitude', 'f8', ('latitude'))
            o_id.variables['latitude'].units = 'degrees'
            o_id.variables['latitude'][:] = lat
            o_id.createVariable('time', 'f8', ('time'))
            o_id.variables['time'].units = 'hours since 1900-1-1 00:00:00'
            o_id.variables['time'].calendar = 'standard'
            o_id.createVariable(var, 'f8', ('time', 'latitude', 'longitude'))
            o_id.variables[var].units = var_units_12h[i]
            t_posn = 0  # Day of simulation
            # Loop over years
            for year in range(year_start, year_end+1):
                print('Year ' + str(year))
                # Open input file for this year
                # First figure out which file it's in
                if var in ['tp', 'sf']:
                    # Always in FC file
                    infile_head_12h = infile_head_12h_1
                elif var in ['ssrd', 'strd']:
                    # Always in ER file
                    infile_head_12h = infile_head_12h_2
                elif var == 'e':
                    # In ER file from 1992-2005, FC file from 2006-2016
                    if year <= 2005:
                        infile_head_12h = infile_head_12h_2
                    else:
                        infile_head_12h = infile_head_12h_1
                infile_name = input_dir + infile_head_12h + str(year) + infile_tail
                id = Dataset(infile_name, 'r')
                # Figure out number of days in this year
                num_days = 365
                if year % 4 == 0:
                    num_days = 366
                # Loop over days
                for t in range(num_days):
                    # Select correct index and save to output file
                    time = id.variables['time'][2*t+j]
                    o_id.variables['time'][t_posn] = time
                    data = id.variables[var][2*t+j,:,:]
                    o_id.variables[var][t_posn,:,:] = data
                    t_posn += 1
            o_id.close()
    

# Command-line interface
if __name__ == "__main__":

    eraint_extend()
