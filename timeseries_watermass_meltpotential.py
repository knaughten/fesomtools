from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from fesom_grid import *
from unesco import *

def timeseries_watermass_meltpotential (mesh_path, output_path, start_year, end_year, log_file):

    # Titles for each sector
    sector_names = ['Filchner-Ronne Ice Shelf Cavity', 'Eastern Weddell Region Cavities', 'Amery Ice Shelf Cavity', 'Australian Sector Cavities', 'Ross Sea Cavities', 'Amundsen Sea Cavities', 'Bellingshausen Sea Cavities', 'Larsen Ice Shelf Cavities', 'All Ice Shelf Cavities']
    num_sectors = len(sector_names)
    # Water masses to consider
    wm_names = ['ISW', 'HSSW', 'LSSW', 'AASW', 'MCDW', 'CDW']
    num_watermasses = len(wm_names)
    # Only consider elements south of 30S
    circumpolar = True
    # Don't make second copies of elements that cross 180E
    cross_180 = False
    # Naming conventions for FESOM output files
    file_head = output_path + 'MK44005.'
    file_tail = '.oce.mean.nc'
    num_years = end_year - start_year + 1
    # Specific heat of seawater (J/K/kg)
    cpw = 4180
    # Coefficients for in-situ freezing point calculation
    a = -0.0575  # Salinity dependence (K/psu)
    b = 0.0901   # Surface freezing point at 0 salinity (C)
    c = 7.61e-4  # Depth dependence (K/m)

    print('Building mesh')
    elements = fesom_grid(mesh_path, circumpolar, cross_180)

    print('Categorising elements into sectors')
    location_flag = zeros([num_sectors, len(elements)])
    for i in range(len(elements)):
        elm = elements[i]
        # Make sure we're actually in an ice shelf cavity
        if elm.cavity:
            # Figure out which sector this ice shelf element falls into
            lon = mean(elm.lon)
            lat = mean(elm.lat)
            if lon >= -85 and lon < -30 and lat < -74:
                # Filchner-Ronne
                location_flag[0,i] = 1           
            elif lon >= -30 and lon < 65:
                # Eastern Weddell region
                location_flag[1,i] = 1       
            elif lon >= 65 and lon < 76:
                # Amery
                location_flag[2,i] = 1       
            elif lon >= 76 and lon < 165 and lat >= -74:
                # Australian sector
                location_flag[3,i] = 1       
            elif (lon >= 155 and lon < 165 and lat < -74) or (lon >= 165) or (lon < -140):
                # Ross Sea
                location_flag[4,i] = 1         
            elif (lon >= -140 and lon < -105) or (lon >= -105 and lon < -98 and lat < -73.1):
                # Amundsen Sea
                location_flag[5,i] = 1          
            elif (lon >= -104 and lon < -98 and lat >= -73.1) or (lon >= -98 and lon < -66 and lat >= -75):
                # Bellingshausen Sea
                location_flag[6,i] = 1
            elif lon >= -66 and lon < -59 and lat >= -74:
                # Larsen Ice Shelves
                location_flag[7,i] = 1
            else:
                print('No region found for lon=',str(lon),', lat=',str(lat))
                break #return
            # All ice shelf elements are in Total Antarctica
            location_flag[8,i] = 1        

    print('Calculating melt potential')
    mp = zeros([num_watermasses, num_sectors, num_years])
    for year in range(start_year, end_year+1):
        print('Processing ' + str(year))
        # Read temperature and salinity
        id = Dataset(file_head + str(year) + file_tail, 'r')
        temp = mean(id.variables['temp'][:,:], axis=0)
        salt = mean(id.variables['salt'][:,:], axis=0)
        id.close()
        # Loop over elements
        for i in range(len(elements)):
            elm = elements[i]
            # Check if we're in an ice shelf cavity
            if elm.cavity:
                # Get area of 2D element
                area = elm.area()
                nodes = [elm.nodes[0], elm.nodes[1], elm.nodes[2]]
                # Loop downward
                while True:
                    if nodes[0].below is None or nodes[1].below is None or nodes[2].below is None:
                        # Reached the bottom
                        break
                    # Calculate average temperature, salinity, depth, and
                    # layer thickness for this 3D triangular prism
                    temp_vals = []
                    salt_vals = []
                    z_vals = []
                    dz_vals = []
                    for n in range(3):
                        temp_vals.append(temp[nodes[n].id])
                        salt_vals.append(salt[nodes[n].id])
                        z_vals.append(nodes[n].depth)
                        temp_vals.append(temp[nodes[n].below.id])
                        salt_vals.append(salt[nodes[n].below.id])
                        z_vals.append(nodes[n].below.depth)
                        dz_vals.append(abs(nodes[n].depth - nodes[n].below.depth))
                        # Get ready for next iteration of loop
                        nodes[n] = nodes[n].below
                    curr_temp = mean(array(temp_vals))
                    curr_salt = mean(array(salt_vals))
                    curr_z = mean(array(z_vals))
                    curr_volume = area*mean(array(dz_vals))
                    # Get surface freezing point at this salinity
                    curr_tfrz = -0.0575*curr_salt + 1.7105e-3*sqrt(curr_salt**3) - 2.155e-4*curr_salt**2
                    # Figure out what water mass this is
                    if curr_temp < curr_tfrz:
                        # ISW
                        wm_key = 0
                    elif curr_salt < 34:
                        # AASW
                        wm_key = 3
                    elif curr_temp > 0:
                        # CDW
                        wm_key = 5
                    elif curr_temp > -1.5:
                        # MCDW
                        wm_key = 4
                    elif curr_salt < 34.5:
                        # LSSW
                        wm_key = 2
                    else:
                        # HSSW
                        wm_key = 1
                    # Integrate melt potential
                    # First need (potential) density
                    curr_rho = unesco(curr_temp, curr_salt, 0)
                    # And in-situ freezing point
                    curr_tfrz_insitu = a*curr_salt + b + c*(-1*curr_z)
                    curr_sectors = 0
                    for sector in range(num_sectors):
                        if location_flag[sector,i] == 1:
                            curr_sectors += 1
                            mp[wm_key, sector, year-start_year] += (curr_temp-curr_tfrz_insitu)*curr_volume*cpw*curr_rho
                    # Should be in exactly 2 sectors (1 + total Antarctica)
                    if curr_sectors != 2:
                        print('Wrong number of sectors for element ' + str(i))

    print('Saving results to log file')
    f = open(log_file, 'w')
    for wm_key in range(num_watermasses):
        for sector in range(num_sectors):
            f.write('Melt potential of ' + wm_names[wm_key] + ' in ' + sector_names[sector] + '(J)\n')
            for t in range(num_years):
                f.write(str(mp[wm_key, sector, t]) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    mesh_path = input("Path to FESOM mesh directory: ")
    output_path = input("Path to FESOM output directory: ")
    start_year = int(input("First year to process: "))
    end_year = int(input("Last year to process: "))
    log_file = input("Path to logfile: ")
    timeseries_watermass_meltpotential(mesh_path, output_path, start_year, end_year, log_file)
    
    
                
                            
        
        
    
