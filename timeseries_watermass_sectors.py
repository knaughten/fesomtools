from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from fesom_grid import *

def timeseries_watermass_sectors (mesh_path, output_path, start_year, end_year, log_file, fig_dir=''):

    # Titles and figure names for each sector
    sector_names = ['Filchner-Ronne Ice Shelf Cavity', 'Eastern Weddell Region Cavities', 'Amery Ice Shelf Cavity', 'Australian Sector Cavities', 'Ross Sea Cavities', 'Amundsen Sea Cavities', 'Bellingshausen Sea Cavities', 'Larsen Ice Shelf Cavities', 'All Ice Shelf Cavities']
    fig_names = ['filchner_ronne_watermass.png', 'eweddell_watermass.png', 'amery_watermass.png', 'australian_watermass.png', 'ross_watermass.png', 'amundsen_watermass.png', 'bellingshausen_watermass.png', 'larsen_watermass.png', 'total_antarctica_watermass.png']
    num_sectors = len(sector_names)
    # Water masses to consider
    wm_names = ['ISW', 'HSSW', 'LSSW', 'AASW', 'MCDW', 'CDW']
    num_watermasses = len(wm_names)
    wm_colours = ['cyan', 'black', 'blue', 'green', 'magenta', 'red']
    # Only consider elements south of 30S
    circumpolar = True
    # Don't make second copies of elements that cross 180E
    cross_180 = False
    # Naming conventions for FESOM output files
    file_head = output_path + 'MK44005.'
    file_tail = '.oce.mean.nc'
    num_years = end_year - start_year + 1

    prev_years = 0
    # Check if the log file exists
    if exists(log_file):
        print('Reading previously calculated values')
        # First just figure out how many years are in the log file
        f = open(log_file, 'r')
        f.readline()
        for line in f:
            try:
                tmp = float(line)
                prev_years += 1
            except(ValueError):
                break
        f.close()
        # Now set up array of water mass proportions in each sector
        percent_watermass = empty([num_watermasses, num_sectors, prev_years+num_years])
        # Fill the first prev_years
        f = open(log_file, 'r')
        f.readline()
        wm_key = 0
        while wm_key < num_watermasses:
            sector = 0
            while sector < num_sectors:
                year = 0
                for line in f:
                    try:
                        percent_watermass[wm_key, sector, year] = float(line)
                        year += 1
                    except(ValueError):
                        break
                sector += 1
            wm_key += 1
        f.close()        
    else:
        # Set up empty array for water mass proportions
        percent_watermass = empty([num_watermasses, num_sectors, num_years])

    print('Building grid')
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

    print('Calculating water mass breakdown')    
    # Loop over years
    for year in range(start_year, end_year+1):
        print('Processing year ' + str(year))
        # Initialise volume of each water mass in each sector
        vol_watermass = zeros([num_watermasses, num_sectors])
        # Read temperature and salinity for this year, annually average
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
                    # Calculate average temperature, salinity, and
                    # layer thickness for this 3D triangular prism
                    temp_vals = []
                    salt_vals = []
                    dz_vals = []
                    for n in range(3):
                        temp_vals.append(temp[nodes[n].id])
                        salt_vals.append(salt[nodes[n].id])
                        temp_vals.append(temp[nodes[n].below.id])
                        salt_vals.append(salt[nodes[n].below.id])
                        dz_vals.append(abs(nodes[n].depth - nodes[n].below.depth))
                        # Get ready for next iteration of loop
                        nodes[n] = nodes[n].below
                    curr_temp = mean(array(temp_vals))
                    curr_salt = mean(array(salt_vals))
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
                    # Integrate its volume for sector(s) the element is in
                    curr_sectors = 0
                    for sector in range(num_sectors):
                        if location_flag[sector,i] == 1:
                            curr_sectors += 1
                            vol_watermass[wm_key, sector] += curr_volume
                    # Should be in exactly 2 sectors (1 + total Antarctica)
                    if curr_sectors != 2:
                        print('Wrong number of sectors for element ' + str(i))
        if year==start_year:
            # Find the total volume of each sector by adding up the volume
            # of each water mass. Only need to do this once because shouldn't
            # change over time.
            vol_sectors = sum(vol_watermass, axis=0)
        # Calculate percentage of each water mass in each sector
        for wm_key in range(num_watermasses):
            for sector in range(num_sectors):
                percent_watermass[wm_key, sector, year-start_year+prev_years] = vol_watermass[wm_key, sector]/vol_sectors[sector]*100

    # Make time axis
    time = list(range(start_year-prev_years, end_year+1))

    print('Plotting')
    # One plot for each sector
    for sector in range(num_sectors):
        fig = figure()
        ax = fig.add_subplot(1,1,1)
        # Loop over water masses
        for wm_key in range(num_watermasses):
            plot(time, percent_watermass[wm_key, sector, :], color=wm_colours[wm_key], label=wm_names[wm_key], linewidth=2)
        xlabel('year')
        ylabel('percent volume')
        xlim([start_year-prev_years, end_year])
        title(sector_names[sector])
        grid(True)
        # Move plot over to make room for legend
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        # Make legend
        ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
        fig.savefig(fig_dir + fig_names[sector])

    print('Saving results to log file')
    f = open(log_file, 'w')
    for wm_key in range(num_watermasses):
        for sector in range(num_sectors):
            f.write(wm_names[wm_key] + 'in ' + sector_names[sector] + '(%)\n')
            for t in range(prev_years+num_years):
                f.write(str(percent_watermass[wm_key, sector, t]) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    mesh_path = input("Path to FESOM mesh directory: ")
    output_path = input("Path to FESOM output directory: ")
    start_year = int(input("First year to process: "))
    end_year = int(input("Last year to process: "))
    log_file = input("Path to logfile: ")
    fig_dir = input("Path to directories to save figures: ")
    timeseries_watermass_sectors(mesh_path, output_path, start_year, end_year, log_file, fig_dir)
    
    
                
                            
        
        
    
