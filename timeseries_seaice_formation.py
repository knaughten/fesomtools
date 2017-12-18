from netCDF4 import Dataset
from numpy import *
from fesom_grid import *

def timeseries_seaice_formation (mesh_path, output_path, start_year, end_year, log_file):

    # Naming conventions for FESOM output files
    file_head = output_path + 'MK44005.'
    file_tail = '.ice.diag.nc'
    num_years = end_year - start_year + 1
    # Parameters for selecting continental shelf
    lat0 = -60
    h0 = 1500
    # Seconds to years conversion
    sec_per_year = 365.25*24*60*60

    print 'Building mesh'
    elements = fesom_grid(mesh_path, circumpolar=True, cross_180=True)

    print 'Selecting continental shelf'
    # Set up an array of area of each element, zero if it's not on the
    # continental shelf
    shelf_areas = zeros(len(elements))
    for i in range(len(elements)):
        elm = elements[i]
        lat = mean(elm.lat)
        bathy = mean(array([(elm.nodes[0].find_bottom()).depth, (elm.nodes[1].find_bottom()).depth, (elm.nodes[2].find_bottom()).depth]))
        if lat < lat0 and bathy < h0 and not elm.cavity:
            shelf_areas[i] = elm.area()

    # Set up array for net sea ice formation on continental shelf
    formation = zeros(num_years)
    for year in range(start_year, end_year+1):
        print 'Processing year ' + str(year)
        id = Dataset(file_head + str(year) + file_tail, 'r')
        # Read thdgr, annually average, and convert from m/s to m/y
        thdgr = mean(id.variables['thdgr'][:,:], axis=0)*sec_per_year
        id.close()
        # Average over elements
        thdgr_elm = zeros(len(elements))
        for i in range(len(elements)):
            elm = elements[i]
            thdgr_elm[i] = mean(array([thdgr[elm.nodes[0].id], thdgr[elm.nodes[1].id], thdgr[elm.nodes[2].id]]))
        # Integrate and convert to thousand km^3/y
        formation[year-start_year] = sum(thdgr_elm*shelf_areas)*1e-12

    print 'Saving results to log file'
    f = open(log_file, 'w')
    f.write('Net sea ice formation on continental shelf (thousand km^3/y):\n')
    for t in range(num_years):
        f.write(str(formation[t]) + '\n')
    f.close()
        

# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    output_path = raw_input("Path to FESOM output directory: ")
    start_year = int(raw_input("First year to process: "))
    end_year = int(raw_input("Last year to process: "))
    log_file = raw_input("Path to logfile: ")
    timeseries_seaice_formation(mesh_path, output_path, start_year, end_year, log_file)
