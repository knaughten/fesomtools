from numpy import *
from netCDF4 import Dataset
from fesom_grid import *

def massloss_percent_winds ():

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/meshB/'
    directory_winds = '/short/y99/kaa561/FESOM/rcp85_M/'
    directory_nowinds = '/short/y99/kaa561/FESOM/rcp85_M_no_wind_anom/'
    file_name = 'annual_avg.forcing.diag.2091.2100.nc'
    # Seconds per year
    sec_per_year = 365.25*24*3600
    # Density of ice in kg/m^3
    rho_ice = 916
    # Sectors to split Antarctica into
    sector_names = ['Filchner-Ronne Ice Shelf', 'Eastern Weddell Region', 'Amery Ice Shelf', 'Australian Sector', 'Ross Sea', 'Amundsen Sea', 'Bellingshausen Sea', 'Larsen Ice Shelves']
    # Number of sectors
    num_sectors = len(sector_names)

    print 'Building mesh'
    elements = fesom_grid(mesh_path, circumpolar=True, cross_180=False)

    print 'Reading data'
    id = Dataset(directory_winds + file_name, 'r')
    ismr_nodes_winds = id.variables['wnet'][0,:]*sec_per_year
    id.close()
    id = Dataset(directory_nowinds + file_name, 'r')
    ismr_nodes_nowinds = id.variables['wnet'][0,:]*sec_per_year
    id.close()
    # Average over elements in ice shelf cavities
    ismr_elm_winds = []
    ismr_elm_nowinds = []
    for elm in elements:
        if elm.cavity:
            ismr_elm_winds.append(mean([ismr_nodes_winds[elm.nodes[0].id], ismr_nodes_winds[elm.nodes[1].id], ismr_nodes_winds[elm.nodes[2].id]]))
            ismr_elm_nowinds.append(mean([ismr_nodes_nowinds[elm.nodes[0].id], ismr_nodes_nowinds[elm.nodes[1].id], ismr_nodes_nowinds[elm.nodes[2].id]]))
    ismr_elm_winds = array(ismr_elm_winds)
    ismr_elm_nowinds = array(ismr_elm_nowinds)

    print 'Integrating mass loss'
    total_massloss_winds = zeros(num_sectors)
    total_massloss_nowinds = zeros(num_sectors)
    # Loop over elements
    i = 0
    for elm in elements:
        if elm.cavity:
            # Figure out which sector this ice shelf element falls into
            # First get average lon and lat across 3 Nodes
            lon = mean(elm.lon)
            lat = mean(elm.lat)
            if lon >= -85 and lon < -30 and lat < -74:
                # Filchner-Ronne
                index = 0            
            elif lon >= -30 and lon < 65:
                # Eastern Weddell region
                index = 1            
            elif lon >= 65 and lon < 76:
                # Amery
                index = 2            
            elif lon >= 76 and lon < 165 and lat >= -74:
                # Australian sector
                index = 3            
            elif (lon >= 155 and lon < 165 and lat < -74) or (lon >= 165) or (lon < -140):
                # Ross Sea
                index = 4            
            elif (lon >= -140 and lon < -105) or (lon >= -105 and lon < -98 and lat < -73.1):
                # Amundsen Sea
                index = 5            
            elif (lon >= -104 and lon < -98 and lat >= -73.1) or (lon >= -98 and lon < -66 and lat >= -75):
                # Bellingshausen Sea
                index = 6            
            elif lon >= -66 and lon < -59 and lat >= -74:
                # Larsen Ice Shelves
                index = 7
            else:
                print 'No region found for lon=',str(lon),', lat=',str(lat)
                break #return
            # Integrate total mass loss in this sector
            total_massloss_winds[index] += ismr_elm_winds[i]*elm.area()*rho_ice*1e-12
            total_massloss_nowinds[index] += ismr_elm_nowinds[i]*elm.area()*rho_ice*1e-12
            i += 1

    # Calculate change in mass loss due to removing winds
    for index in range(num_sectors):
        massloss_winds = total_massloss_winds[index]
        massloss_nowinds = total_massloss_nowinds[index]
        percent_change = (massloss_nowinds - massloss_winds)/massloss_winds*1e2
        print sector_names[index] + ': ' + str(percent_change) + '%'
    # Total Antarctica
    massloss_winds = sum(total_massloss_winds)
    massloss_nowinds = sum(total_massloss_nowinds)
    percent_change = (massloss_nowinds - massloss_winds)/massloss_winds*1e2
    print 'Total Antarctica: ' + str(percent_change) + '%'


# Command-line interface
if __name__ == "__main__":

    massloss_percent_winds()
