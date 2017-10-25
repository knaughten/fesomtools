from netCDF4 import Dataset
from numpy import *
from fesom_grid import *
from monthly_avg import *

def rcp_seaice_extent_change ():

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/meshB/'
    directory_beg = '/short/y99/kaa561/FESOM/highres_spinup/'
    directories = ['/short/y99/kaa561/FESOM/rcp45_M/', '/short/y99/kaa561/FESOM/rcp45_A/', '/short/y99/kaa561/FESOM/rcp85_M/', '/short/y99/kaa561/FESOM/rcp85_A/', '/short/y99/kaa561/FESOM/highres_spinup/']
    file_beg = 'avg.ice.mean.1996.2005.nc'
    file_end = 'avg.ice.mean.2091.2100.nc'
    # Titles
    expt_names = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A', 'CONTROL']
    num_expts = len(directories)
    # Mesh parameters
    circumpolar = True
    cross_180 = False

    print 'Building mesh'
    elements = fesom_grid(mesh_path, circumpolar, cross_180)
    num_elm = len(elements)

    print 'Reading data'
    print '...1996-2005'
    # Calculate monthly averages for September
    aice_nodes_beg = monthly_avg(directory_beg + file_beg, 'area', 8)
    n2d = size(aice_nodes_beg)
    aice_nodes_end = empty([num_expts, n2d])
    for expt in range(num_expts):
        print '...' + expt_names[expt]
        aice_nodes_end[expt,:] = monthly_avg(directories[expt] + file_end, 'area', 8)

    print 'Calculating element-averages'
    aice_beg = empty(num_elm)
    aice_end = empty([num_expts, num_elm])
    # Also save area of each element
    area_elm = empty(num_elm)
    for i in range(num_elm):
        elm = elements[i]
        area_elm[i] = elm.area()
        aice_beg[i] = (aice_nodes_beg[elm.nodes[0].id] + aice_nodes_beg[elm.nodes[1].id] + aice_nodes_beg[elm.nodes[2].id])/3.0
        for expt in range(num_expts):
            aice_end[expt,i] = (aice_nodes_end[expt,elm.nodes[0].id] + aice_nodes_end[expt,elm.nodes[1].id] + aice_nodes_end[expt,elm.nodes[2].id])/3.0

    print 'Sea ice extent:'
    # 1996-2005
    # Select elements with concentration >= 15%
    flag_beg = aice_beg > 0.15
    # Integrate the area of these elements and convert to million km^2
    extent_beg = sum(flag_beg*area_elm)*1e-12
    print '1996-2005: ' + str(extent_beg) + ' million km^2'
    # 2091-2100
    flag_end = aice_end > 0.15
    for expt in range(num_expts):
        extent_end = sum(flag_end[expt,:]*area_elm)*1e-12
        percent_change = (extent_end - extent_beg)/extent_beg*100
        print expt_names[expt] + ': ' + str(extent_end) + ' million km^2; change of ' + str(percent_change) + '%'


# Command-line interface
if __name__ == "__main__":

    rcp_seaice_extent_change()
        
        
    
    

    
