from netCDF4 import Dataset
from numpy import *

def calc_annual_ice_prod (model_dir, start_year, end_year, out_file):

    file_head = model_dir + 'MK44005.'
    file_tail = '.ice.diag.nc'
    sec_per_step = 5*24*60*60

    id = Dataset(file_head + str(start_year) + file_tail, 'r')
    n2d = id.variables['thdgr'].shape[1]
    id.close()
    ice_prod = zeros(n2d)

    for year in range(start_year, end_year+1):
        id = Dataset(file_head + str(year) + file_tail, 'r')
        thdgr = id.variables['thdgr'][:,:]
        id.close()
        index = thdgr < 0
        thdgr[index] = 0
        ice_prod += sum(thdgr, axis=0)*sec_per_step
    ice_prod = ice_prod/(end_year-start_year+1)

    id = Dataset(out_file, 'w')
    id.createDimension('nodes_2d', size(ice_prod))
    id.createVariable('ice_prod', 'f8', ('nodes_2d'))
    id.variables['ice_prod'].units = 'm/y'
    id.variables['ice_prod'][:] = ice_prod
    id.close()


if __name__ == "__main__":

    model_dir = input("Path to FESOM simulation output directory: ")
    start_year =  int(input("First year to process: "))
    end_year = int(input("Last year to process: "))
    out_file = input("Path to desired output file: ")
    calc_annual_ice_prod(model_dir, start_year, end_year, out_file)
    
