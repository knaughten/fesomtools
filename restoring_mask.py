from numpy import *
from unrotate_grid import *

# Make a mask for the 2D FESOM mesh indicating which nodes should get surface
# salinity restoring (1) and which should be left alone (0). The current
# formulation is to have surface salinity restoring everywhere south of 30S
# except for the Antarctic continental shelf (defined as regions south of 60S
# with bathymetry shallower than 1500 m) and ice shelf cavities. Write an ASCII
# file with the mask, similar to cavity_flag_nod2d.out. To use this with FESOM
# switch on use_restoring_mask in namelist.oce. 
# Input:
# mesh_path = path to FESOM mesh directory, containing the files nod2d.out,
#             cavity_flag_nod2d.out, and depth.out. The output file
#             ssrestoring_mask.out will be saved in this directory.
def restoring_mask (mesh_path):

    # Parameters for mask
    # No restoring if lat > nbdry
    nbdry = -30
    # No restoring if lat < lat0 and depth < h0
    lat0 = -60       
    h0 = 1500

    # Read rotated latitude and longitude of 2D nodes
    f = open(mesh_path + 'nod2d.out', 'r')
    f.readline()
    rlon = []
    rlat = []
    for line in f:
        tmp = line.split()
        lon_tmp = float(tmp[1])
        if lon_tmp < -180:
            lon_tmp += 360
        elif lon_tmp > 180:
            lon_tmp -= 360
        rlon.append(lon_tmp)
        rlat.append(float(tmp[2]))
    f.close()
    rlon = array(rlon)
    rlat = array(rlat)
    # Unrotate grid
    glon, glat = unrotate_grid(rlon, rlat)

    # Read cavity flag
    f = open(mesh_path + 'cavity_flag_nod2d.out', 'r')
    cavity = []
    for line in f:
        cavity.append(int(line))
    f.close()
    cavity = array(cavity)

    # Read depth
    f = open(mesh_path + 'depth.out', 'r')
    depth = []
    for line in f:
        depth.append(-float(line))
    f.close()
    depth=array(depth)

    # Make sure all the arrays are the same size
    if size(glat) != size(cavity) or size(glon) != size(depth):
        print "Problem with array sizes"
        return

    # Start with restoring everywhere
    flag = ones(size(cavity))
    # No restoring in cavities
    index = cavity == 1
    flag[index] = 0.0
    # No restoring north of 30S
    index = glat > nbdry
    flag[index] = 0.0
    # No restoring south of 60S if bathymetry shallower than 1500 m
    index = (glat < lat0)*(depth < h0)
    flag[index] = 0.0

    # Write mask to file
    f = open(mesh_path + 'ssrestoring_mask.out', 'w')
    for elm in flag:
        f.write(str(int(elm)) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    restoring_mask(mesh_path)

    

    

    
