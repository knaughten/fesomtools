from netCDF4 import Dataset
from numpy import *

# Based on find_shelf.m by Ralph Timmermann
def shelf ():

    # Paths to RTopo files
    rtopo_shelf_file = '../RTopo105/RTopo105_data.nc' # Contains draft
    rtopo_imask_file = '../RTopo105/RTopo105_aux.nc' # Contains amask
    # Paths to existing 2D mesh files
    nod2d_file = 'nod2d.out'
    elem2d_file = 'elem2d.out'
    # Desired path to output file, containing ice shelf draft at every 2D node
    output_file = 'shelf.out'
    # Ice shelf draft cutoff
    shelf_min = -10
    # Latitude indicies which contain all the ice shelf points in RTopo
    shelf_sbdry = 300
    shelf_nbdry = 1600

    # Read grid and shelf (negative for ice shelves, positive for grounded ice)
    print 'Reading RTopo data'
    id = Dataset(rtopo_shelf_file, 'r')
    lon = id.variables['lon'][:]
    lat = id.variables['lat'][:]
    shelf = id.variables['draft'][:,:]
    id.close()
    num_lon = size(lon)
    num_lat = size(lat)
    # Read mask (0 = open ocean, 1 = grounded ice, 2 = ice shelf, 3 = bare rock)
    id = Dataset(rtopo_imask_file, 'r')
    amask = id.variables['amask'][:]
    id.close()

    # Read nod2d.out
    print 'Reading 2D nodes'
    file = open(nod2d_file, 'r')
    num_nod2d = int(file.readline())
    node_lon = zeros(num_nod2d)
    node_lat = zeros(num_nod2d)
    node_coast = zeros(num_nod2d)
    posn = 0
    for line in file:
        tmp = line.split()
        node_lon[posn] = float(tmp[1])
        node_lat[posn] = float(tmp[2])
        node_coast[posn] = int(tmp[3])
        posn += 1
    file.close()

    # Read elem2d.out
    print 'Reading 2D elements'
    file = open(elem2d_file, 'r')
    num_elem2d = int(file.readline())
    elem2d = zeros([num_elem2d, 3])
    posn = 0
    for line in file:
        tmp = line.split()
        elem2d[posn,0] = float(tmp[0])
        elem2d[posn,1] = float(tmp[1])
        elem2d[posn,2] = float(tmp[2])
        posn +=1
    file.close()

    # Make working mask and shelf arrays
    amask_tmp = copy(amask)
    shelf_tmp = copy(shelf)

    # Extend ice shelf mask into grounded ice, repeating 5 times
    for rep in range(5):
        print 'Extending ice shelf into grounded ice: repetition ' + str(rep+1) + ' of 5'
        for j in range(shelf_sbdry, shelf_nbdry):
            for i in range(num_lon):
                # Continue if this is a grounded ice point, or a very shallow
                # ice shelf point
                if amask[j,i] == 1 or shelf[j,i] > shelf_min:
                    # Look at neighbours 1 cell away in any direction
                    amask_cut = amask[j-1:j+1,i-1:i+1]
                    shelf_cut = shelf[j-1:j+1,i-1:i+1]
                    # Continue if there is no open water and at least one
                    # ice shelf point
                    if all(amask_cut!=0) and any((amask_cut==2)*(shelf_cut<shelf_min)):
                        print 'Modifying index ('+str(j)+','+str(i)+')'
                        # Set this as an ice shelf point in the new mask
                        amask_tmp[j,i] = 2
                        # Set the ice shelf draft to the average of the
                        # ice shelf neighbour points
                        index = nonzero((amask_cut==2)*(shelf_cut<shelf_min))
                        shelf_tmp[j,i] = mean(shelf_cut[index])
        # Save the changes made to amask and shelf
        amask = copy(amask_tmp)
        shelf = copy(shelf_tmp)

    # Fill zero thickness gaps
    print 'Filling gaps of zero thickness in ice shelf'
    for j in range(shelf_sbdry, shelf_nbdry):
        for i in range(num_lon):
            # Continue if shelf > shelf_min (either grounded ice or very shallow
            # ice shelf)
            if shelf[j,i] > shelf_min:
                # Look at neighbours 1 cell away in any direction
                shelf_cut = shelf[j-1:j+1,i-1:i+1]
                if any(shelf_cut<0):
                    print 'Modifying index ('+str(j)+','+str(i)+')'
                    # If there are any ice shelf neighbour points (shelf<0),
                    # set the ice shelf draft to the average of these
                    index = nonzero(shelf_cut<0)
                    shelf_tmp[j,i] = mean(shelf_cut[index])
                else:
                    # If there are no ice shelf neighbour points, set the
                    # ice shelf draft to zero
                    shelf_tmp[j,i] = 0.0
    # Save the changes made to shelf
    shelf = copy(shelf_tmp)

    # Make a binary ice shelf mask
    imask = zeros(shape(amask))
    imask[amask==2] = 1

    # Set up arrays to store ice shelf draft and ice shelf mask at each 2d node
    node_shelf = zeros(num_nod2d)
    node_imask = zeros(num_nod2d)

    # Interpolate from RTopo to 2D nodes
    print 'Interpolating to 2D nodes'
    for n in range(num_nod2d):

        # Find the RTopo indices this node falls between
        jN = nonzero(lat > node_lat[n])[0][0] # First latitude index to the north
        jS = jN - 1 # Last latitude index to the south
        iE = nonzero(lon > node_lon[n])[0][0] # First longitude index to the east
        iW = iE - 1 # Last longitude index to the west

        # Bilinear interpolation in two steps
        # First interpolate to the correct latitude, giving two values at
        # (lon[iW], node_lat) and (lon[iE], node_lat)
        shelfW = (shelf[jN,iW]-shelf[jS,iW])/(lat[jN]-lat[jS])*(node_lat[n]-lat[jS]) + shelf[jS,iW]
        shelfE = (shelf[jN,iE]-shelf[jS,iE])/(lat[jN]-lat[jS])*(node_lat[n]-lat[jS]) + shelf[jS,iE]
        # Now interpolate to the correct longitude
        node_shelf[n] = (shelfE-shelfW)/(lon[iE]-lon[iW])*(node_lon[n]-lon[iW]) + shelfW
        # Similarly for imask
        imaskW = (imask[jN,iW]-imask[jS,iW])/(lat[jN]-lat[jS])*(node_lat[n]-lat[jS]) + imask[jS,iW]
        imaskE = (imask[jN,iE]-imask[jS,iE])/(lat[jN]-lat[jS])*(node_lat[n]-lat[jS]) + imask[jS,iE]
        node_imask[n] = (imaskE-imaskW)/(lon[iE]-lon[iW])*(node_lon[n]-lon[iW]) + imaskW

        # Check how many nonzero shelf points were in the RTopo window
        shelf_cut = shelf[jS:jN+1,iW:iE+1]
        num_is = count_nonzero(shelf_cut<shelf_min)
        if num_is == 0:
            # No such points
            # Look at a bigger window
            shelf_cut = shelf[jS-1:jN+2,iW-1:iE+2]
            num_is = count_nonzero(shelf_cut<shelf_min)
            if num_is == 0:
                # Leave it; this is actually a grounded ice point
                pass
            else:
                # Replace shelf with the average of these ice shelf points
                print 'Node ' + str(n+1) + ' is isolated, fixing'
                index = nonzero(shelf_cut<shelf_min)
                node_shelf[n] = mean(shelf_cut[index])
        elif num_is < 4:
            # Less than 4 such points
            print 'Node ' + str(n+1) + ' is nearly isolated, fixing'
            # Replace shelf with the average of these ice shelf points
            index = nonzero(shelf_cut<shelf_min)
            node_shelf[n] = mean(shelf_cut[index])

    print 'Removing artifacts in interior'
    # Find nodes with imask<0.5 which are not coastal nodes
    index = nonzero((node_imask<0.5)*(node_coast==0))
    # Set shelf to zero here
    node_shelf[index] = 0.0

    # Fill holes in interpolated ice shelf
    print 'Filling holes in ice shelf'
    for n in range(num_nod2d):
        # Find coastal nodes with zero ice shelf
        if node_coast[n] == 1 and node_shelf[n] == 0.0:
            # Find elements which this node is part of
            target_elem = nonzero((elem2d[:,0]==n+1)+(elem2d[:,1]==n+1)+(elem2d[:,2]==n+1))[0]
            # Figure out how many nodes in these elements have nonzero shelf,
            # and integrate those shelf values. The nodes will all be double
            # counted but that's ok.
            num_is = 0
            shelf_sum = 0.0
            for e in target_elem:
                target_nodes = elem2d[e,:]
                for m in target_nodes:
                    if node_shelf[m-1] > 0:
                        num_is += 1
                        shelf_sum += node_shelf[m-1]
            # If there are more than 2 such nodes (actually 1 because they are
            # double-counted), set the current shelf value to the average.
            if num_is > 2:
                print 'Modifying node ' + str(n+1)
                node_shelf[n] = shelf_sum/num_is

    # Write to file
    print 'Writing ' + output_file
    file = open(output_file, 'w')
    for val in node_shelf:
        file.write(str(round(val)) + '\n')
    file.close()


if __name__ == "__main__":

    shelf()

                
        

        
        
                    
                    
                    
                    
                    
                    
        

    
    
    

    
    

    
