from netCDF4 import Dataset
from numpy import *
from unrotate_grid import *
from matplotlib.pyplot import *

def coastal_current ():

    # File paths
    mesh_path = '/short/y99/kaa561/FESOM/mesh/meshB/'
    directory_beg = '/short/y99/kaa561/FESOM/highres_spinup/'
    directories = ['/short/y99/kaa561/FESOM/rcp45_M/', '/short/y99/kaa561/FESOM/rcp45_A/', '/short/y99/kaa561/FESOM/rcp85_M/', '/short/y99/kaa561/FESOM/rcp85_A/', '/short/y99/kaa561/FESOM/highres_spinup/']
    file_beg = 'annual_avg.oce.mean.1996.2005.nc'
    file_end = 'annual_avg.oce.mean.2091.2100.nc'
    # Titles for plotting
    expt_names = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A', 'CONTROL']
    num_expts = len(directories)
    # Colours for plotting
    beg_colour = (0,0,0)
    end_colours = [(0, 0.4, 1), (0.6, 0.2, 1), (0, 0.6, 0), (1, 0, 0.4), (0.6, 0.6, 0.6)]
    # Mesh parameters
    circumpolar = False
    cross_180 = False
    # Spacing of longitude bins
    dlon = 1
    # Parameters for continental shelf selection
    lat0 = -64  # Maximum latitude to consider
    h0 = 2500  # Deepest depth to consider

    # Set up longitude bins
    # Start with edges
    lon_bins = arange(-180, 180+dlon, dlon)
    # Centres for plotting
    lon_centres = 0.5*(lon_bins[:-1] + lon_bins[1:])
    num_bins = size(lon_centres)
    # Set up arrays to store maximum surface speed in each bin
    current_beg = zeros(num_bins)
    current_end = zeros([num_expts, num_bins])

    print 'Building mesh'
    # We only care about nodes, not elements, so don't need to use the
    # fesom_grid function.
    # Read cavity flag for each 2D surface node
    node_cavity = []
    f = open(mesh_path + 'cavity_flag_nod2d.out', 'r')
    for line in f:
        tmp = int(line)
        if tmp == 1:
            node_cavity.append(True)
        elif tmp == 0:
            node_cavity.append(False)
        else:
            print 'Problem with cavity flags'
    f.close()
    # Save the number of 2D nodes
    n2d = len(node_cavity)
    # Read rotated lat and lon for each node, also depth
    f = open(mesh_path + 'nod3d.out', 'r')
    f.readline()
    rlon = []
    rlat = []
    depth = []
    for line in f:
        tmp = line.split()
        lon_tmp = float(tmp[1])
        lat_tmp = float(tmp[2])
        depth_tmp = -1*float(tmp[3])
        if lon_tmp < -180:
            lon_tmp += 360
        elif lon_tmp > 180:
            lon_tmp -= 360
        rlon.append(lon_tmp)
        rlat.append(lat_tmp)
        depth.append(depth_tmp)
    f.close()
    # For lat and lon, only care about the 2D nodes (the first n2d indices)
    rlon = array(rlon[0:n2d])
    rlat = array(rlat[0:n2d])
    depth = array(depth)
    # Unrotate longitude
    lon, lat = unrotate_grid(rlon, rlat)
    # Read lists of which nodes are directly below which
    f = open(mesh_path + 'aux3d.out', 'r')
    max_num_layers = int(f.readline())
    node_columns = zeros([n2d, max_num_layers])
    for n in range(n2d):
        for k in range(max_num_layers):
            node_columns[n,k] = int(f.readline())
    node_columns = node_columns.astype(int)
    f.close()
    # Now figure out the bottom depth of each 2D node
    bottom_depth = zeros(n2d)
    for n in range(n2d):
        node_id = node_columns[n,0] - 1
        for k in range(1, max_num_layers):
            if node_columns[n,k] == -999:
                # Reached the bottom
                break
            node_id = node_columns[n,k] - 1
        # Save the last valid depth
        bottom_depth[n] = depth[node_id]
    # Delete this check later
    if amin(bottom_depth) < 0:
        # Probably accidentally kept a -999
        print 'Problem with bottom depths: minimum ' + amin(bottom_depth)      

    print 'Reading data'
    print '...1996-2005'
    id = Dataset(directory_beg + file_beg)
    # Only save surface nodes
    u_tmp = id.variables['u'][0,:n2d]
    v_tmp = id.variables['v'][0,:n2d]
    id.close()
    # Calculate speed
    speed_beg = sqrt(u_tmp**2 + v_tmp**2)
    # Set up array for speed in following experiments
    speed_end = zeros([num_expts, n2d])
    for expt in range(num_expts):
        print '...' + expt_names[expt]
        id = Dataset(directories[expt] + file_end)
        u_tmp = id.variables['u'][0,:n2d]
        v_tmp = id.variables['v'][0,:n2d]
        id.close()
        speed_end[expt,:] = sqrt(u_tmp**2 + v_tmp**2)

    print 'Selecting coastal current'
    for n in range(n2d):
        # Check if we care about this node
        if lat[n] <= lat0 and bottom_depth[n] <= h0 and not node_cavity[n]:
            # Find longitude bin
            lon_index = nonzero(lon_bins > lon[n])[0][0] - 1
            # Update coastal current speed in this bin if needed
            if speed_beg[n] > current_beg[lon_index]:
                current_beg[lon_index] = speed_beg[n]
            for expt in range(num_expts):
                if speed_end[expt,n] > current_end[expt,lon_index]:
                    current_end[expt,lon_index] = speed_end[expt,n]

    # Calculate percent change at each longitude bin for each experiment
    percent_change = zeros(shape(current_end))
    for expt in range(num_expts):
        percent_change[expt,:] = (current_end[expt,:] - current_beg)/current_beg*100

    print 'Plotting'
    # Coastal current
    fig = figure(figsize=(12,8))
    plot(lon_centres, current_beg, color=beg_colour, label='1996-2005')
    for expt in range(num_expts):
        plot(lon_centres, current_end[expt,:], color=end_colours[expt], label=expt_names[expt])
    grid(True)
    title('Coastal current speed', fontsize=20)
    xlabel('Longitude', fontsize=14)
    ylabel('m/s', fontsize=14)
    xlim([-180, 180])
    legend()
    fig.show()
    # Percent change
    fig = figure(figsize=(12,8))
    for expt in range(num_expts):
        plot(lon_centres, percent_change[expt,:], color=end_colours[expt], label=expt_names[expt])
    grid(True)
    title('Percent change in coastal current speed', fontsize=20)
    xlabel('Longitude', fontsize=14)
    ylabel('%', fontsize=14)
    xlim([-180, 180])
    legend()
    fig.show()

    print 'Mean coastal current, 1996-2005: ' + str(mean(current_beg)) + ' m/s'
    for expt in range(num_expts):
        print 'Mean coastal current, ' + expt_names[expt] + ': ' + str(mean(current_end[expt,:])) + ' m/s, mean percent change of ' + str(mean(percent_change[expt,:]))


# Command-line interface
if __name__ == "__main__":

    coastal_current()
    
        
    
    
