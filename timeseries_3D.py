from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from fesom_grid import *
from unesco import *

# Calculate and plot timeseries of ocean heat content, average salinity, and 
# total kinetic energy (all restricted to the Southern Ocean i.e. south of 30S)
# during a FESOM simulation.
# Takes 32 GB memory on raijin for Kaitlin's low_res mesh, and about 20 minutes
# walltime per year of output (in 5-day averages)
# Input:
# ocn_file = path to output oce.mean.nc, assumed to have 5-day averages
# log_file = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
def timeseries_3D (mesh_path, ocn_file, log_file):

    circumpolar = True   # Only consider elements south of 30S
    cross_180 = False    # Don't make second copies of elements that cross 180E
    days_per_output = 5  # Number of days for each output step
    rhoCp = 4.2e6            # Volumetric heat capacity of seawater (J/K/m^3)
    C2K = 273.15         # Celsius to Kelvin conversion

    ohc = []
    avgsalt = []
    tke = []
    # Check if the log file exists
    if exists(log_file):
        print('Reading previously calculated values')
        f = open(log_file, 'r')
        # Skip the first line (header)
        f.readline()
        for line in f:
            try:
                ohc.append(float(line))
            except(ValueError):
                # Reached the header for the next variable
                break
        for line in f:
            try:
                avgsalt.append(float(line))
            except(ValueError):
                break
        for line in f:
            tke.append(float(line))
        f.close()

    print('Building grid')
    elements = fesom_grid(mesh_path, circumpolar, cross_180)
    # Also read the depth of each node
    f = open(mesh_path + 'nod3d.out', 'r')
    f.readline()
    depth = []
    for line in f:
        tmp = line.split()
        depth.append(float(tmp[3]))
    f.close()
    # Convert to pressure in bar
    press = abs(array(depth))/10.0

    print('Reading data')
    id = Dataset(ocn_file, 'r')
    num_time = id.variables['time'].shape[0]
    temp = id.variables['temp'][:,:]
    salt = id.variables['salt'][:,:]
    u = id.variables['u'][:,:]
    v = id.variables['v'][:,:]
    id.close()

    print('Calculating density')
    rho = unesco(temp, salt, tile(press, (num_time,1)))

    print('Setting up arrays')
    # First calculate volume of each element
    dV_e3d = []
    # Loop over 2D elements
    for elm in elements:
        # Select the three nodes making up this element
        nodes = [elm.nodes[0], elm.nodes[1], elm.nodes[2]]
        # Calculate area of the surface triangle
        area = elm.area()
        # Loop downward through the water column
        while True:
            if nodes[0].below is None or nodes[1].below is None or nodes[2].below is None:
                # We've reached the bottom
                break
            # Calculate volume as area * average depth
            dV_e3d.append(area*(abs(nodes[0].depth - nodes[0].below.depth) + abs(nodes[1].depth - nodes[1].below.depth) + abs(nodes[2].depth - nodes[2].below.depth))/3.0)
            # Update nodes
            for i in range(3):
                nodes[i] = nodes[i].below
    dV_e3d = array(dV_e3d)

    # Set up arrays for timeseries of variables at each 3D element
    temp_e3d = zeros([num_time,size(dV_e3d)])
    salt_e3d = zeros([num_time,size(dV_e3d)])
    rho_e3d = zeros([num_time,size(dV_e3d)])
    u_e3d = zeros([num_time,size(dV_e3d)])
    v_e3d = zeros([num_time,size(dV_e3d)])
    # Loop over 2D elements again
    j = 0
    for elm in elements:
        # Select the three nodes making up this element
        nodes = [elm.nodes[0], elm.nodes[1], elm.nodes[2]]
        # Loop downward through the water column
        while True:
            if nodes[0].below is None or nodes[1].below is None or nodes[2].below is None:
                # We've reached the bottom
                break
            # Value of each variable in this triangular prism is the
            # average of the six vertices
            temp_e3d[:,j] = (temp[:,nodes[0].id] + temp[:,nodes[1].id] + temp[:,nodes[2].id] + temp[:,nodes[0].below.id] + temp[:,nodes[1].below.id] + temp[:,nodes[2].below.id])/6.0
            salt_e3d[:,j] = (salt[:,nodes[0].id] + salt[:,nodes[1].id] + salt[:,nodes[2].id] + salt[:,nodes[0].below.id] + salt[:,nodes[1].below.id] + salt[:,nodes[2].below.id])/6.0
            rho_e3d[:,j] = (rho[:,nodes[0].id] + rho[:,nodes[1].id] + rho[:,nodes[2].id] + rho[:,nodes[0].below.id] + rho[:,nodes[1].below.id] + rho[:,nodes[2].below.id])/6.0
            u_e3d[:,j] = (u[:,nodes[0].id] + u[:,nodes[1].id] + u[:,nodes[2].id] + u[:,nodes[0].below.id] + u[:,nodes[1].below.id] + u[:,nodes[2].below.id])/6.0
            v_e3d[:,j] = (v[:,nodes[0].id] + v[:,nodes[1].id] + v[:,nodes[2].id] + v[:,nodes[0].below.id] + v[:,nodes[1].below.id] + v[:,nodes[2].below.id])/6.0
            # Update nodes
            for i in range(3):
                nodes[i] = nodes[i].below
            j += 1
    
    print('Building timeseries')
    for t in range(num_time):
        # Integrate temp*rhoCp*dV to get OHC
        ohc.append(sum((temp_e3d[t,:]+C2K)*rhoCp*dV_e3d))
        # Average salinity (weighted with rho*dV)
        avgsalt.append(sum(salt_e3d[t,:]*rho_e3d[t,:]*dV_e3d)/sum(rho_e3d[t,:]*dV_e3d))
        # Integrate 0.5*rho*speed^2*dV to get TKE
        tke.append(sum(0.5*rho_e3d[t,:]*(u_e3d[t,:]**2 + v_e3d[t,:]**2)*dV_e3d))

    # Calculate time values
    time = arange(len(ohc))*days_per_output/365.

    print('Plotting ocean heat content')
    clf()
    plot(time, ohc)
    xlabel('Years')
    ylabel('Southern Ocean Heat Content (J)')
    grid(True)
    savefig('ohc.png')

    print('Plotting average salinity')
    clf()
    plot(time, avgsalt)
    xlabel('Years')
    ylabel('Southern Ocean Average Salinity (psu)')
    grid(True)
    savefig('avgsalt.png')

    print('Plotting total kinetic energy')
    clf()
    plot(time, tke)
    xlabel('Years')
    ylabel('Southern Ocean Total Kinetic Energy (J)')
    grid(True)
    savefig('tke.png')

    print('Saving results to log file')
    f = open(log_file, 'w')
    f.write('Southern Ocean Heat Content (J):\n')
    for elm in ohc:
        f.write(str(elm) + '\n')
    f.write('Southern Ocean Average Salinity (psu):\n')
    for elm in avgsalt:
        f.write(str(elm) + '\n')
    f.write('Southern Ocean Total Kinetic Energy (J):\n')
    for elm in tke:
        f.write(str(elm) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    mesh_path = input("Path to FESOM mesh directory: ")
    ocn_file = input("Path to FESOM oce.mean.nc output file: ")
    log_file = input("Path to logfile to save values and/or read previously calculated values: ")
    timeseries_3D(mesh_path, ocn_file, log_file)

    
        

    
