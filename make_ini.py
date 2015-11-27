from netCDF4 import Dataset
from numpy import *

# Make a FESOM initial conditions file using ECCO2 temperature and salinity for
# January 1995.

# Paths to ECCO2 temperature and salinity files
theta_file = 'THETA.1440x720x50.199501.nc'
salt_file = 'SALT.1440x720x50.199501.nc'
# Path to desired initial conditions file
out_file = '../ecco2_199501_ts.out'

# Read ECCO2 grid and data
theta_fid = Dataset(theta_file, 'r')
lon_ecco = theta_fid.variables['LONGITUDE_T'][:]
lat_ecco = theta_fid.variables['LATITUDE_T'][:]
depth_ecco = -theta_fid.variables['DEPTH_T'][:]
theta = theta_fid.variables['THETA'][0,:,:,:]
theta_fid.close()
salt_fid = Dataset(salt_file, 'r')
salt = salt_fid.variables['SALT'][0,:,:,:]
salt_fid.close()

# Fill land masks with constant values 
theta[theta.mask] = -0.5
salt[salt.mask] = 34.5

f=open(out_file, 'w')

# Header
sizes = array([size(lon_ecco), size(lat_ecco), size(depth_ecco)])

print 'Writing grid'
sizes.tofile(f, ' ')
f.write('\n')

# Columns of 5
posn = 0
while True:
    end = min(posn+5, size(lon_ecco))
    data = lon_ecco[posn:end]
    data.tofile(f, ' ')
    f.write('\n')
    if end == size(lon_ecco):
        break
    else:
        posn = posn+5

posn = 0
while True:
    end = min(posn+5, size(lat_ecco))
    data = lat_ecco[posn:end]
    data.tofile(f, ' ')
    f.write('\n')
    if end == size(lat_ecco):
        break
    else:
        posn = posn+5

posn = 0
while True:
    end = min(posn+5, size(depth_ecco))
    data = depth_ecco[posn:end]
    data.tofile(f, ' ')
    f.write('\n')
    if end == size(depth_ecco):
        break
    else:
        posn = posn+5

print 'Setting up temperature'
# Flatten the 3D array
theta_flat = []
for i in range(size(lon_ecco)):
    for j in range(size(lat_ecco)):
        for k in range(size(depth_ecco)):
            theta_flat.append(theta[k,j,i])

print 'Writing temperature'
theta_flat = array(theta_flat)
posn = 0
while True:
    end = min(posn+5, size(theta_flat))
    data = theta_flat[posn:end]
    data.tofile(f, ' ')
    f.write('\n')
    if end == size(theta_flat):
        break
    else:
        posn = posn+5

print 'Setting up salinity'
salt_flat = []
for i in range(size(lon_ecco)):
    for j in range(size(lat_ecco)):
        for k in range(size(depth_ecco)):
            salt_flat.append(salt[k,j,i])

print 'Writing salinity'
salt_flat = array(salt_flat)
posn = 0
while True:
    end = min(posn+5, size(salt_flat))
    data = salt_flat[posn:end]
    data.tofile(f, ' ')
    f.write('\n')
    if end == size(salt_flat):
        break
    else:
        posn = posn+5
f.close()
