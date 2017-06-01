from numpy import *
from numpy.linalg import inv

# Unrotate longitude and latitude on the FESOM grid.
# Input:
# rlon, rlat = 1D arrays of rotated longitude and latitude in degrees
#              (i.e. straight from nod2d.out or nod3d.out), can be either
#              2D or 3D nodes
# Output:
# glon, glat = geographical longitude and latitude, in degrees
def unrotate_grid (rlon, rlat):

    # Grid rotation parameters (grep inside mesh_path if unsure; if they're not
    # mentioned, it's probably not a rotated grid, so set alpha=beta=gamma=0)
    alpha = 50
    beta = 15
    gamma = -90
    deg2rad = pi/180
    rad2deg = 180/pi

    # Unrotate grid
    alpha = alpha*deg2rad
    beta = beta*deg2rad
    gamma = gamma*deg2rad
    # Transformation matrix
    Tm = zeros((3,3))
    Tm[0,0] = cos(gamma)*cos(alpha) - sin(gamma)*cos(beta)*sin(alpha)
    Tm[0,1] = cos(gamma)*sin(alpha) + sin(gamma)*cos(beta)*cos(alpha)
    Tm[0,2] = sin(gamma)*sin(beta)
    Tm[1,0] = -sin(gamma)*cos(alpha) - cos(gamma)*cos(beta)*sin(alpha)
    Tm[1,1] = -sin(gamma)*sin(alpha) + cos(gamma)*cos(beta)*cos(alpha)
    Tm[1,2] = cos(gamma)*sin(beta)
    Tm[2,0] = sin(beta)*sin(alpha)
    Tm[2,1] = -sin(beta)*cos(alpha)
    Tm[2,2] = cos(beta)
    invTm = asarray(inv(matrix(Tm)))
    # Convert to radians
    rlon_rad = rlon*deg2rad
    rlat_rad = rlat*deg2rad
    # Rotated Cartesian coordinates
    x_rot = cos(rlat_rad)*cos(rlon_rad)
    y_rot = cos(rlat_rad)*sin(rlon_rad)
    z_rot = sin(rlat_rad)
    # Geographical Cartesian coordinates
    x_geo = invTm[0,0]*x_rot + invTm[0,1]*y_rot + invTm[0,2]*z_rot
    y_geo = invTm[1,0]*x_rot + invTm[1,1]*y_rot + invTm[1,2]*z_rot
    z_geo = invTm[2,0]*x_rot + invTm[2,1]*y_rot + invTm[2,2]*z_rot
    # Geographical lat-lon
    glat = arcsin(z_geo)
    glon = arctan2(y_geo, x_geo)
    index = nonzero(y_geo*x_geo == 0)
    glon[index] = 0
    # Convert back to degrees
    glon = glon*rad2deg
    glat = glat*rad2deg

    return glon, glat
