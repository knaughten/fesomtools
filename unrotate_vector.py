from numpy import *

# Given a 2D vector on the rotated FESOM grid, unrotate it.
# Input:
# rlon, rlat = rotated longitude and latitude in degrees
# ur, vr = rotated vector components
# Output:
# ug, vg = unrotated vector components
def unrotate_vector (rlon, rlat, ur, vr):

    # Rotation parameters (check your FESOM namelists)
    alpha = 50
    beta = 15
    gamma = -90
    # Degrees to radians conversion factor)
    deg2rad = pi/180

    # Convert everything to radians
    alpha *= deg2rad
    beta *= deg2rad
    gamma *= deg2rad
    rlon *= deg2rad
    rlat *= deg2rad

    rotate_matrix = zeros([3, 3])
    rotate_matrix[0,0] = cos(gamma)*cos(alpha) - sin(gamma)*cos(beta)*sin(alpha)
    rotate_matrix[0,1] = cos(gamma)*sin(alpha) + sin(gamma)*cos(beta)*cos(alpha)
    rotate_matrix[0,2] = sin(gamma)*sin(beta)
    rotate_matrix[1,0] = -sin(gamma)*cos(alpha) - cos(gamma)*cos(beta)*sin(alpha)
    rotate_matrix[1,1] = -sin(gamma)*sin(alpha) + cos(gamma)*cos(beta)*cos(alpha)
    rotate_matrix[1,2] = cos(gamma)*sin(beta)
    rotate_matrix[2,0] = sin(beta)*sin(alpha)
    rotate_matrix[2,1] = -sin(beta)*cos(alpha)
    rotate_matrix[2,2] = cos(beta)

    # Convert to rotated Cartesian coordinates
    xr = cos(rlat)*cos(rlon)
    yr = cos(rlat)*sin(rlon)
    zr = sin(rlat)

    # Convert to geographical Cartesian coordinates
    xg = rotate_matrix[0,0]*xr + rotate_matrix[1,0]*yr + rotate_matrix[2,0]*zr
    yg = rotate_matrix[0,1]*xr + rotate_matrix[1,1]*yr + rotate_matrix[2,1]*zr
    zg = rotate_matrix[0,2]*xr + rotate_matrix[1,2]*yr + rotate_matrix[2,2]*zr

    # Convert to geographical lon-lat coordinates
    glat = arcsin(zg)
    glon = arctan2(yg,xg)

    # Convert vector components to rotated Cartesian space
    vector_xr = -vr*sin(rlat)*cos(rlon) - ur*sin(rlon)
    vector_yr = -vr*sin(rlat)*sin(rlon) + ur*cos(rlon)
    vector_zr = vr*cos(rlat)

    # Convert vector components to geographical Cartesian space
    vector_xg = rotate_matrix[0,0]*vector_xr + rotate_matrix[1,0]*vector_yr + rotate_matrix[2,0]*vector_zr
    vector_yg = rotate_matrix[0,1]*vector_xr + rotate_matrix[1,1]*vector_yr + rotate_matrix[2,1]*vector_zr
    vector_zg = rotate_matrix[0,2]*vector_xr + rotate_matrix[1,2]*vector_yr + rotate_matrix[2,2]*vector_zr

    # Convert vector components to geographical lat-lon space
    vg = -sin(glat)*cos(glon)*vector_xg - sin(glat)*sin(glon)*vector_yg + cos(glat)*vector_zg
    ug = -sin(glon)*vector_xg + cos(glon)*vector_yg

    return ug, vg
    

    

    

    

    

    
