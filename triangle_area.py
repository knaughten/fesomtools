from numpy import *

# Given the longitude and latitude coordinates of a triangular element,
# calculate the area of the triangle in m^2.
# High school trigonometry at its finest.
# Input:
# lon = array of length 3 containing longitude values (degrees east)
# lat = array of length 3 containing latitude values (degrees north)
# Output: area of triangle (real scalar) in m^2
def triangle_area (lon, lat):

    # Constants for spherical coordinates transformations
    r = 6.371e6 # Radius of the Earth in metres
    deg2rad = pi/180.0 # Conversion factor from degrees to radians

    # Convert from lon-lat space to x-y space
    # y = r*lat where lat is converted to radians (arc length formula)
    y1 = r*lat[0]*deg2rad
    y2 = r*lat[1]*deg2rad
    y3 = r*lat[2]*deg2rad
    # x = r*cos(lat)*lon where lat and lon are converted to radians
    # (arc length formula for a circle of radius r*cos(lat))
    x1 = r*cos(lat[0]*deg2rad)*lon[0]*deg2rad
    x2 = r*cos(lat[1]*deg2rad)*lon[1]*deg2rad
    x3 = r*cos(lat[2]*deg2rad)*lon[2]*deg2rad

    # Calculate the length of each side of the triangle (distance formula)
    l12 = sqrt((x2-x1)**2 + (y2-y1)**2)
    l13 = sqrt((x3-x1)**2 + (y3-y1)**2)
    l23 = sqrt((x3-x2)**2 + (y3-y2)**2)

    # Calculate two of the three angles in this triangle (cosine law)
    theta2 = arccos((l12**2 + l23**2 - l13**2)/(2*l12*l23))
    theta3 = arccos((l13**2 + l23**2 - l12**2)/(2*l23*l13))
    if theta2 < 1e-5 or theta3 < 1e-5 or isnan(theta2) or isnan(theta3):
        #print 'Warning: triangle area close to zero. Okay for barycentric interpolation. lon=' + str(lon) + ', lat=' + str(lat)
        return 0.0

    # Consider the triangle split down the middle by a line L1, containing
    # point 1 and intersecting the opposite side (l23) at a right angle.
    # Now we have two right triangles.
    # First find b2 and b3, the distances from points 2 and 3 respectively
    # to the intersection of the L1 and l23.
    b2 = l12*cos(theta2)    
    b3 = l13*cos(theta3)
    # Note that b2 + b3 should equal l23 within machine precision.
    # Make sure nothing went wrong here.
    error1 = abs((l23-b2-b3)/l23)
    if error1 > 0.01:
        print('Greater than 1% error in two values of base')
        exit
    # Similarly, find h2 and h3, which both measure the length of L1.
    h2 = l12*sin(theta2)
    h3 = l13*sin(theta3)
    # Note that h2 should equal h3 within machine precision.
    # Make sure nothing went wrong here.
    error2 = abs((h2-h3)/h2)    
    if error2 > 0.01:
        print('Greater than 1% error in two values of height')
        exit

    # The area of each right triangle is now 0.5*base*height.
    A2 = 0.5*b2*h2
    A3 = 0.5*b3*h3

    # The area of the original triangle is the sum of the two.
    return A2 + A3
    
    

    
