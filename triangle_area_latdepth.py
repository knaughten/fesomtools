from numpy import *

# Given the latitude and depth coordinates of a triangular segement of a
# SideElement, calculate the area of the triangle in m^2.
# High school trigonometry at its finest.
# Input:
# lat = array of length 3 containing latitude values (degrees north)
# depth = array of length 3 containing depth values (metres, negative)
# Output: area of triangle (real scalar) in m^2
def triangle_area_latdepth (lat, depth):

    # Constants for spherical coordinates transformations
    r = 6.371e6 # Radius of the Earth in metres
    deg2rad = pi/180.0 # Conversion factor from degrees to radians

    # Convert from lon-lat space to x-y space
    # y = r*lat where lat is converted to radians (arc length formula)
    y1 = r*lat[0]*deg2rad
    y2 = r*lat[1]*deg2rad
    y3 = r*lat[2]*deg2rad
    # Nothing special about depth
    z1 = depth[0]
    z2 = depth[1]
    z3 = depth[2]

    # Calculate the length of each side of the triangle (distance formula)
    l12 = sqrt((y2-y1)**2 + (z2-z1)**2)
    l13 = sqrt((y3-y1)**2 + (z3-z1)**2)
    l23 = sqrt((y3-y2)**2 + (z3-z2)**2)

    # Calculate two of the three angles in this triangle (cosine law)
    theta2 = arccos((l12**2 + l23**2 - l13**2)/(2*l12*l23))
    theta3 = arccos((l13**2 + l23**2 - l12**2)/(2*l23*l13))

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
    
    

    
