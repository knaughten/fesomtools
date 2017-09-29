# Adapated from https://gamedev.stackexchange.com/questions/21096/what-is-an-efficient-2d-line-segment-versus-triangle-intersection-test

# Determine whether points A and B lie on the same side of the line through
# points C and D (in lon-lat space).
# Input: lon and lat coordinates of points A, B, C, D
# Output: float z such that z > 0 if A and B are on the same side of line CD,
#                           z < 0 if A and B are on opposite sides,
#                           z == 0 if A and/or B are on the line itself.
def side_of_line (lonA, latA, lonB, latB, lonC, latC, lonD, latD):

    z1 = (lonD - lonC)*(latA - latC) - (lonA - lonC)*(latD - latC)
    z2 = (lonD - lonC)*(latB - latC) - (lonB - lonC)*(latD - latC)
    return z1*z2
