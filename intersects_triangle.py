from side_of_line import *

# Adapated from https://gamedev.stackexchange.com/questions/21096/what-is-an-efficient-2d-line-segment-versus-triangle-intersection-test

# Determine whether the line segment between (lon0, lat0) and (lon1, lat1)
# intersects the triangle defined by the Element elm. This includes overlapping
# (i.e. aligns with an edge of the triangle) and touching (i.e. one end of the
# segment touches an edge of the triangle). Returns a boolean.
def intersects_triangle (elm, lon0, lat0, lon1, lat1):

    f1 = side_of_line(lon0, lat0, elm.lon[2], elm.lat[2], elm.lon[0], elm.lat[0], elm.lon[1], elm.lat[1])
    f2 = side_of_line(lon1, lat1, elm.lon[2], elm.lat[2], elm.lon[0], elm.lat[0], elm.lon[1], elm.lat[1])
    f3 = side_of_line(lon0, lat0, elm.lon[0], elm.lat[0], elm.lon[1], elm.lat[1], elm.lon[2], elm.lat[2])
    f4 = side_of_line(lon1, lat1, elm.lon[0], elm.lat[0], elm.lon[1], elm.lat[1], elm.lon[2], elm.lat[2])
    f5 = side_of_line(lon0, lat0, elm.lon[1], elm.lat[1], elm.lon[2], elm.lat[2], elm.lon[0], elm.lat[0])
    f6 = side_of_line(lon1, lat1, elm.lon[1], elm.lat[1], elm.lon[2], elm.lat[2], elm.lon[0], elm.lat[0])
    f7 = side_of_line(elm.lon[0], elm.lat[0], elm.lon[1], elm.lat[1], lon0, lat0, lon1, lat1)
    f8 = side_of_line(elm.lon[1], elm.lat[1], elm.lon[2], elm.lat[2], lon0, lat0, lon1, lat1)

    if (f1 > 0 and f2 > 0 and f3 > 0 and f4 > 0 and f5 > 0 and f6 > 0):
        # Segment is completely inside triangle
        # This is a problem
        print('You need to choose a larger line segment.')
        exit
    elif (f1 < 0 and f2 < 0) or (f3 < 0 and f4 < 0) or (f5 < 0 and f6 < 0) or (f7 > 0 and f8 > 0):
        return False
    else:
        return True    
    

    
