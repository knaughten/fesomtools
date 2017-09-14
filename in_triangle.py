# Determine if the point (lon0, lat0) lies within the triangular Element elm.
# Returns a boolean.
def in_triangle (elm, lon0, lat0):

   alpha = ((elm.lat[1] - elm.lat[2])*(lon0 - elm.lon[2]) + (elm.lon[2] - elm.lon[1])*(lat0 - elm.lat[2]))/((elm.lat[1] - elm.lat[2])*(elm.lon[0] - elm.lon[2]) + (elm.lon[2] - elm.lon[1])*(elm.lat[0] - elm.lat[2]))
   beta = ((elm.lat[2] - elm.lat[0])*(lon0 - elm.lon[2]) + (elm.lon[0] - elm.lon[2])*(lat0 - elm.lat[2]))/((elm.lat[1] - elm.lat[2])*(elm.lon[0] - elm.lon[2]) + (elm.lon[2] - elm.lon[1])*(elm.lat[0] - elm.lat[2]))
   gamma = 1 - alpha - beta

   return alpha >= 0 and beta >= 0 and gamma >= 0
