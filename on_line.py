# Determine if the point (lonA, latA) falls on the line through the points
# (lon0, lat0) and (lon1, lat1). Returns a boolean.
def on_line (lonA, latA, lon0, lat0, lon1, lat1):

    return latA == (lat1-lat0)/(lon1-lon0)*(lonA-lon0) + lat0
