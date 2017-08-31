from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from fesom_grid import *
from unrotate_vector import *
from triangle_area import *

def div_peninsula_res (mesh_path_lr, seasonal_file_lr, mesh_path_hr, seasonal_file_hr):

    # Bounds on regular grid
    lon_min = -80
    lon_max = -30
    lat_min = -80
    lat_max = -60
    # Size of regular grid
    num_lon = 500
    num_lat = 500
    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Bounds on plot (for polar coordinate transformation)
    x_min = -22.5
    x_max = -12
    y_min = 6
    y_max = 15
    # Season names for title
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    # Bound for colour scale (10^-7 m/s)
    bound = 2

    print 'Building FESOM mesh'
    elements_lr = fesom_grid(mesh_path_lr, circumpolar=False, cross_180=False)
    elements_hr = fesom_grid(mesh_path_hr, circumpolar=False, cross_180=False)
    # Read rotated latitude and longitude at each node
    # Low-res
    file = open(mesh_path_lr + 'nod2d.out', 'r')
    file.readline()
    rlon_lr = []
    rlat_lr = []
    for line in file:
        tmp = line.split()
        lon_tmp = float(tmp[1])
        lat_tmp = float(tmp[2])
        if lon_tmp < -180:
            lon_tmp += 360
        elif lon_tmp > 180:
            lon_tmp -= 360
        rlon_lr.append(lon_tmp)
        rlat_lr.append(lat_tmp)
    file.close()
    rlon_lr = array(rlon_lr)
    rlat_lr = array(rlat_lr)
    # High-res
    file = open(mesh_path_hr + 'nod2d.out', 'r')
    file.readline()
    rlon_hr = []
    rlat_hr = []
    for line in file:
        tmp = line.split()
        lon_tmp = float(tmp[1])
        lat_tmp = float(tmp[2])
        if lon_tmp < -180:
            lon_tmp += 360
        elif lon_tmp > 180:
            lon_tmp -= 360
        rlon_hr.append(lon_tmp)
        rlat_hr.append(lat_tmp)
    file.close()
    rlon_hr = array(rlon_hr)
    rlat_hr = array(rlat_hr)

    print 'Reading data'
    # Read uhice and vhice
    id = Dataset(seasonal_file_lr, 'r')
    urhice_nodes_lr = id.variables['uhice'][:,:]
    vrhice_nodes_lr = id.variables['vhice'][:,:]
    id.close()
    id = Dataset(seasonal_file_hr, 'r')
    urhice_nodes_hr = id.variables['uhice'][:,:]
    vrhice_nodes_hr = id.variables['vhice'][:,:]
    id.close()
    print 'Unrotating vectors'
    uhice_nodes_lr, vhice_nodes_lr = unrotate_vector(rlon_lr, rlat_lr, urhice_nodes_lr, vrhice_nodes_lr)
    uhice_nodes_hr, vhice_nodes_hr = unrotate_vector(rlon_hr, rlat_hr, urhice_nodes_hr, vrhice_nodes_hr)

    print 'Setting up regular grid'
    # Start with boundaries
    lon_reg_edges = linspace(lon_min, lon_max, num_lon+1)
    lat_reg_edges = linspace(lat_min, lat_max, num_lat+1)
    # Now get centres
    lon_reg = 0.5*(lon_reg_edges[:-1] + lon_reg_edges[1:])
    lat_reg = 0.5*(lat_reg_edges[:-1] + lat_reg_edges[1:])
    # Also get differentials in lon-lat space
    dlon = lon_reg_edges[1:] - lon_reg_edges[:-1]
    dlat = lat_reg_edges[1:] - lat_reg_edges[:-1]
    # Make 2D versions
    lon_reg_2d, lat_reg_2d = meshgrid(lon_reg, lat_reg)
    dlon_2d, dlat_2d = meshgrid(dlon, dlat)
    # Calculate polar coordinates transformation for plotting
    x_reg = -(lat_reg_2d+90)*cos(lon_reg_2d*deg2rad+pi/2)
    y_reg = (lat_reg_2d+90)*sin(lon_reg_2d*deg2rad+pi/2)
    # Calculate differentials in Cartesian space
    dx = r*cos(lat_reg_2d*deg2rad)*dlon_2d*deg2rad
    dy = r*dlat_2d*deg2rad
    # Set up arrays for uhice and vhice at each resolution
    uhice_reg_lr = zeros([4, num_lat, num_lon])
    vhice_reg_lr = zeros([4, num_lat, num_lon])
    uhice_reg_hr = zeros([4, num_lat, num_lon])
    vhice_reg_hr = zeros([4, num_lat, num_lon])

    print 'Interpolating to regular grid'
    # Low-res
    # For each element, check if a point on the regular lat-lon grid lies
    # within. If so, do barycentric interpolation to that point, at each
    # depth on the regular grid.
    for elm in elements_lr:
        # Don't care about ice shelf cavities
        if not elm.cavity:
            # Check if we are within domain of regular grid
            if amax(elm.lon) > lon_min and amin(elm.lon) < lon_max and amax(elm.lat) > lat_min and amin(elm.lat) < lat_max:
                # Find largest regular longitude value west of Element
                tmp = nonzero(lon_reg > amin(elm.lon))[0]
                if len(tmp) == 0:
                    # Element crosses the western boundary
                    iW = 0
                else:
                    iW = tmp[0] - 1
                # Find smallest regular longitude value east of Element
                tmp = nonzero(lon_reg > amax(elm.lon))[0]
                if len(tmp) == 0:
                    # Element crosses the eastern boundary
                    iE = num_lon
                else:
                    iE = tmp[0]
                # Find largest regular latitude value south of Element
                tmp = nonzero(lat_reg > amin(elm.lat))[0]
                if len(tmp) == 0:
                    # Element crosses the southern boundary
                    jS = 0
                else:
                    jS = tmp[0] - 1
                # Find smallest regular latitude value north of Element
                tmp = nonzero(lat_reg > amax(elm.lat))[0]
                if len(tmp) == 0:
                    # Element crosses the northern boundary
                    jN = num_lat
                else:
                    jN = tmp[0]
                for i in range(iW+1,iE):
                    for j in range(jS+1,jN):
                        # There is a chance that the regular gridpoint at (i,j)
                        # lies within this element
                        lon0 = lon_reg[i]
                        lat0 = lat_reg[j]
                        if in_triangle(elm, lon0, lat0):
                            # Yes it does
                            # Get area of entire triangle
                            area = triangle_area(elm.lon, elm.lat)
                            # Get area of each sub-triangle formed by
                            # (lon0, lat0)
                            area0 = triangle_area([lon0, elm.lon[1], elm.lon[2]], [lat0, elm.lat[1], elm.lat[2]])
                            area1 = triangle_area([lon0, elm.lon[0], elm.lon[2]], [lat0, elm.lat[0], elm.lat[2]])
                            area2 = triangle_area([lon0, elm.lon[0], elm.lon[1]], [lat0, elm.lat[0], elm.lat[1]])
                            # Find fractional area of each
                            cff = [area0/area, area1/area, area2/area]
                            # Loop over seasons
                            for season in range(4):
                                # Find value of uhice and vhice at each Node
                                uhice_vals = []
                                vhice_vals = []
                                for n in range(3):
                                    uhice_vals.append(uhice_nodes_lr[season,elm.nodes[n].id])
                                    vhice_vals.append(vhice_nodes_lr[season,elm.nodes[n].id])
                                # Barycentric interpolation to lon0, lat0
                                uhice_reg_lr[season,j,i] = sum(array(cff)*array(uhice_vals))
                                vhice_reg_lr[season,j,i] = sum(array(cff)*array(vhice_vals))
    # Save land mask: wherever identically zero
    mask_lr = ones([num_lat, num_lon])
    index = uhice_reg_lr[0,:,:] == 0.0
    mask_lr[index] = 0.0
    # Repeat for high-res
    for elm in elements_hr:
        if not elm.cavity:
            if amax(elm.lon) > lon_min and amin(elm.lon) < lon_max and amax(elm.lat) > lat_min and amin(elm.lat) < lat_max:
                tmp = nonzero(lon_reg > amin(elm.lon))[0]
                if len(tmp) == 0:
                    iW = 0
                else:
                    iW = tmp[0] - 1
                tmp = nonzero(lon_reg > amax(elm.lon))[0]
                if len(tmp) == 0:
                    iE = num_lon
                else:
                    iE = tmp[0]
                tmp = nonzero(lat_reg > amin(elm.lat))[0]
                if len(tmp) == 0:
                    jS = 0
                else:
                    jS = tmp[0] - 1
                tmp = nonzero(lat_reg > amax(elm.lat))[0]
                if len(tmp) == 0:
                    jN = num_lat
                else:
                    jN = tmp[0]
                for i in range(iW+1,iE):
                    for j in range(jS+1,jN):
                        lon0 = lon_reg[i]
                        lat0 = lat_reg[j]
                        if in_triangle(elm, lon0, lat0):
                            area = triangle_area(elm.lon, elm.lat)
                            area0 = triangle_area([lon0, elm.lon[1], elm.lon[2]], [lat0, elm.lat[1], elm.lat[2]])
                            area1 = triangle_area([lon0, elm.lon[0], elm.lon[2]], [lat0, elm.lat[0], elm.lat[2]])
                            area2 = triangle_area([lon0, elm.lon[0], elm.lon[1]], [lat0, elm.lat[0], elm.lat[1]])
                            cff = [area0/area, area1/area, area2/area]
                            for season in range(4):
                                uhice_vals = []
                                vhice_vals = []
                                for n in range(3):
                                    uhice_vals.append(uhice_nodes_hr[season,elm.nodes[n].id])
                                    vhice_vals.append(vhice_nodes_hr[season,elm.nodes[n].id])
                                uhice_reg_hr[season,j,i] = sum(array(cff)*array(uhice_vals))
                                vhice_reg_hr[season,j,i] = sum(array(cff)*array(vhice_vals))
    mask_hr = ones([num_lat, num_lon])
    index = uhice_reg_hr[0,:,:] == 0.0
    mask_hr[index] = 0.0

    print 'Calculating divergence'
    # Low-res
    div_lr = ma.empty(shape(uhice_reg_lr))
    # Loop over seasons
    for season in range(4):
        # First calculate the two derivatives
        duhice_dx_lr = zeros([num_lat, num_lon])
        # Forward difference approximation
        duhice_dx_lr[:,:-1] = (uhice_reg_lr[season,:,1:] - uhice_reg_lr[season,:,:-1])/dx[:,:-1]
        # Backward difference for the last row
        duhice_dx_lr[:,-1] = (uhice_reg_lr[season,:,-1] - uhice_reg_lr[season,:,-2])/dx[:,-1]
        dvhice_dy_lr = zeros([num_lat, num_lon])
        dvhice_dy_lr[:-1,:] = (vhice_reg_lr[season,1:,:] - vhice_reg_lr[season,:-1,:])/dy[:-1,:]
        dvhice_dy_lr[-1,:] = (vhice_reg_lr[season,-1,:] - vhice_reg_lr[season,-2,:])/dy[-1,:]
        # Sum for the divergence
        div_lr_tmp = duhice_dx_lr + dvhice_dy_lr
        # Apply land mask
        div_lr[season,:,:] = ma.masked_where(mask_lr==0, div_lr_tmp)
    # High-res
    div_hr = ma.empty(shape(uhice_reg_hr))
    for season in range(4):
        duhice_dx_hr = zeros([num_lat, num_lon])
        duhice_dx_hr[:,:-1] = (uhice_reg_hr[season,:,1:] - uhice_reg_hr[season,:,:-1])/dx[:,:-1]
        duhice_dx_hr[:,-1] = (uhice_reg_hr[season,:,-1] - uhice_reg_hr[season,:,-2])/dx[:,-1]
        dvhice_dy_hr = zeros([num_lat, num_lon])
        dvhice_dy_hr[:-1,:] = (vhice_reg_hr[season,1:,:] - vhice_reg_hr[season,:-1,:])/dy[:-1,:]
        dvhice_dy_hr[-1,:] = (vhice_reg_hr[season,-1,:] - vhice_reg_hr[season,-2,:])/dy[-1,:]
        div_hr_tmp = duhice_dx_hr + dvhice_dy_hr
        div_hr[season,:,:] = ma.masked_where(mask_hr==0, div_hr_tmp)

    # Multiply by 10^7 so colourbar is more readable
    div_lr *= 1e7
    div_hr *= 1e7
    # Set up colour levels
    lev = linspace(-bound, bound, 50)

    print 'Plotting'
    fig = figure(figsize=(19,9))
    for season in range(4):
        # Low-res
        ax = fig.add_subplot(2, 4, season+1, aspect='equal')
        contourf(x_reg, y_reg, div_lr[season,:,:], lev, cmap='RdBu_r', extend='both')
        title(season_names[season], fontsize=24)
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        axis('off')
        if season == 0:
            text(-24, 10, 'low-res', fontsize=24, ha='right')
        # High-res
        ax = fig.add_subplot(2, 4, season+5, aspect='equal')
        img = contourf(x_reg, y_reg, div_hr[season,:,:], lev, cmap='RdBu_r', extend='both')
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        axis('off')
        if season == 0:
            text(-24, 10, 'high-res', fontsize=24, ha='right')
    # Colourbar at bottom
    cbaxes = fig.add_axes([0.35, 0.04, 0.3, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, ticks=arange(-bound, bound+1, 1))
    cbar.ax.tick_params(labelsize=20)
    suptitle(r'FESOM sea ice flux divergence (10$^{-7}$ m/s), 1992-2016 average', fontsize=30)
    subplots_adjust(wspace=0.025,hspace=0.025)
    fig.show()
    fig.savefig('div_peninsula_res.png')
    


def in_triangle (elm, lon0, lat0):

   alpha = ((elm.lat[1] - elm.lat[2])*(lon0 - elm.lon[2]) + (elm.lon[2] - elm.lon[1])*(lat0 - elm.lat[2]))/((elm.lat[1] - elm.lat[2])*(elm.lon[0] - elm.lon[2]) + (elm.lon[2] - elm.lon[1])*(elm.lat[0] - elm.lat[2]))
   beta = ((elm.lat[2] - elm.lat[0])*(lon0 - elm.lon[2]) + (elm.lon[0] - elm.lon[2])*(lat0 - elm.lat[2]))/((elm.lat[1] - elm.lat[2])*(elm.lon[0] - elm.lon[2]) + (elm.lon[2] - elm.lon[1])*(elm.lat[0] - elm.lat[2]))
   gamma = 1 - alpha - beta

   return alpha >= 0 and beta >= 0 and gamma >= 0


# Command-line interface
if __name__ == "__main__":

    mesh_path_lr = raw_input("Path to low-res mesh directory: ")
    seasonal_file_lr = raw_input("Path to low-res seasonal climatology file containing thdgr: ")
    mesh_path_hr = raw_input("Path to high-res mesh directory: ")
    seasonal_file_hr = raw_input("Path to high-res seasonal climatology file containing thdgr: ")
    div_peninsula_res(mesh_path_lr, seasonal_file_lr, mesh_path_hr, seasonal_file_hr)
            

