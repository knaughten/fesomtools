from temp_salt_seasonal import *

# Call temp_salt_seasonal for a predetermined set of longitude slices and depth
# bounds. Save all the figures in a common directory.
# Input:
# file_path = path to FESOM oce.mean.nc climatology file (created using
#             average_years.py) with 5-day averages
# figure_dir = path to directory to store figures
# res_flag = 1 (for low-res) or 2 (for high-res)
def make_zonal_slices (file_path, figure_dir, res_flag):

    # Longitudes to plot
    lon = [-55, -40, -18, 0, 30, 71, 85, 97, 117, 145, 170, -160, -148, -120, -113, -105, -101, -95, -73]
    # Corresponding deepest depths to plot
    depth = [1700, 1500, 500, 1200, 500, 1300, 500, 500, 800, 500, 1000, 800, 400, 800, 800, 1100, 1000, 600, 800]           

    # Build FESOM grid
    if res_flag == 1:
        mesh_path = '/short/y99/kaa561/FESOM/mesh/low_res/'
    elif res_flag == 2:
        mesh_path = '/short/y99/kaa561/FESOM/mesh/high_res/'
    elements = fesom_grid(mesh_path)

    # Loop over longitudes
    for i in range(len(lon)):
        # Figure out how to save longitude in the filename
        if lon[i] < 0:
            fig_name = figure_dir + str(-lon[i]) + 'W.png'
        else:
            fig_name = figure_dir + str(lon[i]) + 'E.png'
        # Make the figure
        temp_salt_seasonal(elements, file_path, file_path, lon[i], -1*depth[i], True, fig_name)


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to output avg.oce.mean file: ")
    figure_dir = raw_input("Directory to store figures: ")
    res_flag = int(raw_input("Low resolution (1) or high resolution (2)? "))
    make_zonal_slices(file_path, figure_dir, res_flag)

    
