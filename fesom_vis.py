from patches import *
from lonlat_plot import *
from zonal_slice_plot import *
from zonal_avg_plot import *

# Command-line interface for FESOM plots

def fesom_vis_lonlat ():

    # Get path to mesh directory
    mesh_path = raw_input("Path to mesh directory: ")    

    # Get file name and variable name
    file_path = raw_input("Path to FESOM output file: ")
    var_name = raw_input("Variable name: ")
    if var_name == 'wnet' or '.ice.' in file_path:
        # Mask ice shelf cavities for all sea ice variables; mask open ocean for
        # ice shelf melt rate
        mask_cavities = True
    else:
        mask_cavities = False

    # Get depth information
    depth_type = raw_input("Single depth (s) or vertical average (v)? ")
    if depth_type == 's':
        depth_input = raw_input("Surface nodes (s), bottom nodes (b), or specific depth (d)? ")
        if depth_input == 's':
            depth_key = 0
            depth = NaN
            depth_bounds = None
        elif depth_input == 'b':
            depth_key = 1
            depth = NaN
            depth_bounds = None
        elif depth_input == 'd':
            depth_key = 3
            depth = float(raw_input("Enter depth (positive, in metres): "))
            depth_bounds = None
    elif depth_type == 'v':
        depth_input = raw_input("Vertical average throughout the entire water column (w) or between two specific depths (d)? ")
        if depth_input == 'w':
            depth_key = 2
            depth = NaN
            depth_bounds = None
        elif depth_input == 'd':
            depth_key = 4
            depth = NaN
            shallow_bound = float(raw_input("Enter shallow depth bound (positive, in metres): "))
            deep_bound = float(raw_input("Enter deep depth bound (positive, in metres): "))
            depth_bounds = [shallow_bound, deep_bound]

    # Get index of time axis in FESOM output file
    tstep = int(raw_input("Timestep number: "))

    # Current domain options are global or circumpolar Antarctic (to 60S)
    domain = raw_input("Global (g) or circumpolar (c)? ")
    if domain == 'c':
        circumpolar = True
    elif domain == 'g':
        circumpolar = False
    else:
        print "Problem with global/circumpolar choice"
        exit()

    # Get colour bounds if necessary
    set_limits = False
    limits = None
    get_bounds = raw_input("Set bounds on colour scale (y/n)? ")
    if get_bounds == 'y':
        set_limits = True
        lower_bound = float(raw_input("Lower bound: "))
        upper_bound = float(raw_input("Upper bound: "))
        limits = [lower_bound, upper_bound]

    # Current options are to save figure as a file, or display in window
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    else:
        print "Problem with save/display choice"
        exit()

    # Build FESOM grid
    elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)

    # Call lonlat_plot
    lonlat_plot(mesh_path, file_path, var_name, depth_key, depth, depth_bounds, tstep, circumpolar, elements, patches, mask_cavities, save, fig_name, set_limits, limits)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            new_grid = False
            while True:
                # Ask for changes to the input parameters; repeat until the user
                # is finished
                changes = raw_input("Enter a parameter to change: (1) mesh path, (2) file path, (3) variable name, (4) depth, (5) timestep number, (6) global/circumpolar, (7) colour bounds, (8) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New mesh
                        mesh_path = raw_input("Path to mesh directory: ")
                        # We will have to make a new grid
                        new_grid = True
                        file_path = raw_input("Path to FESOM output file: ")
                    elif int(changes) == 2:
                        # New FESOM output file
                        file_path = raw_input("Path to FESOM output file: ")
                    elif int(changes) == 3:
                        # New variable name
                        var_name = raw_input("Variable name: ")
                    elif int(changes) == 4:
                        # New depth information
                        depth_type = raw_input("Single depth (s) or vertical average (v)? ")
                        if depth_type == 's':
                            depth_input = raw_input("Surface nodes (s), bottom nodes (b), or specific depth (d)? ")
                            if depth_input == 's':
                                depth_key = 0
                                depth = NaN
                                depth_bounds = None
                            elif depth_input == 'b':
                                depth_key = 1
                                depth = NaN
                                depth_bounds = None
                            elif depth_input == 'd':
                                depth_key = 3
                                depth = float(raw_input("Enter depth (positive, in metres): "))
                                depth_bounds = None
                        elif depth_type == 'v':
                            depth_input = raw_input("Vertical average throughout the entire water column (w) or between two specific depths (d)? ")
                            if depth_input == 'w':
                                depth_key = 2
                                depth = NaN
                                depth_bounds = None
                            elif depth_input == 'd':
                                depth_key = 4
                                depth = NaN
                                shallow_bound = float(raw_input("Enter shallow depth bound (positive, in metres): "))
                                deep_bound = float(raw_input("Enter deep depth bound (positive, in metres): "))
                                depth_bounds = [shallow_bound, deep_bound]
                    elif int(changes) == 5:
                        # New time index
                        tstep = int(raw_input("Timestep number: "))
                    elif int(changes) == 6:
                        # Change from global to circumpolar, or vice versa
                        circumpolar = not circumpolar
                        # We will have to make a new grid
                        new_grid = True
                    elif int(changes) == 7:
                        # New colour bounds
                        set_limits = False
                        limits = None
                        get_bounds = raw_input("Set bounds on colour scale (y/n)? ")
                        if get_bounds == 'y':
                            set_limits = True
                            lower_bound = float(raw_input("Lower bound: "))
                            upper_bound = float(raw_input("Upper bound: "))
                            limits = [lower_bound, upper_bound]
                    elif int(changes) == 8:
                        # Change from display to save, or vice versa
                        save = not save
                    else:
                        print "Invalid option"

                    if var_name == 'wnet' or '.ice.' in file_path:
                        # mask_cavities will be true
                        tmp_mask = True
                    else:
                        tmp_mask = False
                    # Check if value of mask_cavities will change since last time;
                    # if so, we need to make a new grid
                    if tmp_mask != mask_cavities:
                        new_grid = True
                    mask_cavities = tmp_mask

            if save:
                # Get file name for figure
                fig_name = raw_input("File name for figure: ")

            if new_grid:
                # Build a new grid if necessary
                elements, patches = make_patches(mesh_path, circumpolar, mask_cavities)
                
            # Call lonlat_plot
            lonlat_plot(mesh_path, file_path, var_name, depth_key, depth, depth_bounds, tstep, circumpolar, elements, patches, mask_cavities, save, fig_name, set_limits, limits)

        else:
            break


def fesom_vis_latdepth ():

    # Get path to mesh directory
    mesh_path = raw_input("Path to mesh directory: ")

    # Get file name and variable name
    file_path = raw_input("Path to FESOM output file: ")
    var_name = raw_input("Variable name: ")

    # Get index of time axis in FESOM output file
    tstep = int(raw_input("Timestep number: "))

    action = raw_input("Zonal slice (s) or zonal average (a)? ")
    if action == 's':
        avg = False
        # Get longitude for the zonal slice
        lon0 = float(raw_input("Longitude in degrees (positive east, negative west): "))
    elif action == 'a':
        avg = True
        # Get longitude bounds for the zonal average
        lon_min = float(raw_input("Minimum longitude (positive east, negative west): "))
        lon_max = float(raw_input("Maximum longitude (positive east, negative west): "))
    else:
        print 'Problem with zonal slice/average choice'
        exit()
    

    # Get cutoff depth and convert to negative
    depth_min = -1*float(raw_input("Deepest depth to plot (positive, metres): "))

    # Get colour bounds if necessary
    set_limits = False
    limits = None
    get_bounds = raw_input("Set bounds on colour scale (y/n)? ")
    if get_bounds == 'y':
        set_limits = True
        lower_bound = float(raw_input("Lower bound: "))
        upper_bound = float(raw_input("Upper bound: "))
        limits = [lower_bound, upper_bound]

    # Current options are to save figure as a file, or display in window
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    else:
        print "Problem with save/display choice"
        exit()

    # Make the figure
    if avg:
        zonal_avg_plot(mesh_path, file_path, var_name, tstep, lon_min, lon_max, depth_min, save, fig_name, set_limits, limits)
    else:
        zonal_slice_plot(mesh_path, file_path, var_name, tstep, lon0, depth_min, save, fig_name, set_limits, limits)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask to changes to the input parameters; repeat until the user
                # is finished
                changes = raw_input("Enter a parameter to change: (1) mesh path, (2) file path, (3) variable name, (4) timestep, (5) longitude, (6) cutoff depth, (7) colour bounds, (8) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New mesh
                        mesh_path = raw_input("Path to mesh directory: ")
                    elif int(changes) == 2:
                        # New FESOM output file
                        file_path = raw_input("Path to FESOM output file: ")
                    elif int(changes) == 3:
                        # New variable name
                        var_name = raw_input("Variable name: ")
                    elif int(changes) == 4:
                        # New time index 
                        tstep = int(raw_input("Timestep number: "))
                    elif int(changes) == 5:
                        action = raw_input("Zonal slice (s) or zonal average (a)? ")
                        if action == 's':
                            avg = False
                            # Get longitude for the zonal slice
                            lon0 = float(raw_input("Longitude in degrees (positive east, negative west): "))
                        elif action == 'a':
                            avg = True
                            # Get longitude bounds for the zonal average
                            lon_min = float(raw_input("Minimum longitude (positive east, negative west): "))
                            lon_max = float(raw_input("Maximum longitude (positive east, negative west): "))
                        else:
                            print 'Problem with zonal slice/average choice'
                            exit()
                    elif int(changes) == 6:
                        # New cutoff depth
                        depth_min = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
                    elif int(changes) == 7:
                        # New colour bounds
                        set_limits = False
                        limits = None
                        get_bounds = raw_input("Set bounds on colour scale (y/n)? ")
                        if get_bounds == 'y':
                            set_limits = True
                            lower_bound = float(raw_input("Lower bound: "))
                            upper_bound = float(raw_input("Upper bound: "))
                            limits = [lower_bound, upper_bound]
                    elif int(changes) == 8:
                        # Change from display to save, or vice versa
                        save = not save
                    else:
                        print "Invalid option"
            if save:
                # Get file name for figure
                fig_name = raw_input("File name for figure: ")

            # Make the figure
            if avg:
                zonal_avg_plot(mesh_path, file_path, var_name, tstep, lon_min, lon_max, depth_min, save, fig_name, set_limits, limits)
            else:
                zonal_slice_plot(mesh_path, file_path, var_name, tstep, lon0, depth_min, save, fig_name, set_limits, limits)
        else:
            break


if __name__ == "__main__":

    # Choose orientation of plot
    orientation = int(raw_input("Lon-lat plot (1) or lat-depth plot (2)? "))
    if orientation == 1:
        fesom_vis_lonlat()
    elif orientation == 2:
        fesom_vis_latdepth()
    else:
        print "Problem with horizontal/vertical choice"
        exit()

    
            

                
                
                        



    

