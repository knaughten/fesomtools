from patches import *
from global_plot import *
from circumpolar_plot import *

# Command-line interface for FESOM plots

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

# Current domain options are global or circumpolar Antarctic (to 60S)
domain = raw_input("Global (g) or circumpolar (c)? ")
if domain == 'c':
    circumpolar = True
elif domain == 'g':
    circumpolar = False
else:
    print "Problem with global/circumpolar choice"
    exit()

# Get file name and variable name
file_path = raw_input("Path to FESOM output file: ")
var_name = raw_input("Variable name: ")
if var_name == 'wnet' or '.ice.' in file_path:
    # Mask ice shelf cavities for all sea ice variables; mask open ocean for
    # ice shelf melt rate
    mask_cavities = True
else:
    mask_cavities = False

# Get index of time axis in FESOM output file
tstep = int(raw_input("Timestep number: "))

# Build FESOM grid
elements, patches = make_patches(circumpolar, mask_cavities)

# Call circumpolar_plot or global_plot depending on domain
if circumpolar:
    circumpolar_plot(file_path, var_name, tstep, elements, patches, mask_cavities, save, fig_name)
else:
    global_plot(file_path, var_name, tstep, elements, patches, mask_cavities, save, fig_name)

# Repeat until the user wants to exit
while True:
    repeat = raw_input("Make another plot (y/n)? ")
    if repeat == 'y':
        new_grid = False
        while True:
            # Ask for changes to the input parameters; repeat until the user
            # is finished
            changes = raw_input("Enter a parameter to change: (1) save/display, (2) global/circumpolar, (3) file path, (4) variable name, (5) timestep number; or enter to continue: ")
            if len(changes) == 0:
                # No more changes to parameters
                break
            else:
                if int(changes) == 1:
                    # Change from display to save, or vice versa
                    save = not save
                elif int(changes) == 2:
                    # Change from global to circumpolar, or vice versa
                    circumpolar = not circumpolar
                    # We will have to make a new grid
                    new_grid = True
                elif int(changes) == 3:
                    # New FESOM output file
                    file_path = raw_input("Path to FESOM output file: ")
                elif int(changes) == 4:
                    # New variable name
                    var_name = raw_input("Variable name: ")
                elif int(changes) == 5:
                    # New time index
                    tstep = int(raw_input("Timestep number: "))
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

        if new_grid:
            # Build a new grid if necessary
            elements, patches = make_patches(circumpolar, mask_cavities)
        if save:
            # Get file name for figure
            fig_name = raw_input("File name for figure: ")        

        # Call circumpolar_plot or global_plot depending on domain
        if circumpolar:
            circumpolar_plot(file_path, var_name, tstep, elements, patches, mask_cavities, save, fig_name)
        else:
            global_plot(file_path, var_name, tstep, elements, patches, mask_cavities, save, fig_name)
    else:
        break
            

                
                
                        



    

