from numpy import *
from matplotlib.pyplot import *

def transport_3pt ():

    # Paths to RCP experiment directories
    directory_head = '/short/y99/kaa561/FESOM/'
    rcp_expt = ['rcp45_M_highres/', 'rcp45_A_highres/', 'rcp85_M_highres/', 'rcp85_A_highres/']
    num_rcps = len(rcp_expt)
    # Titles for plot
    rcp_titles = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A']
    # Colours for plotting
    rcp_colours = ['blue', 'cyan', 'green', 'magenta']
    # Paths to control experiment directories
    control_expt = 'highres_spinup/'
    # Title for plot
    control_title = 'CONTROL'
    # Colour for plotting
    control_colour = 'black'
    # Years to plot
    year_start = 1992
    rcp_year_start = 2006
    year_end = 2100

    # Subpolar gyre parameters
    subpolar_log = 'subpolar_gyres.log'
    # Drake Passage transport parameters
    dpt_log = 'dpt.log'
    # Spinup years to discard (first 2 repetitions of 1992-2005)
    control_skipyears = 28
    # Output steps per year
    peryear = 365/5

    # Build time axis: 1 point per year
    num_years = year_end - year_start + 1
    num_years_rep3 = rcp_year_start - year_start
    time = arange(year_start, year_end+1)

    # Drake Passage transport, annual averages
    dpt = empty([num_rcps+1, num_years])
    # Loop over RCP experiments
    for expt in range(num_rcps):
        # Read logfile
        dpt_tmp = []
        f = open(directory_head + rcp_expt[expt] + dpt_log)
        f.readline()
        for line in f:
            dpt_tmp.append(float(line))
        f.close()
        # Calculate annual averages
        for year in range(num_years):
            dpt[expt,year] = mean(array(dpt_tmp[peryear*year:peryear*(year+1)]))
    # Read control experiment
    # Read logfile
    dpt_tmp = []
    f = open(directory_head + control_expt + dpt_log)
    f.readline()
    for line in f:
        dpt_tmp.append(float(line))
    f.close()
    # Throw away first 2 repetitions
    dpt_tmp = dpt_tmp[control_skipyears*peryear:]
    # Calculate annual averages
    for year in range(num_years):
        dpt[-1,year] = mean(array(dpt_tmp[peryear*year:peryear*(year+1)]))

    # Weddell Sea Gyre and Ross Sea Gyre transport, already annually averaged
    ws_trans = empty([num_rcps+1, num_years])
    rs_trans = empty([num_rcps+1, num_years])
    # First read control experiment (correct years already selected)
    f = open(directory_head + control_expt + subpolar_log, 'r')
    f.readline()
    t = 0
    for line in f:
        try:
            ws_trans[-1, t] = float(line)
            t += 1
        except(ValueError):
            break
    t = 0
    for line in f:
        rs_trans[-1, t] = float(line)
        t += 1
    f.close()
    # Loop over RCP experiments
    for expt in range(num_rcps):
        # Copy the beginning of the control
        ws_trans[expt, 0:num_years_rep3] = ws_trans[-1, 0:num_years_rep3]
        rs_trans[expt, 0:num_years_rep3] = rs_trans[-1, 0:num_years_rep3]
        # Now read logfile for the rest of it
        f = open(directory_head + rcp_expt[expt] + subpolar_log, 'r')
        f.readline()
        t = num_years_rep3
        for line in f:
            try:
                ws_trans[expt, t] = float(line)
                t += 1
            except(ValueError):
                break
        t = num_years_rep3
        for line in f:
            rs_trans[expt, t] = float(line)
            t += 1
        f.close()

    # Plot
    fig = figure(figsize=(15,6))
    gs = GridSpec(1,3)
    gs.update(left=0.05, right=0.95, bottom=0.15, top=0.92)
    # Drake Passage transport
    ax = subplot(gs[0,0])
    # One line for each RCP
    for expt in range(num_rcps):
        ax.plot(time, dpt[expt,:], color=rcp_colours[expt], linewidth=1.8)
    # One line for control experiment
    ax.plot(time, dpt[-1,:], color=control_colour, linewidth=1.8)
    title('Drake Passage Transport', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel('Sv', fontsize=14)
    xlim([year_start, year_end])
    grid(True)
    # Weddell Sea Gyre transport
    ax = subplot(gs[0,1])
    for expt in range(num_rcps):
        ax.plot(time, ws_trans[expt,:], color=rcp_colours[expt], linewidth=1.8)
    ax.plot(time, ws_trans[-1,:], color=control_colour, linewidth=1.8)
    title('Weddell Sea Gyre Transport', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel('Sv', fontsize=14)
    xlim([year_start, year_end])
    grid(True)
    # Ross Sea Gyre transport
    ax = subplot(gs[0,2])
    for expt in range(num_rcps):
        # Add labels now
        ax.plot(time, rs_trans[expt,:], color=rcp_colours[expt], label=rcp_titles[expt], linewidth=1.8)
    ax.plot(time, rs_trans[-1,:], color=control_colour, label=control_title, linewidth=1.8)
    title('Ross Sea Gyre Transport', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel('Sv', fontsize=14)
    xlim([year_start, year_end])
    grid(True)
    subplots_adjust(wspace=0.15)
    # Add legend at bottom
    ax.legend(bbox_to_anchor=(0.5,-0.09), ncol=5, fontsize=14)
    fig.show()
    fig.savefig('transport.png')


# Command-line interface
if __name__ == "__main__":

    transport_3pt()
    
    
            
