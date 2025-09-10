from numpy import *
from matplotlib.pyplot import *

def transport_3pt ():

    # Paths to RCP experiment directories
    directory_head = '/short/y99/kaa561/FESOM/'
    rcp_expt = ['rcp45_M/', 'rcp45_A/', 'rcp85_M/', 'rcp85_A/']
    num_rcps = len(rcp_expt)
    # Titles for plot
    rcp_titles = ['RCP 4.5 M', 'RCP 4.5 A', 'RCP 8.5 M', 'RCP 8.5 A']
    # Colours for plotting
    rcp_colours = [(0, 0.4, 1), (0.6, 0.2, 1), (0, 0.6, 0), (1, 0, 0.4)]
    # Orders to plot
    zorders = [4, 2, 5, 3, 1]
    # Paths to control experiment directories
    control_expt = 'highres_spinup/'
    # Title for plot
    control_title = 'CONTROL'
    # Colour for plotting
    control_colour = (0,0,0)
    # Years to plot
    year_start = 1992
    rcp_year_start = 2006
    year_end = 2100
    # Trends in transport for each scenario
    #dpt_trends = ['-0.14 %/y', '-0.18 %/y', '-0.14 %/y', '-0.16 %/y', '-0.06 %/y']
    #wsg_trends = ['n/a', '+0.07 %/y', 'n/a', '-0.15 %/y', '-0.11 %/y']
    #rsg_trends = ['n/a', '+0.10 %/y', '+0.18 %/y', '+0.15 %/y', 'n/a']

    # Subpolar gyre parameters
    subpolar_log = 'subpolar_gyres.log'
    # Drake Passage transport parameters
    dpt_log = 'dpt.log'
    # Spinup years to discard (first 2 repetitions of 1992-2005)
    control_skipyears = 28
    # Output steps per year
    peryear = 365//5

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
    fig = figure(figsize=(15,5))
    gs = GridSpec(1,3)
    gs.update(left=0.05, right=0.95, bottom=0.17, top=0.92)
    # Drake Passage transport
    ax = subplot(gs[0,0])
    # One line for each RCP
    for expt in range(num_rcps):
        ax.plot(time, dpt[expt,:], color=rcp_colours[expt], linewidth=1.8, zorder=zorders[expt])
    # One line for control experiment
    ax.plot(time, dpt[-1,:], color=control_colour, linewidth=1.8, zorder=zorders[-1])
    title('a) Drake Passage Transport', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel('Sv', fontsize=14)
    xlim([year_start, year_end])
    #ylim([125, 170])
    grid(True)
    # Text labels showing trends
    #for expt in range(num_rcps):
        #text(1997, 146-expt*2.5, dpt_trends[expt], color=rcp_colours[expt], fontsize=13)
    #text(1997, 146-num_rcps*2.5, dpt_trends[-1], color=control_colour, fontsize=13)
    # Weddell Sea Gyre transport
    ax = subplot(gs[0,1])
    for expt in range(num_rcps):
        ax.plot(time, ws_trans[expt,:], color=rcp_colours[expt], linewidth=1.8, zorder=zorders[expt])
    ax.plot(time, ws_trans[-1,:], color=control_colour, linewidth=1.8, zorder=zorders[-1])
    title('b) Weddell Sea Gyre Transport', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel('Sv', fontsize=14)
    xlim([year_start, year_end])
    #ylim([50, 92])
    grid(True)
    #for expt in range(num_rcps):
        #text(1997, 89-expt*2.4, wsg_trends[expt], color=rcp_colours[expt], fontsize=13)
    #text(1997, 89-num_rcps*2.4, wsg_trends[-1], color=control_colour, fontsize=13)
    # Ross Sea Gyre transport
    ax = subplot(gs[0,2])
    for expt in range(num_rcps):
        # Add labels now
        ax.plot(time, rs_trans[expt,:], color=rcp_colours[expt], label=rcp_titles[expt], linewidth=1.8, zorder=zorders[expt])
    ax.plot(time, rs_trans[-1,:], color=control_colour, label=control_title, linewidth=1.8, zorder=zorders[-1])
    title('c) Ross Sea Gyre Transport', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel('Sv', fontsize=14)
    xlim([year_start, year_end])
    #ylim([33,60])
    grid(True)
    #for expt in range(num_rcps):
        #text(1997, 58.5-expt*1.5, rsg_trends[expt], color=rcp_colours[expt], fontsize=13)
        #text(1997, 58.5-num_rcps*1.5, rsg_trends[-1], color=control_colour, fontsize=13)
    subplots_adjust(wspace=0.15)
    # Add legend at bottom
    ax.legend(bbox_to_anchor=(0.5,-0.09), ncol=5, fontsize=14)
    fig.show()
    fig.savefig('transport.png')


# Command-line interface
if __name__ == "__main__":

    transport_3pt()
    
    
            
