from numpy import *
from matplotlib.pyplot import *
from unesco import *

def watermass_schematic ():

    min_salt = 33.7
    max_salt = 34.8
    min_temp = -2.4
    max_temp = 0.5
    salt_ticks = [34, 34.5]
    salt_tick_labels = ['34', '34.5']
    temp_ticks = [tfreeze(min_salt), -1.5, 0]
    temp_tick_labels = [r'$T_f$', '-1.5', '0']

    # Get potential density
    salt_vals = linspace(min_salt, max_salt, num=500)
    temp_vals = linspace(min_temp, max_temp, num=500)
    salt_2d, temp_2d = meshgrid(salt_vals, temp_vals)
    density = unesco(temp_2d, salt_2d, zeros(shape(temp_2d)))-1000

    fig = figure(figsize=(8,6))
    gs = GridSpec(1,1)
    gs.update(left=0.2, right=0.9, bottom=0.2, top=0.9)
    ax = subplot(gs[0,0])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    # First contour potential density
    contour(salt_vals, temp_vals, density, colors='black', linewidth=1, linestyles='dotted')
    # Draw all the lines
    plot([min_salt, max_salt], [tfreeze(min_salt), tfreeze(max_salt)], color='black', linewidth=2, linestyle='dashed')
    arrow(max_salt*0.999, tfreeze(max_salt*0.999), max_salt*0.001, tfreeze(max_salt)-tfreeze(max_salt*0.999), head_width=0.09, head_length=0.025, fc='k', ec='k', length_includes_head=True, clip_on=False)
    arrow(min_salt*1.001, tfreeze(min_salt*1.001), -0.001*min_salt, tfreeze(min_salt)-tfreeze(min_salt*1.001), head_width=0.09, head_length=0.025, fc='k', ec='k', length_includes_head=True, clip_on=False)
    plot([34, 34], [tfreeze(34), max_temp], color='black', linewidth=2)
    arrow(34, max_temp*0.95, 0, max_temp*0.05, head_width=0.025, head_length=0.1, fc='k', ec='k', length_includes_head=True, clip_on=False)
    plot([34, max_salt], [-1.5, -1.5], color='black', linewidth=2)
    arrow(max_salt*0.999, -1.5, max_salt*0.001, 0, head_width=0.09, head_length=0.025, fc='k', ec='k', length_includes_head=True, clip_on=False)
    plot([34.5, 34.5], [tfreeze(34.5), -1.5], color='black', linewidth=2)
    plot([34, max_salt], [0, 0], color='black', linewidth=2)
    arrow(max_salt*0.999, 0, max_salt*0.001, 0, head_width=0.09, head_length=0.025, fc='k', ec='k', length_includes_head=True, clip_on=False)
    # Label all the sections
    text(34.25, -2.1, 'ISW', fontsize=30, ha='center', va='center')
    text(33.85, -0.75, 'AASW', fontsize=30, ha='center', va='center')
    text(34.25, -1.7, 'LSSW', fontsize=30, ha='center', va='center')
    text(34.65, -1.7, 'HSSW', fontsize=30, ha='center', va='center')
    text(34.4, -0.75, 'MCDW', fontsize=30, ha='center', va='center')
    text(34.4, 0.25, 'CDW', fontsize=30, ha='center', va='center')
    # Set up axis labels
    ax.set_xticks(salt_ticks)
    ax.set_xticklabels(salt_tick_labels, fontsize=27)
    ax.set_yticks(temp_ticks)
    ax.set_yticklabels(temp_tick_labels, fontsize=27)
    xlabel('Salinity (psu)', fontsize=27)
    ylabel(r'Temperature ($^{\circ}$C)', fontsize=27)
    xlim([min_salt, max_salt])
    ylim([min_temp, max_temp])
    fig.show()
    fig.savefig('watermass_schematic.png')

    
# For the given salinity, return the surface freezing point as seen by the
# FESOM sea ice code.
def tfreeze (salt):
    return -0.0575*salt + 1.7105e-3*sqrt(salt**3) - 2.155e-4*salt**2


# Command-line interface
if __name__ == "__main__":

    watermass_schematic()
