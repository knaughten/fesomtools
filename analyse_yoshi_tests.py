# Example of a simple script to call fesom_vis routines manually
# This makes bottom water lon-lat plots of temperature and salinity, as well as
# lat-depth plots of temperature and salinity at 10 degree longitude
# increments, at the last timestep of each of six simulations (3 different 
# meshes x 2 different mixing schemes)

from patches import *
from circumpolar_plot import *
from zonal_slice_plot import *

# Bottom water lon-lat plots

#  Original mesh
#elements, patches = make_patches('/short/y99/kaa561/FESOM/mesh/fesom_grid_4.4.4/', True)
#   Original mixing
#    Temperature
#print 'Original mesh, original mixing, bottom water temperature'
#circumpolar_plot('/short/y99/kaa561/FESOM/yoshi_tests/origsigma_origmixing/output/MK44005.1995.oce.mean.nc', 'temp', 1, NaN, 3, elements, patches, False, True, 'origsigma_origmixing_bottomtemp.pdf')
#    Salinity
#print 'Original mesh, original mixing, bottom water salinity'
#circumpolar_plot('/short/y99/kaa561/FESOM/yoshi_tests/origsigma_origmixing/output/MK44005.1995.oce.mean.nc', 'salt', 1, NaN, 3, elements, patches, False, True, 'origsigma_origmixing_bottomsalt.pdf')
#   New mixing
#    Temperature
#print 'Original mesh, new mixing, bottom water temperature'
#circumpolar_plot('/short/y99/kaa561/FESOM/yoshi_tests/origsigma_newmixing/output/MK44005.1995.oce.mean.nc', 'temp', 1, NaN, 3, elements, patches, False, True, 'origsigma_newmixing_bottomtemp.pdf')
#    Salinity
#print 'Original mesh, new mixing, bottom water salinity'
#circumpolar_plot('/short/y99/kaa561/FESOM/yoshi_tests/origsigma_newmixing/output/MK44005.1995.oce.mean.nc', 'salt', 1, NaN, 3, elements, patches, False, True, 'origsigma_newmixing_bottomsalt.pdf')

#  Shallow sigma mesh
#elements, patches = make_patches('/short/y99/kaa561/FESOM/yoshi_grids/ice2sea5.0.37/', True)
#   Original mixing
#    Temperature
#print 'Shallow sigma mesh, original mixing, bottom water temperature'
#circumpolar_plot('/short/y99/kaa561/FESOM/yoshi_tests/shallowsigma_origmixing/output/MK44005.1995.oce.mean.nc', 'temp', 1, NaN, 3, elements, patches, False, True, 'shallowsigma_origmixing_bottomtemp.pdf')
#    Salinity
#print 'Shallow sigma mesh, original mixing, bottom water salinity'
#circumpolar_plot('/short/y99/kaa561/FESOM/yoshi_tests/shallowsigma_origmixing/output/MK44005.1995.oce.mean.nc', 'salt', 1, NaN, 3, elements, patches, False, True, 'shallowsigma_origmixing_bottomsalt.pdf')
#   New mixing
#    Temperature
#print 'Shallow sigma mesh, new mixing, bottom water temperature'
#circumpolar_plot('/short/y99/kaa561/FESOM/yoshi_tests/shallowsigma_newmixing/output/MK44005.1995.oce.mean.nc', 'temp', 1, NaN, 3, elements, patches, False, True, 'shallowsigma_newmixing_bottomtemp.pdf')
#    Salinity
#print 'Shallow sigma mesh, new mixing, bottom water salinity'
#circumpolar_plot('/short/y99/kaa561/FESOM/yoshi_tests/shallowsigma_newmixing/output/MK44005.1995.oce.mean.nc', 'salt', 1, NaN, 3, elements, patches, False, True, 'shallowsigma_newmixing_bottomsalt.pdf')

#  Deep sigma mesh
#elements, patches = make_patches('/short/y99/kaa561/FESOM/yoshi_grids/ice2sea5.0.38_deep/', True)
#   Original mixing
#    Temperature
#print 'Deep sigma mesh, original mixing, bottom water temperature'
#circumpolar_plot('/short/y99/kaa561/FESOM/yoshi_tests/deepsigma_origmixing/output/MK44005.1995.oce.mean.nc', 'temp', 1, NaN, 3, elements, patches, False, True, 'deepsigma_origmixing_bottomtemp.pdf')
#    Salinity
#print 'Deep sigma mesh, original mixing, bottom water salinity'
#circumpolar_plot('/short/y99/kaa561/FESOM/yoshi_tests/deepsigma_origmixing/output/MK44005.1995.oce.mean.nc', 'salt', 1, NaN, 3, elements, patches, False, True, 'deepsigma_origmixing_bottomsalt.pdf')
#   New mixing
#    Temperature
#print 'Deep sigma mesh, new mixing, bottom water temperature'
#circumpolar_plot('/short/y99/kaa561/FESOM/yoshi_tests/deepsigma_newmixing/output/MK44005.1995.oce.mean.nc', 'temp', 1, NaN, 3, elements, patches, False, True, 'deepsigma_newmixing_bottomtemp.pdf')
#    Salinity
#print 'Deep sigma mesh, new mixing, bottom water salinity'
#circumpolar_plot('/short/y99/kaa561/FESOM/yoshi_tests/deepsigma_newmixing/output/MK44005.1995.oce.mean.nc', 'salt', 1, NaN, 3, elements, patches, False, True, 'deepsigma_newmixing_bottomsalt.pdf')

# Zonal slice plots

#  Original mesh
#   Original mixing
#    Temperature
for lon in range(-180, 180, 10):
    if lon < 0:
        fig_name = 'origsigma_origmixing_tempslice_'+str(lon)+'W.pdf'
    else:
        fig_name = 'origsigma_origmixing_tempslice_'+str(lon)+'E.pdf'
    print 'Original mesh, original mixing, temperature at longitude '+str(lon)
    zonal_slice_plot('/short/y99/kaa561/FESOM/mesh/fesom_grid_4.4.4/', '/short/y99/kaa561/FESOM/yoshi_tests/origsigma_origmixing/output/MK44005.1995.oce.mean.nc', 'temp', 3, lon, -2000, True, fig_name)
#    Salinity
for lon in range(-180, 180, 10):
    if lon < 0:
        fig_name = 'origsigma_origmixing_saltslice_'+str(lon)+'W.pdf'
    else:
        fig_name = 'origsigma_origmixing_saltslice_'+str(lon)+'E.pdf'
    print 'Original mesh, original mixing, salinity at longitude '+str(lon)
    zonal_slice_plot('/short/y99/kaa561/FESOM/mesh/fesom_grid_4.4.4/', '/short/y99/kaa561/FESOM/yoshi_tests/origsigma_origmixing/output/MK44005.1995.oce.mean.nc', 'salt', 3, lon, -2000, True, fig_name)
#   New mixing
#    Temperature
for lon in range(-180, 180, 10):
    if lon < 0:
        fig_name = 'origsigma_newmixing_tempslice_'+str(lon)+'W.pdf'
    else:
        fig_name = 'origsigma_newmixing_tempslice_'+str(lon)+'E.pdf'
    print 'Original mesh, new mixing, temperature at longitude '+str(lon)
    zonal_slice_plot('/short/y99/kaa561/FESOM/mesh/fesom_grid_4.4.4/', '/short/y99/kaa561/FESOM/yoshi_tests/origsigma_newmixing/output/MK44005.1995.oce.mean.nc', 'temp', 3, lon, -2000, True, fig_name)
#    Salinity
for lon in range(-180, 180, 10):
    if lon < 0:
        fig_name = 'origsigma_newmixing_saltslice_'+str(lon)+'W.pdf'
    else:
        fig_name = 'origsigma_newmixing_saltslice_'+str(lon)+'E.pdf'
    print 'Original mesh, new mixing, salinity at longitude '+str(lon)
    zonal_slice_plot('/short/y99/kaa561/FESOM/mesh/fesom_grid_4.4.4/', '/short/y99/kaa561/FESOM/yoshi_tests/origsigma_newmixing/output/MK44005.1995.oce.mean.nc', 'salt', 3, lon, -2000, True, fig_name)

#  Shallow sigma mesh
#   Original mixing
#    Temperature
for lon in range(-180, 180, 10):
    if lon < 0:
        fig_name = 'shallowsigma_origmixing_tempslice_'+str(lon)+'W.pdf'
    else:
        fig_name = 'shallowsigma_origmixing_tempslice_'+str(lon)+'E.pdf'
    print 'Shallow sigma mesh, original mixing, temperature at longitude '+str(lon)
    zonal_slice_plot('/short/y99/kaa561/FESOM/yoshi_grids/ice2sea5.0.37/', '/short/y99/kaa561/FESOM/yoshi_tests/shallowsigma_origmixing/output/MK44005.1995.oce.mean.nc', 'temp', 3, lon, -2000, True, fig_name)
#    Salinity
for lon in range(-180, 180, 10):
    if lon < 0:
        fig_name = 'shallowsigma_origmixing_saltslice_'+str(lon)+'W.pdf'
    else:
        fig_name = 'shallowsigma_origmixing_saltslice_'+str(lon)+'E.pdf'
    print 'Shallow sigma mesh, original mixing, salinity at longitude '+str(lon)
    zonal_slice_plot('/short/y99/kaa561/FESOM/yoshi_grids/ice2sea5.0.37/', '/short/y99/kaa561/FESOM/yoshi_tests/shallowsigma_origmixing/output/MK44005.1995.oce.mean.nc', 'salt', 3, lon, -2000, True, fig_name)
#   New mixing
#    Temperature
for lon in range(-180, 180, 10):
    if lon < 0:
        fig_name = 'shallowsigma_newmixing_tempslice_'+str(lon)+'W.pdf'
    else:
        fig_name = 'shallowsigma_newmixing_tempslice_'+str(lon)+'E.pdf'
    print 'Shallow sigma mesh, new mixing, temperature at longitude '+str(lon)
    zonal_slice_plot('/short/y99/kaa561/FESOM/yoshi_grids/ice2sea5.0.37/', '/short/y99/kaa561/FESOM/yoshi_tests/shallowsigma_newmixing/output/MK44005.1995.oce.mean.nc', 'temp', 3, lon, -2000, True, fig_name)
#    Salinity
for lon in range(-180, 180, 10):
    if lon < 0:
        fig_name = 'shallowsigma_newmixing_saltslice_'+str(lon)+'W.pdf'
    else:
        fig_name = 'shallowsigma_newmixing_saltslice_'+str(lon)+'E.pdf'
    print 'Shallow sigma mesh, new mixing, salinity at longitude '+str(lon)
    zonal_slice_plot('/short/y99/kaa561/FESOM/yoshi_grids/ice2sea5.0.37/', '/short/y99/kaa561/FESOM/yoshi_tests/shallowsigma_newmixing/output/MK44005.1995.oce.mean.nc', 'salt', 3, lon, -2000, True, fig_name)

#  Deep sigma mesh
#   Original mixing
#    Temperature
for lon in range(-180, 180, 10):
    if lon < 0:
        fig_name = 'deepsigma_origmixing_tempslice_'+str(lon)+'W.pdf'
    else:
        fig_name = 'deepsigma_origmixing_tempslice_'+str(lon)+'E.pdf'
    print 'Deep sigma mesh, original mixing, temperature at longitude '+str(lon)
    zonal_slice_plot('/short/y99/kaa561/FESOM/yoshi_grids/ice2sea5.0.38_deep/', '/short/y99/kaa561/FESOM/yoshi_tests/deepsigma_origmixing/output/MK44005.1995.oce.mean.nc', 'temp', 3, lon, -2000, True, fig_name)
#    Salinity
for lon in range(-180, 180, 10):
    if lon < 0:
        fig_name = 'deepsigma_origmixing_saltslice_'+str(lon)+'W.pdf'
    else:
        fig_name = 'deepsigma_origmixing_saltslice_'+str(lon)+'E.pdf'
    print 'Deep sigma mesh, original mixing, salinity at longitude '+str(lon)
    zonal_slice_plot('/short/y99/kaa561/FESOM/yoshi_grids/ice2sea5.0.38_deep/', '/short/y99/kaa561/FESOM/yoshi_tests/deepsigma_origmixing/output/MK44005.1995.oce.mean.nc', 'salt', 3, lon, -2000, True, fig_name)

#   New mixing
#    Temperature
for lon in range(-180, 180, 10):
    if lon < 0:
        fig_name = 'deepsigma_newmixing_tempslice_'+str(lon)+'W.pdf'
    else:
        fig_name = 'deepsigma_newmixing_tempslice_'+str(lon)+'E.pdf'
    print 'Deep sigma mesh, new mixing, temperature at longitude '+str(lon)
    zonal_slice_plot('/short/y99/kaa561/FESOM/yoshi_grids/ice2sea5.0.38_deep/', '/short/y99/kaa561/FESOM/yoshi_tests/deepsigma_newmixing/output/MK44005.1995.oce.mean.nc', 'temp', 3, lon, -2000, True, fig_name)
#    Salinity
for lon in range(-180, 180, 10):
    if lon < 0:
        fig_name = 'deepsigma_newmixing_saltslice_'+str(lon)+'W.pdf'
    else:
        fig_name = 'deepsigma_newmixing_saltslice_'+str(lon)+'E.pdf'
    print 'Deep sigma mesh, new mixing, salinity at longitude '+str(lon)
    zonal_slice_plot('/short/y99/kaa561/FESOM/yoshi_grids/ice2sea5.0.38_deep/', '/short/y99/kaa561/FESOM/yoshi_tests/deepsigma_newmixing/output/MK44005.1995.oce.mean.nc', 'salt', 3, lon, -2000, True, fig_name)
