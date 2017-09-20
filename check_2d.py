from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.patches import Polygon
from fesom_grid_2d import *

def check_2d (mesh_path):

    lat_max = -60+90

    elements = fesom_grid_2d(mesh_path, True)
    elements_global = fesom_grid_2d(mesh_path, False)
    patches = []
    shelf_patches = []
    ocean_patches = []
    for elm in elements:
        coord = transpose(vstack((elm.x, elm.y)))
        if elm.cavity:
            shelf_patches.append(Polygon(coord, True, linewidth=0.))
        else:
            ocean_patches.append(Polygon(coord, True, linewidth=0.))
        patches.append(Polygon(coord, True, linewidth=0.))
    patches_global = []
    for elm in elements_global:
        coord = transpose(vstack((elm.x, elm.y)))
        patches_global.append(Polygon(coord, True, linewidth=0.))

    node_shelf = []
    f = open(mesh_path + 'shelf.out', 'r')
    for line in f:
        node_shelf.append(-1*float(line))
    f.close()

    node_depth = []
    f = open(mesh_path + 'depth.out', 'r')
    for line in f:
        node_depth.append(-1*float(line))

    elm_shelf = []
    elm_depth = []
    elm_res = []
    for elm in elements:
        if elm.cavity:
            elm_shelf.append(mean(array([node_shelf[elm.nodes[0].id], node_shelf[elm.nodes[1].id], node_shelf[elm.nodes[2].id]])))
        elm_depth.append(mean(array([node_depth[elm.nodes[0].id], node_depth[elm.nodes[1].id], node_depth[elm.nodes[2].id]])))
        elm_res.append(sqrt(elm.area())*1e-3)

    elm_res_global = []
    for elm in elements_global:
        elm_res_global.append(sqrt(elm.area())*1e-3)

    x_reg, y_reg = meshgrid(linspace(-lat_max, lat_max, num=100), linspace(-lat_max, lat_max, num=100))
    land_square = zeros(shape(x_reg))

    fig = figure(figsize=(128, 96))
    ax = fig.add_subplot(1,1,1, aspect='equal')
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(shelf_patches, cmap='jet')
    img.set_array(array(elm_shelf))
    img.set_edgecolor('face')
    img.set_clim(vmin=0, vmax=2300)
    ax.add_collection(img)
    overlay = PatchCollection(ocean_patches, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    axis('off')
    title('Ice shelf draft (m)', fontsize=240)
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=160)
    fig.savefig('shelf.png')

    fig = figure(figsize=(128, 96))
    ax = fig.add_subplot(1,1,1, aspect='equal')
    contourf(x_reg, y_reg, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches, cmap='jet')
    img.set_array(array(elm_depth))
    img.set_edgecolor('face')
    img.set_clim(vmin=0, vmax=2800)
    ax.add_collection(img)
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    axis('off')
    title('Bathymetry (m)', fontsize=240)
    cbar = colorbar(img, extend='max')
    cbar.ax.tick_params(labelsize=160)
    fig.savefig('depth.png')

    fig = figure(figsize=(128, 96))
    ax = fig.add_subplot(1,1,1, aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(array(elm_res))
    img.set_edgecolor('face')
    img.set_clim(vmin=0, vmax=10)
    ax.add_collection(img)
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    axis('off')
    title('Horizontal resolution (km)', fontsize=240)
    cbar = colorbar(img, extend='max')
    cbar.ax.tick_params(labelsize=160)
    fig.savefig('res.png')

    fig = figure(figsize=(16,8))
    ax = fig.add_subplot(1,1,1)
    img = PatchCollection(patches_global, cmap='jet')
    img.set_array(array(elm_res_global))
    img.set_edgecolor('face')
    img.set_clim(vmin=0, vmax=225)
    ax.add_collection(img)
    xlim([-180, 180])
    ylim([-90, 90])
    ax.get_xaxis().set_ticks(arange(-120,120+1,60))
    ax.get_yaxis().set_ticks(arange(-60,60+1,30))
    title('Horizontal resolution (km)', fontsize=30)
    cbar = colorbar(img, extend='max')
    cbar.ax.tick_params(labelsize=20)
    fig.savefig('res_global.png')

if __name__ == "__main__":

    mesh_path = raw_input("Path to FESOM mesh directory: ")
    check_2d(mesh_path)
