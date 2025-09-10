from zonal_ts_before_after import *

def rcp_zonal_slices ():

    lon0 = [-59, -40, -1, 71, 140, 180, -159, -104, -83]
    lat_min = [-83, -82, -74, -73, -67, -85, -85, -75.5, -74]
    lat_max = [-70, -75, -68, -66, -64, -75, -73, -71, -69]
    num_slices = len(lon0)

    for rcp in ['45', '85']:
        for model in ['M', 'A']:
            print('Processing RCP ' + rcp[0] + '.' + rcp[1] + ' ' + model)
            for i in range(num_slices):
                if lon0[i] < 0:
                    lon_string = str(-lon0[i]) + 'W'
                else:
                    lon_string = str(lon0[i]) + 'E'
                print('...' + lon_string)
                zonal_ts_before_after(lon0[i], lat_min[i], lat_max[i], rcp, model, True, lon_string + '_rcp' + rcp + '_' + model + '.png')


if __name__ == "__main__":

    rcp_zonal_slices()
