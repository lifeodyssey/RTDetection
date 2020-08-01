from matplotlib.colors import ListedColormap
import pyresample
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib import colors
import math
from scipy.interpolate import RectSphereBivariateSpline


# # from multiprocess import Pool
# warnings.filterwarnings("ignore")  # 用来去掉Warning,
# os.environ['PROJ_LIB'] = '/Users/zhenjia/opt/anaconda3/share/proj'
# import pandas as pd
# import matplotlib.pyplot as plt
#
# # os.chdir('I:/Nagoya University/Project/Seto/MODIS')
# # path='I:/Nagoya University/Project/Seto/MODIS/'
# os.chdir('/Users/zhenjia/Desktop/Project/Seto/MODIS/l1/l2_lac')
# path = '/Users/zhenjia/Desktop/Project/Seto/MODIS/l1/l2_lac/'
# datalist = glob.glob('2008*.nc')
# # 20110715 Uwajima
# # 20180723 Osaka
# minlat = 32.5
# minlon = 130
# maxlat = 35
# maxlon = 136
# # area of full seto-inland sea
# # # # for a1 in range(1):
# # minlat = 32.32
# # minlon = 131.39
# # maxlat = 33.48
# # maxlon = 132.77
#
# # # bun go si i to
# # minlat = 33.88
# # minlon = 132.85
# # maxlat = 34.73
# # maxlon = 133.837
# # middle
# # minlat = 34.18
# # minlon = 134.00
# # maxlat = 34.88
# # maxlon = 135.59
# # # Osaka
# L = len(datalist)
# initime = datetime.datetime.now()
# #
# # # plt.subplot(2,3,3)
# # # with concurrent.futures.ProcessPoolExecutor() as executor:
# # a1 = np.arange(len(datalist))
# #
# #
# # #
# #
# #
# # def main(a1):
# file = datalist[0]
# filename = path + file
# nc_file = nc4.Dataset(filename, 'r')
# print(nc_file.time_coverage_end)
# lon = nc_file.groups['navigation_data'].variables['longitude'][:]
# lat = nc_file.groups['navigation_data'].variables['latitude'][:]
# variables = nc_file.groups['geophysical_data'].variables
#
# x = np.arange(minlon, maxlon, 0.01)  # 1 km grid,
# y = np.arange(maxlat, minlat, -0.01)
#
# for i in variables:
#     var = variables[i][:]
#     np.where(var <= 0, var, np.nan)
#     if i != 'l2_flags':
#         var_re, grid = gs.swath_resampling(var, lon, lat, x, y, 1000)  # 1 km grid
#         # var_re=var_re.filled()
#         if np.ma.is_mask(var_re):
#             var_re = var_re.filled(np.nan)
#             var_re[var_re == -32767.0] = np.nan
#         variables[i] = var_re
#     else:
#         # var_re = var_re.filled()
#         variables[i] = var
#
# lons = grid.lons
# lats = grid.lats
# Rrs_412 = variables['Rrs_412']
# Rrs_443 = variables['Rrs_443']
# Rrs_469 = variables['Rrs_469']
# Rrs_488 = variables['Rrs_488']
# Rrs_531 = variables['Rrs_531']
# Rrs_547 = variables['Rrs_547']
# Rrs_555 = variables['Rrs_555']
# Rrs_645 = variables['Rrs_645']
# Rrs_667 = variables['Rrs_667']
# Rrs_678 = variables['Rrs_678']
# Rrs_748 = variables['Rrs_748']
# # F0 = [172.912, 187.622, 205.878, 194.933, 185.747, 186.539, 183.869, 157.811, 152.255, 148.052, 128.065]
# # nlw412 = Rrs_412 * F0[0]
# # nlw443 = Rrs_443 * F0[1]
# # nlw469 = Rrs_469 * F0[2]
# # nlw488 = Rrs_488 * F0[3]
# # nlw531 = Rrs_531 * F0[4]
# # nlw547 = Rrs_547 * F0[5]
# # nlw555 = Rrs_555 * F0[6]
# # nlw645 = Rrs_645 * F0[7]
# # nlw667 = Rrs_667 * F0[8]
# # nlw678 = Rrs_678 * F0[9]
# # nlw748 = Rrs_748 * F0[10]
# # # nlw412 = variables['nLw_412']
# # # nlw443 = variables['nLw_443']
# # # nlw469 = variables['nLw_469']
# # # nlw488 = variables['nLw_488']
# # # nlw531 = variables['nLw_531']
# # # nlw547 = variables['nLw_547']
# # # nlw555 = variables['nLw_555']
# # # nlw645 = variables['nLw_645']
# # # nlw667 = variables['nLw_667']
# # # nlw678 = variables['nLw_678']
# # # nlw748 = variables['nLw_748']
# #
# # # Rrs_412 = variables['Rrs_412'].filled(np.nan)
# # # Rrs_443 = variables['Rrs_443'].filled(np.nan)
# # # Rrs_469 = variables['Rrs_469'].filled(np.nan)
# # # Rrs_488 = variables['Rrs_488'].filled(np.nan)
# # # Rrs_531 = variables['Rrs_531'].filled(np.nan)
# # # Rrs_547 = variables['Rrs_547'].filled(np.nan)
# # # Rrs_555 = variables['Rrs_555'].filled(np.nan)
# # # Rrs_645 = variables['Rrs_645'].filled(np.nan)
# # # Rrs_667 = variables['Rrs_667'].filled(np.nan)
# # # Rrs_678 = variables['Rrs_678'].filled(np.nan)
# # # Rrs_748 = variables['Rrs_748'].filled(np.nan)
# # # nlw412 = variables['nLw_412'].filled(np.nan)
# # # nlw443 = variables['nLw_443'].filled(np.nan)
# # # nlw469 = variables['nLw_469'].filled(np.nan)
# # # nlw488 = variables['nLw_488'].filled(np.nan)
# # # nlw531 = variables['nLw_531'].filled(np.nan)
# # # nlw547 = variables['nLw_547'].filled(np.nan)
# # # nlw555 = variables['nLw_555'].filled(np.nan)
# # # nlw645 = variables['nLw_645'].filled(np.nan)
# # # nlw667 = variables['nLw_667'].filled(np.nan)
# # # nlw678 = variables['nLw_678'].filled(np.nan)
# # # nlw748 = variables['nLw_748'].filled(np.nan)
# Rrs = np.array([
#     Rrs_412,
#     Rrs_443,
#     Rrs_469,
#     Rrs_488,
#     Rrs_531,
#     Rrs_547,
#     Rrs_555,
#     Rrs_645,
#     Rrs_667,
#     Rrs_678])
# # nflh = nlw678 - (70 / 81) * nlw667 - (11 / 81) * nlw748
# chl = variables['chlor_a']
#
# wt = np.ma.array(np.zeros(np.shape(chl)))
# wt[chl > 10] = 3
# wt[(chl < 10) & (chl > 5)] = 2
# wt[(chl < 5)] = 1
# wt[chl.mask] = 0
#

def plot_WT_image(sds: np.ma.array, lon: np.ndarray, lat: np.ndarray, title: str = None,
                  lon_range: list = None, lat_range: list = None, save_image: str = None,
                  dpi: int = 400):
    if len(lon.shape) == 1:
        print('MeshGridding...')
        lon, lat = np.meshgrid(lon, lat)

    lon_0 = (lon.min() + lon.max()) / 2
    lat_0 = (lat.min() + lat.max()) / 2

    print(f'Lat: [{lat.min():.3f}, {lat.max():.3f}] | '
          f'Lon: [{lon.min():.3f}, {lon.max():.3f}] | '
          f'SDS: [{sds.min():.3f}, {sds.max():.3f}]')

    if (lon_range is not None) and (lat_range is not None):
        m = Basemap(llcrnrlon=min(lon_range), llcrnrlat=min(lat_range),
                    urcrnrlon=max(lon_range), urcrnrlat=max(lat_range),
                    resolution='f', lon_0=lon_0, lat_0=lat_0, projection='merc')
    else:
        m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(),
                    urcrnrlon=lon.max(), urcrnrlat=lat.max(),
                    resolution='f', lon_0=lon_0, lat_0=lat_0, projection='merc')
    x2d, y2d = m(lon, lat)

    fig = plt.figure(figsize=(12, 12 * m.aspect))
    ax = fig.add_axes([0.08, 0.1, 0.7, 0.7], facecolor='w')

    if (lon_range is not None) and (lat_range is not None):
        parallels = np.arange(min(lat_range), max(lat_range), 0.5)
        meridians = np.arange(min(lon_range), max(lon_range), 0.5)
    else:
        parallels = np.arange(lat.min(), lat.max(), 0.5)
        meridians = np.arange(lon.min(), lon.max(), 0.5)
        # 0 nodata
        # 1 other water
        # 2 other phytoplankton
        # 3k,bloom
    labels = ['No Data',
              'Other Water',
              'Other Phytoplankton',
              'K.mikimotoi Bloom']
    colors = ['w',
              'deepskyblue',
              'g',
              'r'
              ]
    cmap = ListedColormap(colors)
    p = m.pcolor(x2d, y2d, sds, cmap=cmap)

    if title is not None:
        plt.title(title, fontsize=18)

    plt.plot([], [], label=labels[0], color=colors[0])
    plt.plot([], [], label=labels[1], color=colors[1])
    plt.plot([], [], label=labels[2], color=colors[2])
    plt.plot([], [], label=labels[3], color=colors[3])
    plt.legend(loc='lower left')

    m.drawcoastlines()
    m.drawcountries()
    m.fillcontinents(color='lightgray')
    m.drawrivers()

    m.drawmeridians(meridians, fontsize=10, linewidth=0.25, dashes=[7, 15],
                    color='k', labels=[1, 0, 1, 1])
    m.drawparallels(parallels, fontsize=10, dashes=[7, 15],
                    linewidth=0.3, color='k', labels=[1, 1, 0, 1])
    # plt.show()

    if save_image is not None:
        plt.savefig(save_image, dpi=dpi, facecolor='w', edgecolor='w', orientation='portrait')
        # plt.show()
        plt.close()


def swath_resampling(src_data: np.ma.array, src_lon: np.array, src_lat: np.array,
                     trg_lon: np.array, trg_lat: np.array, search_radius: float):
    if len(trg_lon.shape) == 1:
        grid_def = pyresample.geometry.SwathDefinition(*np.meshgrid(trg_lon, trg_lat))
    else:
        grid_def = pyresample.geometry.SwathDefinition(lons=trg_lon, lats=trg_lat)

    # source grid with original swath data
    # if len(src_lon.shape) == 1:
    #     swath_def = pyresample.geometry.SwathDefinition(*np.meshgrid(src_lon, src_lat,sparse=True))
    # else:
    #     swath_def = pyresample.geometry.SwathDefinition(lons=src_lon, lats=src_lat)

    swath_def = pyresample.geometry.SwathDefinition(lons=src_lon, lats=src_lat)
    # resample (here we use nearest. Bilinear, gaussian and custom defined methods are available)
    # for more, visit https://pyresample.readthedocs.io/en/latest/
    result = pyresample.kd_tree.resample_nearest(swath_def, src_data, grid_def, epsilon=0.5,
                                                 fill_value=np.nan, radius_of_influence=search_radius)
    return result, grid_def


# @concurrent
def plot_geo_image(sds: np.ma.array, lon: np.ndarray, lat: np.ndarray, log10: bool = True, title: str = None,
                   label: str = None,
                   caxis: list = None, lon_range: list = None, lat_range: list = None, save_image: str = None,
                   dpi: int = 400):
    if len(lon.shape) == 1:
        print('MeshGridding...')
        lon, lat = np.meshgrid(lon, lat)

    lon_0 = (lon.min() + lon.max()) / 2
    lat_0 = (lat.min() + lat.max()) / 2

    print(f'Lat: [{lat.min():.3f}, {lat.max():.3f}] | '
          f'Lon: [{lon.min():.3f}, {lon.max():.3f}] | '
          f'SDS: [{sds.min():.3f}, {sds.max():.3f}]')

    if (lon_range is not None) and (lat_range is not None):
        m = Basemap(llcrnrlon=min(lon_range), llcrnrlat=min(lat_range),
                    urcrnrlon=max(lon_range), urcrnrlat=max(lat_range),
                    resolution='f', lon_0=lon_0, lat_0=lat_0, projection='merc')
    else:
        m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(),
                    urcrnrlon=lon.max(), urcrnrlat=lat.max(),
                    resolution='f', lon_0=lon_0, lat_0=lat_0, projection='merc')
    x2d, y2d = m(lon, lat)

    fig = plt.figure(figsize=(12, 12 * m.aspect))
    ax = fig.add_axes([0.08, 0.1, 0.7, 0.7], facecolor='white')
    # changed to facecolor 8 October 2019

    if (lon_range is not None) and (lat_range is not None):
        parallels = np.arange(min(lat_range), max(lat_range), 3)
        meridians = np.arange(min(lon_range), max(lon_range), 4)
    else:
        parallels = meridians = None

    if caxis is not None:
        cmn, cmx = min(caxis), max(caxis)
    else:
        cmn, cmx = sds.min(), sds.max()
    # m.drawparallels(parallels, fontsize=10, linewidth=0.25, dashes=[7, 15],
    #                  color='k', labels=[1, 0, 1, 1])
    # m.drawmeridians(meridians, fontsize=10, dashes=[7, 15],
    #                  linewidth=0.3, color='k', labels=[1, 1, 0, 1])
    # ncl = 150
    # if log10 is True:
    #    norm = colors.LogNorm(vmin=cmn, vmax=cmx)
    # else:
    #   bounds = np.linspace(cmn, cmx, ncl)
    #   norm = colors.BoundaryNorm(boundaries=bounds, ncolors=ncl)

    p = m.pcolor(x2d, y2d, sds, vmin=cmn, vmax=cmx, cmap='jet')

    if title is not None:
        plt.title(title, fontsize=24)

    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes('vertical', size="3%", pad=0.05)
    # cax = plt.axes([cmn, 0, cmx])  # setup colorbar axes

    cb = m.colorbar(p, location="right", size="5%", pad=0.1)  # draw colorbar
    if label is not None:
        cb.set_label("%s" % label)
    plt.sca(ax)  # make the original axes current again
    plt.clim(cmn, cmx)
    unit='Elevation to the sea level'
    cb.set_label(unit, rotation=270, labelpad=10.0, fontsize=10)
    cb.ax.tick_params(labelsize=10)

    m.drawcoastlines()
    m.drawcountries()
    m.fillcontinents()
    # plt.show()

    if save_image is not None:
        plt.savefig(save_image, dpi=dpi, facecolor='w', edgecolor='w', orientation='portrait')

        # plt.show()
        plt.close()


def creategrid(min_lon, max_lon, min_lat, max_lat, cell_size_deg, mesh=False):
    # Output grid within geobounds and specifice cell size
    # cell_size_deg should be in decimal degrees’’’

    min_lon = math.floor(min_lon)
    max_lon = math.ceil(max_lon)
    min_lat = math.floor(min_lat)
    max_lat = math.ceil(max_lat)
    lon_num = (max_lon - min_lon) / cell_size_deg
    lat_num = (max_lat - min_lat) / cell_size_deg
    grid_lons = np.zeros(lon_num)  # fill with lon_min
    grid_lats = np.zeros(lat_num)  # fill with lon_max
    grid_lons = grid_lons + (np.assary(range(lon_num)) * cell_size_deg)
    grid_lats = grid_lats + (np.assary(range(lat_num)) * cell_size_deg)
    grid_lons, grid_lats = np.meshgrid(grid_lons, grid_lats)
    grid_lons = np.ravel(grid_lons)
    grid_lats = np.ravel(grid_lats)
    # if mesh = True:
    # grid_lons = grid_lons
    # grid_lats = grid_lats
    return grid_lons, grid_lats


def geointerp(lats, lons, data, grid_size_deg, mesh=False):
    # We want to interpolate it to a global x-degree grid’’’
    deg2rad = np.pi / 180.
    new_lats = np.linspace(grid_size_deg, 180, 180 / grid_size_deg)
    new_lons = np.linspace(grid_size_deg, 360, 360 / grid_size_deg)
    new_lats_mesh, new_lons_mesh = np.meshgrid(new_lats * deg2rad, new_lons * deg2rad)
    # We need to set up the interpolator object’’’
    lut = RectSphereBivariateSpline(lats * deg2rad, lons * deg2rad, data)
    # Finally we interpolate the data. The RectSphereBivariateSpline
    # object only takes 1-D arrays as input, therefore we need to do some reshaping.’’’
    new_lats = new_lats_mesh.ravel()
    new_lons = new_lons_mesh.ravel()
    data_interp = lut.ev(new_lats, new_lons)
    if mesh == True:
        data_interp = data_interp.reshape((360 / grid_size_deg,
                                           180 / grid_size_deg)).T
        return new_lats / deg2rad, new_lons / deg2rad, data_interp


def plot_chl(sds: np.ma.array, lon: np.ndarray, lat: np.ndarray, log10: bool = True, title: str = None,
             label: str = 'CHL [mg/m$^3$]',
             caxis: list = None, lon_range: list = None, lat_range: list = None, save_image: str = None,
             dpi: int = 400):
    if len(lon.shape) == 1:
        print('MeshGridding...')
        lon, lat = np.meshgrid(lon, lat)

    lon_0 = (lon.min() + lon.max()) / 2
    lat_0 = (lat.min() + lat.max()) / 2

    print(f'Lat: [{lat.min():.3f}, {lat.max():.3f}] | '
          f'Lon: [{lon.min():.3f}, {lon.max():.3f}] | '
          f'SDS: [{sds.min():.3f}, {sds.max():.3f}]')

    if (lon_range is not None) and (lat_range is not None):
        m = Basemap(llcrnrlon=min(lon_range), llcrnrlat=min(lat_range),
                    urcrnrlon=max(lon_range), urcrnrlat=max(lat_range),
                    resolution='f', lon_0=lon_0, lat_0=lat_0, projection='merc')
    else:
        m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(),
                    urcrnrlon=lon.max(), urcrnrlat=lat.max(),
                    resolution='f', lon_0=lon_0, lat_0=lat_0, projection='merc')
    x2d, y2d = m(lon, lat)

    fig = plt.figure(figsize=(8, 8 * m.aspect))
    ax = fig.add_axes([0.08, 0.1, 0.7, 0.7], facecolor='w')
    # changed to facecolor 8 October 2019

    if (lon_range is not None) and (lat_range is not None):
        parallels = np.linspace(min(lat_range), max(lat_range), 0.5)
        meridians = np.linspace(min(lon_range), max(lon_range), 0.5)
    else:
        parallels = np.arange(lat.min(), lat.max(), 0.5)
        meridians = np.arange(lon.min(), lon.max(), 0.5)

    if caxis is not None:
        cmn, cmx = min(caxis), max(caxis)
    else:
        cmn, cmx = sds.min(), sds.max()

    ncl = 150
    if log10 is True:
        norm = colors.LogNorm(vmin=cmn, vmax=cmx)
    else:
        bounds = np.linspace(cmn, cmx, ncl)
        norm = colors.BoundaryNorm(boundaries=bounds, ncolors=ncl)

    p = m.pcolor(x2d, y2d, sds, norm=norm, cmap=plt.cm.jet)

    if title is not None:
        plt.title(title)

    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes('vertical', size="3%", pad=0.05)
    cax = plt.axes([0.85, 0.1, 0.05, 0.7])  # setup colorbar axes

    cb = plt.colorbar(format='%5.2f', cax=cax)  # draw colorbar
    if label is not None:
        cb.set_label("%s" % label)
    plt.sca(ax)  # make the original axes current again
    plt.clim(cmn, cmx)
    m.fillcontinents(color='lightgrey')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmeridians(meridians, fontsize=10, linewidth=0.25, dashes=[7, 15],
                    color='k', labels=[1, 0, 1, 1])
    m.drawparallels(parallels, fontsize=10, dashes=[7, 15],
                    linewidth=0.3, color='k', labels=[1, 1, 0, 1])
    m.drawrivers()
    # plt.show()

    if save_image is not None:
        plt.savefig(save_image, dpi=dpi, facecolor='w', edgecolor='w', orientation='portrait')
        # plt.show()
        plt.close()


def QAAv6MODIS(Rrs):
    # acording to Lee,QAAv6
    # write by Zhou,20191130
    # Input data need to be an arrary that contain 10 bands Rrs of MODIS,from short wavelength to long wavelength in a certain station
    # Output is a tuple, first array is aph,second array is bbp
    # use as import QAAV6
    # B1: 412 nm 0
    # B2: 443 nm 1
    # B3: 469 nm 2
    # B4: 488 nm 3
    # B5: 531 nm 4
    # B6: 547 nm 5
    # B7: 555 nm 6
    # B8: 645 nm 7
    # B9: 667 nm 8
    # B10: 678 nm 9

    Lambda = np.array([412, 443, 469, 488, 531, 547, 555, 645, 667, 678])
    nbands = np.shape(Lambda)[0]

    IOPw = np.array([[0.003344468, 0.004572564],
                     [0.00244466, 0.00635],
                     [0.001910803, 0.010483637],
                     [0.001609567, 0.014361745],
                     [0.00111757, 0.043747657],
                     [0.000983055, 0.053262848],
                     [0.000923288, 0.0595],
                     [0.000482375, 0.325],
                     [0.00041731, 0.433497423],
                     [0.00038884, 0.457440162]])

    #     if(Rrs=np.nan):
    #     return np.nan
    # else:
    # bbw from Morel (1974).aw  from Pope and Fry (1997)
    bbp = np.ones(10)
    adg = np.ones(10)
    if (np.nan in Rrs):
        bbp[:] = np.nan
    else:
        bw = IOPw[:, 0]  # backscaterring of pure water
        aw = IOPw[:, 1]  # absorption of pure water
        rrs = Rrs / (0.52 + 1.7 * Rrs)
        g0 = 0.089
        g1 = 0.1245
        u = (-g0 + ((g0 ** 2) + 4 * g1 * rrs) ** 0.5) / (2 * g1)

        aph = np.ones(10)  # adg is the absorption of CDOM and NAP
        if Rrs[6] < 0.0015:  # select 555 as reference
            r = 550
            p1 = (rrs[1] + rrs[3])
            p2 = rrs[6] + 5 * (((rrs[8]) ** 2)) / (rrs[3])
            x = np.log10(p1 / p2)
            ar = aw[6] + np.power(10, (-1.146 - 1.366 * x - 0.469 * (x ** 2)))  # step 2
            bbpr = ((u[6] * ar) / (1 - u[6])) - bw[6]  # step3
        else:
            r = 670
            p1 = Rrs[8] / (Rrs[1] + Rrs[3])
            p2 = 0.39 * (p1 ** 1.14)
            ar = (aw[8]) + p2  # step2
            bbpr = (u[8] * ar / (1 - (u[8])) - (bw[8]))  # step3
        eta = 2 * (1 - 1.2 * np.exp(-0.9 * (rrs[1] / rrs[6])))  # step4

        zeta = 0.74 + 0.2 / (0.8 + rrs[1] / rrs[6])  # step 7&8
        S = 0.015 + 0.002 / (0.6 + rrs[1] / rrs[6])
        xi = np.exp(S * (442.5 - 415.5))
        for i in range(nbands):
            bbp[i] = bbpr * np.power(r / Lambda[i], eta)  # step5
        a = ((1 - u) * (bw + bbp)) / u  # step6
        for i in range(nbands):
            ag443 = ((a[0] - zeta * a[1]) / (xi - zeta)) - ((aw[0] - zeta * aw[1]) / (xi - zeta))
            adg[i] = ag443 * np.exp(-S * (Lambda[i] - 443))
            aph[i] = a[i] - adg[i] - aw[i]
        return bbp, a
