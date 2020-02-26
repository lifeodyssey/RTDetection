import pyresample
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib import colors
import math
from scipy.interpolate import RectSphereBivariateSpline

def swath_resampling(src_data: np.ma.array, src_lon: np.array, src_lat: np.array,
                     trg_lon: np.array, trg_lat: np.array, search_radius: float):
    if len(trg_lon.shape) == 1:
        grid_def = pyresample.geometry.SwathDefinition(*np.meshgrid(trg_lon, trg_lat))
    else:
        grid_def = pyresample.geometry.SwathDefinition(lons=trg_lon, lats=trg_lat)

    # source grid with original swath data
    swath_def = pyresample.geometry.SwathDefinition(lons=src_lon, lats=src_lat)

    # resample (here we use nearest. Bilinear, gaussian and custom defined methods are available)
    # for more, visit https://pyresample.readthedocs.io/en/latest/
    result = pyresample.kd_tree.resample_nearest(swath_def, src_data, grid_def, epsilon=0.5,
                                                 fill_value=None, radius_of_influence=search_radius)
    return result, grid_def


def plot_geo_image(sds: np.ma.array, lon: np.ndarray, lat: np.ndarray, log10: bool = True, title: str = None,
                   label: str = None,
                   caxis: list = None, lon_range: list = None, lat_range: list = None, save_image: str = None,
                   dpi: int = 100):
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
                    resolution='i', lon_0=lon_0, lat_0=lat_0, projection='tmerc')
    else:
        m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(),
                    urcrnrlon=lon.max(), urcrnrlat=lat.max(),
                    resolution='i', lon_0=lon_0, lat_0=lat_0, projection='tmerc')
    x2d, y2d = m(lon, lat)

    fig = plt.figure(figsize=(8, 8 * m.aspect))
    ax = fig.add_axes([0.08, 0.1, 0.7, 0.7], facecolor='white')
    # changed to facecolor 8 October 2019

    if (lon_range is not None) and (lat_range is not None):
        parallels = np.linspace(min(lat_range), max(lat_range), 4)
        meridians = np.linspace(min(lon_range), max(lon_range), 3)
    else:
        parallels = meridians = None

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

    m.drawcoastlines()
    m.drawcountries()
    plt.show()

    if save_image is not None:
        plt.savefig(save_image, dpi=dpi, facecolor='w', edgecolor='w', orientation='portrait')
        plt.show()
        plt.close()
# code as is provided by the JAXA SGLI Team
def creategrid(min_lon, max_lon, min_lat, max_lat, cell_size_deg, mesh=False):
#Output grid within geobounds and specifice cell size
#cell_size_deg should be in decimal degrees’’’

    min_lon = math.floor(min_lon)
    max_lon = math.ceil(max_lon)
    min_lat = math.floor(min_lat)
    max_lat = math.ceil(max_lat)
    lon_num = (max_lon - min_lon)/cell_size_deg
    lat_num = (max_lat - min_lat)/cell_size_deg
    grid_lons = np.zeros(lon_num) # fill with lon_min
    grid_lats = np.zeros(lat_num) # fill with lon_max
    grid_lons = grid_lons + (np.assary(range(lon_num))*cell_size_deg)
    grid_lats = grid_lats + (np.assary(range(lat_num))*cell_size_deg)
    grid_lons, grid_lats = np.meshgrid(grid_lons, grid_lats)
    grid_lons = np.ravel(grid_lons)
    grid_lats = np.ravel(grid_lats)
    #if mesh = True:
    # grid_lons = grid_lons
    # grid_lats = grid_lats
    return grid_lons, grid_lats
def geointerp(lats,lons,data,grid_size_deg, mesh=False):
#We want to interpolate it to a global x-degree grid’’’
    deg2rad = np.pi/180.
    new_lats = np.linspace(grid_size_deg, 180, 180/grid_size_deg)
    new_lons = np.linspace(grid_size_deg, 360, 360/grid_size_deg)
    new_lats_mesh, new_lons_mesh = np.meshgrid(new_lats*deg2rad, new_lons*deg2rad)
    #We need to set up the interpolator object’’’
    lut = RectSphereBivariateSpline(lats*deg2rad, lons*deg2rad, data)
    #Finally we interpolate the data. The RectSphereBivariateSpline
    #object only takes 1-D arrays as input, therefore we need to do some reshaping.’’’
    new_lats = new_lats_mesh.ravel()
    new_lons = new_lons_mesh.ravel()
    data_interp = lut.ev(new_lats,new_lons)
    if mesh == True:
        data_interp = data_interp.reshape((360 / grid_size_deg,
                                           180 / grid_size_deg)).T
        return new_lats / deg2rad, new_lons / deg2rad, data_interp