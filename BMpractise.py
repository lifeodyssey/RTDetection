from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
# setup Lambert Conformal basemap.
m = Basemap(projection='mill',
            llcrnrlat = 32,
            llcrnrlon = 130,
            urcrnrlat = 33.4,
            urcrnrlon = 130.8,
            resolution='f')
# draw coastlines.
m.drawcoastlines()
# draw a boundary around the map, fill the background.
# this background will end up being the ocean color, since
# the continents will be drawn on top.
m.drawmapboundary()
# fill continents, set lake color same as ocean color.
m.fillcontinents()
m.drawcountries()
m.drawrivers(color='#0000ff')
parallels=np.arange(32,33.4,0.2)
m.drawparallels(parallels,fontsize=10,linewidth=0.25,dashes=[7,15],
                color='k',labels=[1,0,1,1])
##dash：http://bbs.06climate.com/home.php?mod=space&uid=8934&do=blog&id=3284 参考
##linewidth：线的宽度
meridians = np.arange(130,130.8,0.2)
m.drawmeridians(meridians,fontsize=10,dashes=[7,15],
                linewidth=0.3, color='k',labels=[1,1,0,1])
plt.show()