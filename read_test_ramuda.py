import numpy as np
import matplotlib.pyplot as plt
import pgutil
from datetime import datetime
import os

param = "dc_dev"

filename = "%s/bin/test_ramuda.bin" % param
fp = open(filename, "rb")
cc = np.fromfile(fp, dtype=np.float64).reshape(-1, 128, 128)
pg = pgutil.Playground()
datmax = cc.max()
datmin = cc.min()
color = (cc - datmin) / (datmax - datmin)
print(cc.max(), cc.min())
print(color.shape)
thres = 0.0
now = datetime.now().strftime('%s') #e.g.'1600835026'
path = '%s/output/test_ramuda/%s' % (param, now)
os.makedirs(path, exist_ok=True)
while pg.run():
    ind = (pg.totaltick) % cc.shape[0]

    pg.transform_blit_cmap(cc[ind], datmin, datmax, cmap="Reds")

    #    Max              Min                                             Background
    # r = (255 * color[ind] + 255 * (1.0 - color[ind])) * (cc[ind] > thres) + 255 * (1 - (cc[ind] > thres))
    # g = (0 * color[ind] + 255 * (1.0 - color[ind])) * (cc[ind] > thres) + 255 * (1 - (cc[ind] > thres))
    # b = (0 * color[ind] + 255 * (1.0 - color[ind])) * (cc[ind] > thres) + 255 * (1 - (cc[ind] > thres))
    # drawing = np.dstack((r, g, b))     # dstack:2次元配列→3次元配列
    # pg.transform_blit_3d(drawing)
    nn = 100 * ind
    if ind % 1 == 0:
        pg.take_screenshot("%s/ramuda_%d.png" % (path, nn))
fp.close()
