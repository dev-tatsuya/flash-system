import numpy as np
import matplotlib.pyplot as plt
import pgutil

fp = open("dc_dev/bin/test_yVa.bin", "rb")
c = np.fromfile(fp, dtype=np.float64)
cc = c.reshape(-1, 128, 128)
pg = pgutil.Playground()
datmax = cc.max()
datmin = cc.min()
color = (cc - datmin) / (datmax - datmin)
print(cc.max(), cc.min())
print(color.shape)
thres = 0

now = datetime.now().strftime('%s') #e.g.'1600835026'
path = 'dc_dev/output/test_yVa/%s' % now
os.makedirs(path, exist_ok=True)

while pg.run():
    ind = (pg.totaltick) % cc.shape[0]

    pg.transform_blit_cmap(cc[ind], 0.015000000109873806, 0.01500000102695368, cmap="gray")

    nn = 2000 * ind
    if ind % 1 == 0:
        pg.take_screenshot("%s/yVa_%d.png" % (path, nn))
fp.close()
