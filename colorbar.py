import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime

now = datetime.now().strftime('%s')
path = "colorbar/%s" %now
os.makedirs(path, exist_ok=True)

colorbarMap = {
    "c":"binary",
    "D":"Blues",
    "delT":"Purples",
    "I":"binary",
    "sig":"Greens",
    "T":"Reds",
    "V":"Oranges",
    "W":"binary",
    "yVa":"gray",
}

for element, cmap in colorbarMap.items():

    fp = open("dc_dev/bin/test_%s.bin" % element, "rb")
    cc = np.fromfile(fp, dtype=np.float64).reshape(-1, 128, 128)

    maximum = cc.max()
    minimum = cc.min()
    print("%s: max=%s, min=%s" %(element, maximum, minimum))

    fig = plt.figure()
    plt.imshow(cc[0], cmap=cmap)
    plt.xticks([])
    plt.yticks([])
    plt.clim(minimum, maximum)
    plt.colorbar(ticks=[])

    fig.savefig("%s/%s.png" %(path, element))

    fp.close()
