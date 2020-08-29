import numpy as np
import pgutil

pg = pgutil.Playground()
ph = np.tile(np.linspace(0, 1, pg.size[0]), (pg.size[1], 1))

print(ph.shape)
while pg.run():
    pg.transform_blit_cmap(ph, 0, 1, cmap="binary")