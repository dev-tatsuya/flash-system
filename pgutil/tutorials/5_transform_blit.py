# coding: UTF-8

import pgutil
import numpy as np

pg = pgutil.Playground()

# x, y座標を表す2次元配列を生成する
y, x = np.meshgrid(np.linspace(0, 4 * np.pi, 128), np.linspace(0, 4 * np.pi, 128))

while pg.run():
    # 描画用2次元配列の生成
    z = np.cos(x + pg.totaltick * 0.1) * np.sin(y + pg.totaltick * 0.2)

    # 2次元配列を描画
    pg.transform_blit(z)