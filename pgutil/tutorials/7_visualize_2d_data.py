# coding: UTF-8

import pgutil
import numpy as np

pg = pgutil.GLground()

# 描画される線の太さを太くする
pg.set_line_width(3)

# x, y座標を表す2次元配列を生成する
y, x = np.meshgrid(np.linspace(0, 4 * np.pi, 128), np.linspace(0, 4 * np.pi, 128))

while pg.run():
    # 描画用2次元配列の生成
    z = np.cos(x + pg.totaltick * 0.01) * np.sin(y + pg.totaltick * 0.02)

    # 2次元配列を色付きのメッシュで描画
    pg.draw_mesh_colored(z)
    # 2次元配列を陰影処理のなされた平面で描画
    pg.draw_surface(z)
