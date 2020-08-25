# coding: UTF-8

import pgutil
import numpy as np

# 空間分割数の指定
n = 32

pg = pgutil.GLground()

# x, y, z座標を表す3次元配列を生成する
y, x, z = np.meshgrid(np.linspace(0, 2 * np.pi, n), np.linspace(0, 2 * np.pi, n), np.linspace(0, 2 * np.pi, n))


while pg.run():
    # 描画用3次元配列の生成
    w = np.cos(x + pg.totaltick * 0.01) * np.cos(y + pg.totaltick * 0.02) * np.cos(z + pg.totaltick * 0.03)

    # 0.5より大きい部分を不透明としてマーチングキューブ法により描画
    pg.draw_marching_cubes(w, 0.5)
    # -0.5より小さい部分を不透明としてボクセルデータとして描画
    pg.draw_voxel(w < -0.5)
    # 断面がマウスで動かせるヒートマップとして描画
    pg.draw_tex_cube(w, mouse=True)
