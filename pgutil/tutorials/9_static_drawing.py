# coding: UTF-8

import pgutil
import numpy as np

# 空間分割数の指定
n = 256

pg = pgutil.GLground()

# x, y, z座標を表す3次元配列を生成する
y, x, z = np.meshgrid(np.linspace(0, 2 * np.pi, n), np.linspace(0, 2 * np.pi, n), np.linspace(0, 2 * np.pi, n))
# 描画用3次元配列の生成
w = np.cos(x) * np.cos(y) * np.cos(z)

# マーチングキューブ法のデータを準備
mc_verts, mc_norms = pg.gen_marching_data(w, 0.5)
# ボクセル描画用のデータを準備
vx_verts, vx_colors = pg.gen_voxel_data(w < -0.5)
# 3次元ヒートマップ描画用のデータを準備
pg.bind_tex_cube(w)

while pg.run():
    # マーチングキューブデータの描画
    pg.draw_with_norm(mc_verts, pgutil.GL_TRIANGLES, mc_norms)
    # ボクセルデータの描画
    pg.draw_with_color(vx_verts, pgutil.GL_QUADS, vx_colors)
    # 3次元ヒートマップデータの描画
    pg.draw_tex_cube(mouse=True)
