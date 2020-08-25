# coding: UTF-8

import pgutil
import numpy as np

pg = pgutil.GLground()

# 描画される点のサイズを大きくする
pg.set_point_size(3)

# -1から1のランダムな数で頂点の座標を表す配列を作成
verts = np.random.rand(128, 3) * 2 - 1

while pg.run():
    # 頂点の座標に点を描画
    pg.draw(verts)
    # 頂点同士を線で接続
    pg.draw(verts, pgutil.GL_LINE_STRIP, (1., 0., 0.,))
