# coding: UTF-8

import pgutil

pg = pgutil.Playground()

# 毎フレーム画面を塗りつぶすように設定
pg.fill = True

while pg.run():
    # 直線の描画
    pg.draw.line((0, 255, 0), (300, 300), (400, 450), 3)
    # 長方形の描画
    pg.draw.rect((0, 0, 255), (100, 300, 100, 150))
    # マウスカーソルの位置に円を描画
    pg.draw.circle((255, 0, 0), pg.mouse, 80)
