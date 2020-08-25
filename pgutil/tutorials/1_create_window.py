# coding: UTF-8

import pgutil

# Playgroundインスタンスの作成
pg = pgutil.Playground()

# メインループ
# 1フレームに1回実行される
# フレームレートは60fpsに設定されているため1秒間に60回以上実行されることはない
while pg.run():
    # 何もしない
    pass
