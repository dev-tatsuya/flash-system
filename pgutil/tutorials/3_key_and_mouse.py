# coding: UTF-8

import pgutil

pg = pgutil.Playground()


# デコレータを用いた特定のキーへの機能の割り当て
@pg.keyfunc("a")
def key_pressed0():
    print("A is pressed.")


def key_pressed1():
    print("B is pressed.")


# 登録用の関数を用いた特定のキーへの機能の割り当て
pg.set_key_func("b", key_pressed1)


# デコレータを用いたマウスのボタンへの機能の割り当て
# 登録した関数は引数としてマウスがクリックされた座標を受け取れるようになっている必要がある.
@pg.mousedownfunc(1)
def mouse_clicked(pos):
    print(f"Mouse clicked at {pos}")

@pg.screenmodefunc()
def screen_changed():
    print(f"Current screen size is {pg.size}")

while pg.run():
    pass
