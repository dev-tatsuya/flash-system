# coding: UTF-8

import pygame.draw as pyd
import pygame.gfxdraw as gfx


class DrawWrapper:
    """pygame.drawモジュール内の関数の呼び出し用の薄いラッパークラス.

    このクラスはpygame.drawモジュール以下に含まれる基本図形の描画用の関数の呼び出しの補助するためののクラスです.
    このクラスに含まれる関数は殆どの場合、pygame.drawモジュール内の関数のSurfaceオブジェクトの指定を省略できるようにしたものです.
    このため、それぞれの関数の詳しい説明は http://westplain.sakuraweb.com/translate/pygame/Draw.cgi を参照してください.
    このクラス内の関数はPlaygroundのインスタンスを生成後、Playground.draw以下の名前空間から参照することができます.

    """

    def __init__(self, pg):
        self.pg = pg

    def aaline(self, color, startpos, endpos, blend=1):
        """アンチエイリアス処理のなされた直線を描画します.

        Parameters
        ----------
        color : tuple of int
            描画に使用される色を指定します.
        startpos : tuple of float
            描画される直線の始点を指定します.
        endpos : tuple of float
            描画される直線の終点を指定します.
        blend : int or bool
            Trueに設定された場合、描画位置のピクセルに上書きはされず , 半透明の状態で描画されます.

        Returns
        -------
        rect : pygame.Rect
            描画によって影響を受けた範囲を示すpygame.Rectオブジェクトを返します.

        """

        return pyd.aaline(self.pg.screen, color, startpos, endpos, blend)

    def aalines(self, color, closed, pointlist, blend=1):
        """アンチエイリアス処理のなされた複数の連結された直線を描画します.

        Parameters
        ----------
        color : tuple of int
            描画に使用される色を指定します.
        closed : bool
            Trueに設定された場合, 始点と終点を直線でつないで描画します.
        pointlist : array-like
            頂点の座標を格納した2次元配列.
            配列の2次元目の要素数は2で, 少なくとも2点以上の座標を設定する必要があります.
        blend : int or bool
            Trueに設定された場合、描画位置のピクセルに上書きはされず , 半透明の状態で描画されます.

        Returns
        -------
        rect : pygame.Rect
            描画によって影響を受けた範囲を示すpygame.Rectオブジェクトを返します.

        """

        return pyd.aalines(self.pg.screen, color, closed, pointlist, blend)

    def arc(self, color, Rect, start_angle, stop_angle, width=1):
        """楕円の弧を描画します.

        Parameters
        ----------
        color : tuple of int
            描画に使用される色を指定します.
        Rect : pygame.Rect or tuple of float
            楕円に外接する長方形を示すpygame.Rectオブジェクトもしくは座標を格納したタプル.
            タプルは(x_min, y_min, width, height)の順番で値を格納する.
        start_angle : float
            弧の始点に対応する角度.
            単位はラジアン.
        stop_angle : float
            弧の終点に対応する角度.
            単位はラジアン.
        width : float
            描画される弧の太さ.

        Returns
        -------
        rect : pygame.Rect
            描画によって影響を受けた範囲を示すpygame.Rectオブジェクトを返します.

        """

        return pyd.arc(self.pg.screen, color, Rect, start_angle, stop_angle, width)

    def circle(self, color, pos, radius, width=0):
        """円を描画します.

        Parameters
        ----------
        color : tuple of int
            描画に使用される色を指定します.
        pos : tuple of float
            円の中心の座標を格納したタプル.
        radius : float
            描画される円の半径.
        width : float
            描画される円の太さ.
            0を指定した場合は内側を塗りつぶして描画する.

        Returns
        -------
        rect : pygame.Rect
            描画によって影響を受けた範囲を示すpygame.Rectオブジェクトを返します.

        """

        return pyd.circle(self.pg.screen, color, pos, radius, width)

    def ellipse(self, color, Rect, width=0):
        """楕円を描画します.

        Parameters
        ----------
        color : tuple of int
            描画に使用される色を指定します.
        Rect : pygame.Rect or tuple of float
            楕円に外接する長方形を示すpygame.Rectオブジェクトもしくは座標を格納したタプル.
            タプルは(x_min, y_min, width, height)の順番で値を格納する.
        width : float
            描画される円の太さ.
            0を指定した場合は内側を塗りつぶして描画する.

        Returns
        -------
        rect : pygame.Rect
            描画によって影響を受けた範囲を示すpygame.Rectオブジェクトを返します.

        """

        return pyd.ellipse(self.pg.screen, color, Rect, width)

    def line(self, color, start_pos, end_pos, width=1):
        """直線を描画します.

        Parameters
        ----------
        color : tuple of int
            描画に使用される色を指定します.
        startpos : tuple of float
            描画される直線の始点を指定します.
        endpos : tuple of float
            描画される直線の終点を指定します.
        width : float
            描画される直線の太さ.

        Returns
        -------
        rect : pygame.Rect
            描画によって影響を受けた範囲を示すpygame.Rectオブジェクトを返します.

        """

        return pyd.line(self.pg.screen, color, start_pos, end_pos, width)

    def lines(self, color, closed, pointlist, width=1):
        """複数の連結された直線を描画します.

        Parameters
        ----------
        color : tuple of int
            描画に使用される色を指定します.
        closed : bool
            Trueに設定された場合, 始点と終点を直線でつないで描画します.
        pointlist : array-like
            頂点の座標を格納した2次元配列.
            配列の2次元目の要素数は2で, 少なくとも2点以上の座標を設定する必要があります.
        width : float
            描画される直線の太さ.

        Returns
        -------
        rect : pygame.Rect
            描画によって影響を受けた範囲を示すpygame.Rectオブジェクトを返します.

        """

        return pyd.lines(self.pg.screen, color, closed, pointlist, width)

    def polygon(self, color, pointlist, width=0):
        """複数の頂点からなる多角形を描画します.

        Parameters
        ----------
        color : tuple of int
            描画に使用される色を指定します.
        pointlist : array-like
            頂点の座標を格納した2次元配列.
            配列の2次元目の要素数は2で, 少なくとも3点以上の座標を設定する必要があります.
        width : float
            描画される線の太さ.
            0を指定した場合は内側を塗りつぶして描画する.

        Returns
        -------
        rect : pygame.Rect
            描画によって影響を受けた範囲を示すpygame.Rectオブジェクトを返します.

        """

        return pyd.polygon(self.pg.screen, color, pointlist, width)

    def rect(self, color, Rect, width=0):
        """長方形を描画します.

        Parameters
        ----------
        color : tuple of int
            描画に使用される色を指定します.
        Rect : pygame.Rect or tuple of float
            長方形を示すpygame.Rectオブジェクトもしくは座標を格納したタプル.
            タプルは(x_min, y_min, width, height)の順番で値を格納する.
        width : float
            描画される円の太さ.
            0を指定した場合は内側を塗りつぶして描画する.

        Returns
        -------
        rect : pygame.Rect
            描画によって影響を受けた範囲を示すpygame.Rectオブジェクトを返します.

        """

        return pyd.rect(self.pg.screen, color, Rect, width)


class GfxWrapper:
    """pygame.gfxdrawモジュール内の関数の呼び出し用の薄いラッパークラス.

    このクラスはpygame.gfxdrawモジュール以下に含まれる基本図形の描画用の関数の呼び出しの補助するためののクラスです.
    このクラスに含まれる関数は殆どの場合、pygame.gfxdrawモジュール内の関数のSurfaceオブジェクトの指定を省略できるようにしたものです.
    このため、それぞれの関数の詳しい説明は http://westplain.sakuraweb.com/translate/pygame/Gfxdraw.cgi を参照してください.
    このクラス内の関数はPlaygroundのインスタンスを生成後、Playground.gfxdraw以下の名前空間から参照することができます.

    """

    def __init__(self, pg):
        self.pg = pg

    def aacircle(self, x, y, r, color):
        """アンチエイリアス処理のなされた円を描画します.

        Parameters
        ----------
        x : int
            円の中心のx座標.
        y : int
            円の中心のy座標.
        r : int
            描画される円の半径.
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.aacircle(self.pg.screen, x, y, r, color)

    def aaellipse(self, x, y, rx, ry, color):
        """アンチエイリアス処理のなされた楕円を描画します.

        Parameters
        ----------
        x : int
            楕円の左上のx座標.
        y : int
            楕円の左上のy座標.
        rx : int
            描画される楕円の幅
        ry : int
            描画される楕円の高さ
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.aaellipse(self.pg.screen, x, y, rx, ry, color)

    def aapolygon(self, points, color):
        """複数の頂点からなるアンチエイリアス処理のなされた多角形を描画します.

        Parameters
        ----------
        points : array-like
            頂点の座標を格納した2次元配列.
            配列の2次元目の要素数は2で, 少なくとも3点以上の座標を設定する必要があります.
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.aapolygon(self.pg.screen, points, color)

    def aatrigon(self, x1, y1, x2, y2, x3, y3, color):
        """アンチエイリアス処理のなされた三角形を描画します.

        Parameters
        ----------
        x1 : int
            三角形の1つめの頂点のx座標.
        y1 : int
            三角形の1つめの頂点のy座標.
        x2 : int
            三角形の2つめの頂点のx座標.
        y2 : int
            三角形の2つめの頂点のy座標.
        x3 : int
            三角形の3つめの頂点のx座標.
        y3 : int
            三角形の3つめの頂点のy座標.
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.aatrigon(self.pg.screen, x1, y1, x2, y2, x3, y3, color)

    def arc(self, x, y, r, start, end, color):
        """円の弧を描画します.

        Parameters
        ----------
        x : int
            円の中心のx座標.
        y : int
            円の中心のy座標.
        r : int
            円の半径.
        start : float
            弧の始点に対応する角度.
            単位はラジアン.
        end : float
            弧の終点に対応する角度.
            単位はラジアン.
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.arc(self.pg.screen, x, y, r, start, end, color)

    def bezier(self, points, steps, color):
        """ベジェ曲線を描画します.

        Parameters
        ----------
        pointlist : array-like
            ベジェ曲線の制御点の座標を格納した2次元配列.
            配列の2次元目の要素数は2で, 少なくとも2点以上の座標を設定する必要があります.
        steps : tuple of int
            ベジェ曲線の次数を指定します.
        closed : bool
            Trueに設定された場合, 始点と終点を直線でつないで描画します.

        """

        return gfx.bezier(self.pg.screen, points, steps, color)

    def box(self, *args, **kwargs):
        """長方形を描画します.

        """

        return gfx.box(self.pg.screen, *args, **kwargs)

    def circle(self, x, y, r, color):
        """円を描画します.

        Parameters
        ----------
        x : int
            円の中心のx座標.
        y : int
            円の中心のy座標.
        r : int
            描画される円の半径.
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.circle(self.pg.screen, x, y, r, color)

    def ellipse(self, x, y, rx, ry, color):
        """楕円を描画します.

        Parameters
        ----------
        x : int
            楕円の左上のx座標.
        y : int
            楕円の左上のy座標.
        rx : int
            描画される楕円の幅
        ry : int
            描画される楕円の高さ
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.ellipse(self.pg.screen, x, y, rx, ry, color)

    def filled_circle(self, x, y, r, color):
        """内部が塗りつぶされた円を描画します.

        Parameters
        ----------
        x : int
            円の中心のx座標.
        y : int
            円の中心のy座標.
        r : int
            描画される円の半径.
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.filled_circle(self.pg.screen, x, y, r, color)

    def filled_ellipse(self, x, y, rx, ry, color):
        """内部が塗りつぶされた楕円を描画します.

        Parameters
        ----------
        x : int
            楕円の左上のx座標.
        y : int
            楕円の左上のy座標.
        rx : int
            描画される楕円の幅
        ry : int
            描画される楕円の高さ
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.filled_ellipse(self.pg.screen, x, y, rx, ry, color)

    def filled_polygon(self, points, color):
        """複数の頂点からなる内部が塗りつぶされた多角形を描画します.

        Parameters
        ----------
        points : array-like
            頂点の座標を格納した2次元配列.
            配列の2次元目の要素数は2で, 少なくとも3点以上の座標を設定する必要があります.
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.filled_polygon(self.pg.screen, points, color)

    def filled_trigon(self, x1, y1, x2, y2, x3, y3, color):
        """内部が塗りつぶされた三角形を描画します.

        Parameters
        ----------
        x1 : int
            三角形の1つめの頂点のx座標.
        y1 : int
            三角形の1つめの頂点のy座標.
        x2 : int
            三角形の2つめの頂点のx座標.
        y2 : int
            三角形の2つめの頂点のy座標.
        x3 : int
            三角形の3つめの頂点のx座標.
        y3 : int
            三角形の3つめの頂点のy座標.
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.filled_trigon(self.pg.screen, x1, y1, x2, y2, x3, y3, color)

    def hline(self, x1, x2, y, color):
        """x軸に平行な直線を描画します.

        Parameters
        ----------
        x1 : int
            直線の始点のx座標.
        x2 : int
            直線の終点のx座標.
        y : int
            直線のy座標.
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.hline(self.pg.screen, x1, x2, y, color)

    def line(self, x1, y1, x2, y2, color):
        """直線を描画します.

        Parameters
        ----------
        x1 : int
            直線の始点のx座標.
        y1 : int
            直線の始点のy座標.
        x2 : int
            直線の終点のx座標.
        y2 : int
            直線の終点のy座標.
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.line(self.pg.screen, x1, y1, x2, y2, color)

    def pie(self, x, y, r, start, end, color):
        """内部の塗りつぶされた円弧を描画します.

        Parameters
        ----------
        x : int
            円の中心のx座標.
        y : int
            円の中心のy座標.
        r : int
            円の半径.
        start : float
            弧の始点に対応する角度.
            単位はラジアン.
        end : float
            弧の終点に対応する角度.
            単位はラジアン.
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.pie(self.pg.screen, x, y, r, start, end, color)

    def pixel(self, x, y, color):
        """1ピクセルの点を描画します.

        Parameters
        ----------
        x : int
            描画されるピクセルのx座標.
        y : int
            描画されるピクセルのy座標.
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.pixel(self.pg.screen, x, y, color)

    def polygon(self, points, color):
        """複数の頂点からなる多角形を描画します.

        Parameters
        ----------
        points : array-like
            頂点の座標を格納した2次元配列.
            配列の2次元目の要素数は2で, 少なくとも3点以上の座標を設定する必要があります.
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.polygon(self.pg.screen, points, color)

    def rectangle(self, rect, color):
        """長方形を描画します.

        Parameters
        ----------
        rect : pygame.Rect or tuple of float
            長方形を示すpygame.Rectオブジェクトもしくは座標を格納したタプル.
            タプルは(x_min, y_min, width, height)の順番で値を格納する.
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.rectangle(self.pg.screen, rect, color)

    def textured_polygon(self, points, texture, tx, ty):
        """テクスチャ処理のなされた複数の頂点からなる多角形を描画します.

        """
        return gfx.textured_polygon(self.pg.screen, points, texture, tx, ty)

    def trigon(self, x1, y1, x2, y2, x3, y3, color):
        """三角形を描画します.

        Parameters
        ----------
        x1 : int
            三角形の1つめの頂点のx座標.
        y1 : int
            三角形の1つめの頂点のy座標.
        x2 : int
            三角形の2つめの頂点のx座標.
        y2 : int
            三角形の2つめの頂点のy座標.
        x3 : int
            三角形の3つめの頂点のx座標.
        y3 : int
            三角形の3つめの頂点のy座標.
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.trigon(self.pg.screen, x1, y1, x2, y2, x3, y3, color)

    def vline(self, x, y1, y2, color):
        """y軸に平行な直線を描画します.

        Parameters
        ----------
        x : int
            直線のx座標.
        y1 : int
            直線の始点のy座標.
        y2 : int
            直線の終点のy座標.
        color : tuple of int
            描画に使用される色を指定します.

        """

        return gfx.vline(self.pg.screen, x, y1, y2, color)
