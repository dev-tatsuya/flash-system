# coding: UTF-8

import pygame
import pgutil._keyutil as keyutil
import pgutil._drawutil as drawutil
import sys
import numpy as np
from time import ctime
from OpenGL.GL import *
from OpenGL.GLU import *
import matplotlib.cm as cm
from skimage import measure


def _get_surf_array(data, maxinput, col):
    temp = data / float(maxinput)
    r = (temp * col[0]).astype(np.uint8)
    g = (temp * col[1]).astype(np.uint8)
    b = (temp * col[2]).astype(np.uint8)
    return np.dstack((r, g, b))


def _blit_array(screen, array, maxinput, mininput, col):
    pygame.surfarray.pixels3d(screen)[:] = _get_surf_array(array - mininput, maxinput - mininput, col)


def _blit_array_3d(screen, array):
    pygame.surfarray.pixels3d(screen)[:] = np.clip(array, 0, 255).astype(np.uint8)


def _gen_mesh_index(shape):
    x, y = shape
    index = []
    for j in range(y):
        for i in range(x):
            if j % 2 == 0:
                index.append(j * x + i)
            else:
                index.append(j * x + x - i - 1)
    if y % 2 == 0:
        for i in range(x):
            for j in range(y):
                if i % 2 == 0:
                    index.append(y * (x - 1) + i - j * x)
                else:
                    index.append(i + j * x)
    else:
        for i in range(x):
            for j in range(y):
                if i % 2 == 0:
                    index.append(x * y - 1 - i - j * x)
                else:
                    index.append(x - 1 - i + j * x)

    yf, xf = np.meshgrid(np.linspace(-1, 1, shape[0]), np.linspace(-1, 1, shape[1]))

    draw_cahe = np.vstack((xf.ravel(), np.zeros(shape[0] * shape[1]), yf.ravel()))
    return np.array(index), draw_cahe


def _get_edges(array):
    d = array.shape[0]
    w = 2 / d
    offs = w * np.array(
        [[[1, 0, 0], [1, 0, 1], [0, 0, 1]], [[0, 1, 0], [0, 1, 1], [0, 0, 1]], [[1, 0, 0], [1, 1, 0], [0, 1, 0]]])

    results = []
    x, y, z = np.meshgrid(np.linspace(-1, 1, array.shape[2] + 1, dtype=np.float32),
                          np.linspace(-1, 1, array.shape[1] + 1, dtype=np.float32),
                          np.linspace(-1, 1, array.shape[0] + 1, dtype=np.float32))
    for i in range(3):
        temp = np.zeros((array.shape[0] + 1, array.shape[1] + 1, array.shape[2] + 1), np.int8)
        if i == 0:
            temp[:-1, :-1, :-1] = array[:]
            temp[1:, :-1, :-1] -= array[:]
        elif i == 1:
            temp[:-1, :-1, :-1] = array[:]
            temp[:-1, 1:, :-1] -= array[:]
        else:
            temp[:-1, :-1, :-1] = array[:]
            temp[:-1, :-1, 1:] -= array[:]
        for j in [-1, 1]:
            mask = temp == j
            start = np.dstack([x[mask], y[mask], z[mask]])[0]
            second = start + offs[i, 0]
            third = start + offs[i, 1]
            fourth = start + offs[i, 2]
            results.append(
                np.hstack((start, second, third, fourth)).reshape((start.shape[0] * 4, 3)).astype(np.float32))
    return results


def _get_tex_array(text, font):
    # img = surfarray.array3d(font.render(text, False, (255, 255, 255))).astype(np.uint8)[:, :, 0]
    img = surfarray.array3d(_render_text(text, font)).astype(np.uint8)[:, :, 0]
    return np.dstack((img, img, img, img))


def _load_texture(text, font):
    tex_num = glGenTextures(1)
    im = _get_tex_array(text, font)
    glBindTexture(GL_TEXTURE_2D, tex_num)
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, im.shape[1], im.shape[0], 0, GL_RGBA, GL_UNSIGNED_BYTE, im)
    return tex_num, im.shape[:2]


def _render_text(text, font):
    lines = text.split("\n")
    surfs = [font.render(l, False, (255, 255, 255)) for l in lines]
    width = surfs[0].get_width()
    height = surfs[0].get_height()
    res = Surface((width, height * len(surfs)))
    for i in range(len(surfs)):
        res.blit(surfs[i], (0, i * height))
    return res


def _gen_cube(x_min, x_max, y_min, y_max, z_min, z_max):
    lims = np.array([x_min, x_max, y_min, y_max, z_min, z_max])

    points = np.zeros((8, 3))
    ind0 = np.arange(8)
    points[:, 0] = lims[ind0 // 1 % 2 + 0]
    points[:, 1] = lims[ind0 // 2 % 2 + 2]
    points[:, 2] = lims[ind0 // 4 % 2 + 4]
    inds = [0, 2, 3, 1, 4, 5, 7, 6, 0, 1, 5, 4, 2, 6, 7, 3, 0, 4, 6, 2, 1, 5, 7, 3]
    result = points[inds]
    return result, (result / 2) + 0.5


class __PgBase:
    """PlaygroundとGLgroundのベースとなるクラス.

    このクラスのメンバ関数はPlaygroundとGLgroundの両方から呼び出すことができます.

    Parameters
    ----------
    w_ : int, optional
        ウィンドウの幅,
        デフォルト値は500.
    h_ : int, optional
        ウィンドウの高さ,
        デフォルト値は500.
    name : str, optional
        ウィンドウのタイトル,
        デフォルト値は"PgBase".

    Attributes
    ----------
    size : tuple of int
        現在のウィンドウのサイズ.
    col : tuple of int
        いくつかの描画関数で使用されるデフォルトの色.
    back_col : tuple of int
        現在の背景色, PlaygroundとGLgroundの両方で有効.
    screen : pygame.Surface
        ウィンドウのpygame.Surfaceオブジェクト. pygameの機能を直接呼び出すときに使用する.
    ticking : bool
        spaceキーでTrue/Falseが切替可能なフラグ.
    totaltick : int
        ウィンドウ生成時からのフレーム数. tickingがFalseのときはカウントが増加しない.
    debug_flag : bool
        F3キーでTrue/Falseが切替可能なフラグ. Trueの間タイトルバーにfpsが表示される.
    mouse : tuple of int
        現在のマウスの位置を表すタプル.


    """

    def __init__(self, w_=500, h_=500, name="PgBase"):

        self.initial_size = (w_, h_)
        self.size = self.initial_size
        self.center = (self.initial_size[0] // 2, self.initial_size[1] // 2)

        self.col = (0, 255, 255)
        self.back_col = (0, 0, 0)

        pygame.init()

        self.fullscreen_flag = False
        self.display_size = (pygame.display.Info().current_w, pygame.display.Info().current_h)
        self.screen = self._get_screen(self.fullscreen_flag)

        self.title = name
        pygame.display.set_caption(self.title)

        self.clock = pygame.time.Clock()

        self.totaltick = -1
        self.ticking = True
        self.debug_flag = False
        self.running = True
        self.exit_flag = True

        self.mouse = (0,0)
        self.mouseX, self.mouseY = self.mouse

        self.mousedownfuncs = {}
        self.mouseupfuncs = {}

        self.keyfuncs = {}

        self.screenmodefuncs = []

        self.debug_list = []

        self.cmap_cache = {}

        self.screenshot_flag = False

        self.font = font.Font(font.match_font("arial"), 200)

    def _update(self):
        """Internal function for updating pygame window.

        Needs to be overridden.
        """
        pass

    def _process_event(self):
        """Internal function for processing pygame events.
        """
        for event in pygame.event.get():
            if event.type == QUIT or (event.type == KEYDOWN and event.key == K_ESCAPE):
                self.running = False
            elif event.type == KEYDOWN and event.key in self.keyfuncs:
                self.keyfuncs[event.key]()
            elif event.type == MOUSEBUTTONDOWN and event.button in self.mousedownfuncs:
                self.mousedownfuncs[event.button](event.pos)
            elif event.type == MOUSEBUTTONUP and event.button in self.mouseupfuncs:
                self.mouseupfuncs[event.button](event.pos)
            elif event.type == KEYDOWN and event.key == K_SPACE:
                self.ticking = not self.ticking
            elif event.type == KEYDOWN and event.key == K_F2:
                self.screenshot_flag = True
            elif event.type == KEYDOWN and event.key == K_F3:
                self.debug_flag = not self.debug_flag
                if self.debug_flag:
                    pygame.display.set_caption(self.title + "  fps : " + str(round(self.get_fps() * 10) / 10))
                else:
                    pygame.display.set_caption(self.title)
            elif event.type == KEYDOWN and event.key == K_F11:
                self.fullscreen_flag = not self.fullscreen_flag
                if self.fullscreen_flag:
                    self.size = self.display_size
                    self.screen = self._get_screen(self.fullscreen_flag)
                else:
                    self.size = self.initial_size
                    self.screen = self._get_screen(self.fullscreen_flag)
                self.center = (self.size[0] // 2, self.size[1] // 2)
                for func in self.screenmodefuncs:
                    func()
            elif event.type == MOUSEMOTION:
                self.mouse = event.pos
                self.mouseX, self.mouseY = self.mouse

    def _get_screen(self, isfull=False):
        """Internal function for generating pygame screen surface.

        Needs to be overridden.

        Parameters
        ----------
        isfull : bool, optional
            If set to True, the screen is initialized as full screen surface.
            Default value is False.
        """
        return None

    def _get_cmap(self, cmap):
        if not cmap in self.cmap_cache:
            if cmap == "AVS":
                self.cmap_cache[cmap] = np.array(
                    [cm.get_cmap("gist_rainbow")(i / 255 * 16 / 21) for i in range(256)]).astype(np.float64)[:,
                                        :-1]
            elif cmap == "AVS_r":
                self.cmap_cache[cmap] = np.array(
                    [cm.get_cmap("gist_rainbow")(i / 255 * 16 / 21) for i in range(256)]).astype(np.float64)[:,
                                        :-1][::-1]
            else:
                self.cmap_cache[cmap] = np.array([cm.get_cmap(cmap)(i / 255) for i in range(256)]).astype(np.float64)[:,
                                        :-1]

        return self.cmap_cache[cmap]

    def get_fps(self):
        """現在のfps(frames per second)を返します.

        Returns
        -------
        fps : float
            現在のfps.

        """
        return self.clock.get_fps()

    def set_mouse_down_func(self, button, func):
        """マウスのボタンが押されたときに呼び出される関数を設定します.

        Parameters
        ----------
        button : int
            マウスのボタン番号
            左クリックが1, 右クリックが3. 他のボタンも使用可能.
        func : function
            マウスのボタンが押されたときに呼び出される関数.

        """
        self.mousedownfuncs[button] = func

    def set_mouse_up_func(self, button, func):
        """マウスのボタンが離されたときに呼び出される関数を設定します.

        Parameters
        ----------
        button : int
            マウスのボタン番号.
            左クリックが1, 右クリックが3. 他のボタンも使用可能.
        func : function
            マウスのボタンが離されたときに呼び出される関数.

        """
        self.mouseupfuncs[button] = func

    def set_key_func(self, key, func):
        """特定のキーが押されたときに呼び出される関数を設定します.

        Parameters
        ----------
        button : int or str
            キー番号もしくはキーの名前. キー番号はpgutil.K_F1などのようにして取得可能. キーの名前は"F1"などとして設定可能。
        func : function
            キーが押されたときに呼び出される関数.

        """
        self.keyfuncs[keyutil.get_key(key)] = func

    def mousedownfunc(self, button):
        """マウスのボタンが押されたときに呼び出される関数を設定するデコレータです.

        Parameters
        ----------
        button : int
            マウスのボタン番号
            左クリックが1, 右クリックが3. 他のボタンも使用可能.

        """

        def register(func):
            self.mousedownfuncs[button] = func
            return func

        return register

    def mouseupfunc(self, button):
        """マウスのボタンが離されたときに呼び出される関数を設定するデコレータです.

        Parameters
        ----------
        button : int
            マウスのボタン番号
            左クリックが1, 右クリックが3. 他のボタンも使用可能.

        """

        def register(func):
            self.mouseupfuncs[button] = func
            return func

        return register

    def keyfunc(self, key):
        """特定のキーが押されたときに呼び出される関数を設定するデコレータです.

        Parameters
        ----------
        button : int or str
            キー番号もしくはキーの名前. キー番号はpgutil.K_F1などのようにして取得可能. キーの名前は"F1"などとして設定可能。
        func : function
            キーが押されたときに呼び出される関数.

        """

        def register(func):
            self.keyfuncs[keyutil.get_key(key)] = func
            return func

        return register

    def set_screen_mode_func(self, func):
        """スクリーンモードが変更されたときに呼び出される関数を設定します.

        通常、スクリーンモードが切り替わるのはF11キーが押されて全画面モードと通常モードが切り替わったときです.
        対象の関数はこの関数で設定されたタイミングにも一度呼び出されます.

        Parameters
        ----------
        func : function
            スクリーンモードが変更されたときに呼び出される関数

        """
        self.screenmodefuncs.append(func)
        func()

    def screenmodefunc(self):
        """スクリーンモードが変更されたときに呼び出される関数を設定するデコレータです.

        通常、スクリーンモードが切り替わるのはF11キーが押されて全画面モードと通常モードが切り替わったときです.
        対象の関数はこの関数で設定されたタイミングにも一度呼び出されます.

        Parameters
        ----------
        func : function
            スクリーンモードが変更されたときに呼び出される関数

        """

        def register(func):
            self.screenmodefuncs.append(func)
            func()
            return func

        return register

    def take_screenshot(self, name=None):
        """PNG形式でスクリーンショットを撮影します.

        Parameters
        ----------
        name : str, optional
            保存するファイルの名前.
            指定しない場合は現在の日付と時刻が使用される.

        """
        name = name or ctime().replace(":", "_").replace(" ", "_")
        name = name.replace(".png", "").replace(".PNG", "")
        pygame.image.save(self.screen, name + ".png")

    def run(self):
        """画面の更新とpygameのイベントの処理を行います.

        Returns
        -------
        flag : bool
            ウィンドウが閉じられていない場合, Trueを返す.

        """
        self._update()
        if not self.running:
            return False
        self._process_event()
        return True


class Playground(__PgBase):
    """2次元データの可視化用のクラス.

    Parameters
    ----------
    w_ : int, optional
        ウィンドウの幅,
        デフォルト値は500.
    h_ : int, optional
        ウィンドウの高さ,
        デフォルト値は500.
    name : str, optional
        ウィンドウのタイトル,
        デフォルト値は"Playground".

    Attributes
    ----------
    fill : bool
        ウィンドウを毎フレーム塗りつぶすかどうか. 塗りつぶす場合はback_colに指定された色で塗りつぶされる.
    """

    def __init__(self, w_=500, h_=500, name="Playground"):
        super().__init__(w_, h_, name)

        self.fill = False
        self.surf_cache = {}

        self.draw = drawutil.DrawWrapper(self)
        self.gfxdraw = drawutil.GfxWrapper(self)

    def _update(self):
        """Internal function for updating pygame window.
        """
        if not self.running:
            pygame.quit()
            if self.exit_flag:
                sys.exit()
            return

        if self.screenshot_flag:
            self.take_screenshot()
            self.screenshot_flag = False
        self.clock.tick(60)
        pygame.display.update()
        if self.fill:
            self.screen.fill(self.back_col)
        if self.ticking:
            self.totaltick += 1

        if self.debug_flag:
            pygame.display.set_caption(
                self.title + "  fps : " + str(round(self.get_fps() * 10) / 10) + "".join(self.debug_list))
        self.debug_list = []

    def transform_blit(self, array, minin=None, maxin=None, smallsurf=None):
        """2次元配列を画像として表示する.

        値の小さい部分は暗く, 値の大きい部分は明るく表示される. 表示にはPlayground.colに保存された色が使用される.
        入力された配列は自動的に伸縮され, ウィンドウサイズに一致するようにして描画される.

        Parameters
        ----------
        array : array_like
            表示する配列.
            2次元配列のみ表示が可能.
        minin : float, optional
            表示の際に最小として表示される値.
            これより小さい値は黒として表示される.
            指定がない場合は入力された配列の最小値が使用される.
        maxin : float, optional
            表示の際に最大として表示される値.
            これより大きい値はPlayground.colと同じ色で表示される.
            指定がない場合は入力された配列の最大値が使用される.
        smallsurf : pygame.Surface, optional
            自動伸縮を行う際に内部的に使用されるpygame.Surfaceオブジェクト.
            指定がない場合は自動的に生成され, 2回目以降の描画のためにキャッシュされる.
            殆どの場合, 変更の必要はありません.

        """
        if minin is not None or maxin is not None:
            array = np.clip(array, minin, maxin)

        ma = maxin if maxin is not None else array.max()
        mi = minin if minin is not None else array.min()

        if smallsurf is None:
            if not array.shape in self.surf_cache:
                self.surf_cache[array.shape] = pygame.Surface(array.shape)

            _blit_array(self.surf_cache[array.shape], array, ma, mi, self.col)
            pygame.transform.scale(self.surf_cache[array.shape], self.size, self.screen)
        else:
            _blit_array(smallsurf, array, ma, mi, self.col)
            pygame.transform.scale(smallsurf, self.size, self.screen)

    def transform_blit_smooth(self, array, minin=None, maxin=None, smallsurf=None):
        """2次元配列を画像として表示する.

        値の小さい部分は暗く, 値の大きい部分は明るく表示される. 表示にはPlayground.colに保存された色が使用される.
        入力された配列は自動的に伸縮され, ウィンドウサイズに一致するようにして描画される.
        伸縮の際には見た目がなめらかになるよう自動的に補完処理がなされる.

        Parameters
        ----------
        array : array_like
            表示する配列.
            2次元配列のみ表示が可能.
        minin : float, optional
            表示の際に最小として表示される値.
            これより小さい値は黒として表示される.
            指定がない場合は入力された配列の最小値が使用される.
        maxin : float, optional
            表示の際に最大として表示される値.
            これより大きい値はPlayground.colと同じ色で表示される.
            指定がない場合は入力された配列の最大値が使用される.
        smallsurf : pygame.Surface, optional
            自動伸縮を行う際に内部的に使用されるpygame.Surfaceオブジェクト.
            指定がない場合は自動的に生成され, 2回目以降の描画のためにキャッシュされる.
            殆どの場合, 変更の必要はありません.

        """
        if minin is not None or maxin is not None:
            array = np.clip(array, minin, maxin)

        ma = maxin if maxin is not None else array.max()
        mi = minin if minin is not None else array.min()

        if smallsurf is None:
            if not array.shape in self.surf_cache:
                self.surf_cache[array.shape] = pygame.Surface(array.shape)

            _blit_array(self.surf_cache[array.shape], array, ma, mi, self.col)
            pygame.transform.smoothscale(self.surf_cache[array.shape], self.size, self.screen)
        else:
            _blit_array(smallsurf, array, ma, mi, self.col)
            pygame.transform.smoothscale(smallsurf, self.size, self.screen)

    def transform_blit_3d(self, array, smallsurf=None):
        """3次元配列を画像として表示する.

        3次元配列をRGBのピクセルデータとして解釈して表示する.
        入力された配列は自動的に伸縮され, ウィンドウサイズに一致するようにして描画される.

        Parameters
        ----------
        array : array_like
            表示する配列.
            3次元配列のみ表示が可能.
            3次元目の要素数は3で, すべての値は0から255の範囲に入っている必要がある.
        smallsurf : pygame.Surface, optional
            自動伸縮を行う際に内部的に使用されるpygame.Surfaceオブジェクト.
            指定がない場合は自動的に生成され, 2回目以降の描画のためにキャッシュされる.
            殆どの場合, 変更の必要はありません.

        """

        if smallsurf is None:
            arrayshape = (array.shape[0], array.shape[1])
            if not arrayshape in self.surf_cache:
                self.surf_cache[arrayshape] = pygame.Surface(arrayshape)

            _blit_array_3d(self.surf_cache[arrayshape], array)
            pygame.transform.scale(self.surf_cache[arrayshape], self.size, self.screen)
        else:
            _blit_array_3d(smallsurf, array)
            pygame.transform.scale(smallsurf, self.size, self.screen)

    def transform_blit_cmap(self, array, minin=None, maxin=None, cmap=None, smallsurf=None):
        """2次元配列をカラー画像として表示する.

        2次元配列をカラーマップに従って描画し、画像として表示する.
        入力された配列は自動的に伸縮され, ウィンドウサイズに一致するようにして描画される.

        Parameters
        ----------
        array : array_like
            表示する配列.
            2次元配列のみ表示が可能.
        minin : float, optional
            表示の際に最小として表示される値.
            これより小さい値は黒として表示される.
            指定がない場合は入力された配列の最小値が使用される.
        maxin : float, optional
            表示の際に最大として表示される値.
            これより大きい値はPlayground.colと同じ色で表示される.
            指定がない場合は入力された配列の最大値が使用される.
        cmap : str, optional
            描画に使用されるカラーマップ.
            matplotlibで使用可能なカラーマップ全てと"AVS", "AVS_r"が使用可能.
            指定がない場合は"vivid"が選択される.
        smallsurf : pygame.Surface, optional
            自動伸縮を行う際に内部的に使用されるpygame.Surfaceオブジェクト.
            指定がない場合は自動的に生成され, 2回目以降の描画のためにキャッシュされる.
            殆どの場合, 変更の必要はありません.

        """
        if minin is not None or maxin is not None:
            array = np.clip(array, minin, maxin)

        ma = maxin if maxin is not None else array.max()
        mi = minin if minin is not None else array.min()

        cmap_array = self._get_cmap(cmap)

        self.transform_blit_3d(255 * cmap_array[(255 * (array - mi) / (ma - mi)).astype(np.uint8)], smallsurf)

    def _get_screen(self, isfull=False):
        """Internal function for generating pygame screen surface.

        Parameters
        ----------
        isfull : bool, optional
            If set to True, the screen is initialized as full screen surface.
            Default value is False.
        """
        return pygame.display.set_mode(self.size, pygame.FULLSCREEN) if isfull else pygame.display.set_mode(self.size)

    def blit_to_screen(self, array):
        """2次元配列を使用してウィンドウのピクセルデータを操作する.

        [0-255]の値を持つ整数配列を使用してウィンドウのピクセルデータを直接設定する.
        入力された配列はウィンドウのピクセルのRGBすべてに代入され、モノクロ画像として描画される.
        配列は自動的に伸縮されないため, 配列サイズとウィンドウサイズは一致している必要がある.

        Parameters
        ----------
        array : array_like
            ウィンドウサイズ(Playground.size)と形状の等しい配列.
            2次元配列で, すべての値は0から255の範囲に入っている必要がある.

        """
        pygame.surfarray.pixels3d(self.screen)[:] = np.dstack((array, array, array))

    def blit_to_screen_3d(self, array):
        """3次元配列を使用してウィンドウのピクセルデータを操作する.

        [0-255]の値を持つ整数配列を使用してウィンドウのピクセルデータを直接設定する.
        入力された配列はウィンドウのピクセルのRGBそれぞれに代入され、カラー画像として描画される.
        配列は自動的に伸縮されないため, 配列サイズとウィンドウサイズは一致している必要がある.

        Parameters
        ----------
        array : array_like
            ウィンドウサイズ(Playground.size)と1, 2次元目の要素数が等しい配列.
            3次元配列で, すべての値は0から255の範囲に入っている必要がある.

        """
        pygame.surfarray.pixels3d(self.screen)[:] = array[:]


class GLground(__PgBase):
    """3次元データの可視化用のクラス.

    Parameters
    ----------
    w_ : int, optional
        ウィンドウの幅,
        デフォルト値は500.
    h_ : int, optional
        ウィンドウの高さ,
        デフォルト値は500.
    name : str, optional
        ウィンドウのタイトル,
        デフォルト値は"GLground".
    antialiasing : bool, optional
        アンチエイリアス処理を行うかどうか.
        デフォルト値はFalse.
    ortho : bool, optional
        正射影による描画を行うかどうか.
        デフォルト値はFalseで視点投影により描画される.

    Attributes
    ----------
    ax : float
        x軸を軸としたカメラの回転角. 矢印キーで操作可能なほか直接代入することも可能.
    ay : float
        y軸を軸としたカメラの回転角. 矢印キーで操作可能なほか直接代入することも可能.
    az : float
        z軸を軸としたカメラの回転角.　矢印キーで操作可能なほか直接代入することも可能.
    lx : float
        カメラのx座標. shuft+矢印キーで操作可能なほか直接代入することも可能.
    ly : float
        カメラのy座標. shuft+矢印キーで操作可能なほか直接代入することも可能.
    lz : float
        カメラのz座標.　shuft+矢印キーで操作可能なほか直接代入することも可能.
    cube_flag : bool
        各軸の頂点の座標が+-1の立方体を描画するかどうか. デフォルトはTrue.
    rotating : bool
        描画範囲を自動回転させるかどうか. デフォルトはFalse.
    rot_vel : float
        描画範囲の自動回転の速度.
    """

    def __init__(self, w_=500, h_=500, name="GLground", antialiasing=False, ortho=False):
        self.antialiasing = antialiasing
        self.ortho = ortho

        self.linewidth = 1.
        self.pointsize = 1.

        super().__init__(w_, h_, name)

        self.ax = 0.0
        self.ay = 0.0
        self.az = 0.0
        self.lx = 0.0
        self.ly = 0.0
        self.lz = 0.0

        self.cube_flag = True
        self.rotating = False
        self.rot_vel = 1.0

        self.mesh_cache = {}
        self.surf_cache = {}
        self.text_cache = {}

        self.voxelcol = 0.5 + 0.5 * np.eye(3)

        self.light_initialized = False

        self.cube_tex = glGenTextures(1)

        self.last_bind_func = None

        self.set_screen_mode_func(self.clear_tex_cache)
        self.set_screen_mode_func(self._reload_3d_tex)

    def _get_screen(self, isfull=False):
        """Internal function for generating pygame screen surface.

        Parameters
        ----------
        isfull : bool, optional
            If set to True, the screen is initialized as full screen surface.
            Default value is False.
        """

        if self.antialiasing:
            pygame.display.gl_set_attribute(GL_MULTISAMPLEBUFFERS, 1)
            pygame.display.gl_set_attribute(GL_MULTISAMPLESAMPLES, 4)

        scr = pygame.display.set_mode(self.size,
                                      pygame.FULLSCREEN | pygame.OPENGL | pygame.DOUBLEBUF) if isfull else pygame.display.set_mode(
            self.size, pygame.OPENGL | pygame.DOUBLEBUF)

        if self.antialiasing:
            glEnable(GL_MULTISAMPLE)
            glEnable(GL_SAMPLE_ALPHA_TO_COVERAGE)
            glSampleCoverage(0.5, GL_FALSE)

        glEnable(GL_ALPHA_TEST)
        glEnable(GL_DEPTH_TEST)
        glViewport(0, 0, *self.size)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()

        if self.ortho:
            limit = 1.5
            ratio = self.size[0] / self.size[1]
            if ratio < 1.:
                glOrtho(-limit, limit, -limit / ratio, limit / ratio, -64., 64.)
            elif ratio > 1.:
                glOrtho(-limit * ratio, limit * ratio, -limit, limit, -64., 64.)
            else:
                glOrtho(-limit, limit, -limit, limit, -64., 64.)

            glMatrixMode(GL_MODELVIEW)
            glLoadIdentity()
        else:
            gluPerspective(60.0, self.size[0] / self.size[1], 0.01, 100.0)
            glMatrixMode(GL_MODELVIEW)

        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glEnable(GL_BLEND)

        glLineWidth(self.linewidth)
        glPointSize(self.pointsize)

        return scr

    def _update(self):
        """Internal function for updating pygame window.
        """
        if not self.running:
            pygame.quit()
            if self.exit_flag:
                sys.exit()
            return

        if self.screenshot_flag:
            self.take_screenshot()
            self.screenshot_flag = False

        self.clock.tick(60)
        glFlush()
        pygame.display.flip()

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()
        gluLookAt(0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)
        glTranslated(self.lx, self.ly, self.lz)
        glRotated(self.ax, 1.0, 0.0, 0.0)
        glRotated(self.ay, 0.0, 1.0, 0.0)
        glRotated(self.az, 0.0, 0.0, 1.0)

        glClearColor(*list(map(lambda c: c / 255., self.back_col)), 0.0)
        if self.cube_flag:
            # noinspection PyTypeChecker
            self.draw_cube(1, list(map(lambda c: c / 255., self.col)))
        if self.debug_flag:
            self.draw_axis()

        if self.ticking:
            self.totaltick += 1

        if self.debug_flag:
            pygame.display.set_caption(
                self.title + "  fps : " + str(round(self.get_fps() * 10) / 10) + "".join(self.debug_list))
        self.debug_list = []

    def _process_event(self):
        """Internal function for processing pygame events.
        """
        super()._process_event()
        if pygame.key.get_pressed()[K_DOWN]:
            if pygame.key.get_pressed()[K_LSHIFT]:
                self.ly -= 0.01
            else:
                self.ax += 2.0
        if pygame.key.get_pressed()[K_UP]:
            if pygame.key.get_pressed()[K_LSHIFT]:
                self.ly += 0.01
            else:
                self.ax -= 2.0
        if pygame.key.get_pressed()[K_RIGHT]:
            if pygame.key.get_pressed()[K_LSHIFT]:
                self.lx += 0.01
            else:
                self.ay += 2.0
        if pygame.key.get_pressed()[K_LEFT]:
            if pygame.key.get_pressed()[K_LSHIFT]:
                self.lx -= 0.01
            else:
                self.ay -= 2.0

        if pygame.key.get_pressed()[K_PAGEUP]:
            self.lz += 0.01
        if pygame.key.get_pressed()[K_PAGEDOWN]:
            self.lz -= 0.01

        if self.ortho:
            glMatrixMode(GL_PROJECTION)
            glLoadIdentity()
            limit = -self.lz + 1.5
            ratio = self.size[0] / self.size[1]
            if ratio < 1.:
                glOrtho(-limit, limit, -limit / ratio, limit / ratio, -64., 64.)
            elif ratio > 1.:
                glOrtho(-limit * ratio, limit * ratio, -limit, limit, -64., 64.)
            else:
                glOrtho(-limit, limit, -limit, limit, -64., 64.)
            glMatrixMode(GL_MODELVIEW)

        if self.rotating and self.ticking:
            self.ay += self.rot_vel
            self.ay %= 360

    def set_line_width(self, width):
        """描画される直線の太さを設定する.

        Parameters
        ----------
        width : float
             描画に使用される直線の太さ.
        """
        self.linewidth = width
        glLineWidth(self.linewidth)

    def set_point_size(self, size):
        """描画される点の大きさを設定する.

        Parameters
        ----------
        size : float
             描画に使用される点の大きさ.
        """
        self.pointsize = size
        glPointSize(self.pointsize)

    def draw(self, points, mode=GL_POINTS, color=None):
        """頂点配列を指定した描画モードで描画する.

        Parameters
        ----------
        points : allay_like
            頂点の座標を保持する2次元配列.
            2次元目の要素数は3である必要がある.
        mode : OpenGL.constant.IntConstant
            描画モードを表す定数.
            pgutil.GL_POINTS, pgutil.GL_LINES, pgutil.GL_LINE_STRIP,
            pgutil.TRIANGLES, pgutil.TRIANGLE_STRIP, pgutil.GL_QUADS
            や他のOpenGLの描画モードが使用可能.
            デフォルト値はpgutil.GL_POINTSで頂点を点として描画する.
        color : tuple of float, optional
            描画に使用される色.
            [0-1]の値を持つ要素数3または4のタプルが設定可能.
            指定がない場合はPlayground.colに保存された色が使用される.

        """
        if color is None:
            color = list(map(lambda c: c / 255., self.col))
        glEnableClientState(GL_VERTEX_ARRAY)
        if len(color) == 3:
            glColor3dv(color)
        elif len(color) == 4:
            glColor4dv(color)
        glVertexPointerf(points)
        glDrawArrays(mode, 0, len(points))

    def draw_with_color(self, points, mode, colors):
        """頂点配列を指定した描画モードで描画する.

        それぞれの頂点は異なる色で描画される.

        Parameters
        ----------
        points : allay_like
            頂点の座標を保持する2次元配列.
            2次元目の要素数は3である必要がある.
        mode : OpenGL.constant.IntConstant
            描画モードを表す定数.
            pgutil.GL_POINTS, pgutil.GL_LINES, pgutil.GL_LINE_STRIP,
            pgutil.TRIANGLES, pgutil.TRIANGLE_STRIP, pgutil.GL_QUADS
            や他のOpenGLの描画モードが使用可能.
            デフォルト値はpgutil.GL_POINTSで頂点を点として描画する.
        colors : allay_like
            描画に使用する色の情報を格納した配列.
            1番目の次元の要素数はpointsに渡す配列と同じでなければならない.
            すべての値は[0-1]の範囲に入っている必要があり、2番目の次元の要素数は3または4でなければならない.

        """
        glEnableClientState(GL_COLOR_ARRAY)
        glEnableClientState(GL_VERTEX_ARRAY)

        glColorPointerf(colors)
        glVertexPointerf(points)
        glDrawArrays(mode, 0, len(points))

        glDisableClientState(GL_COLOR_ARRAY)
        glDisableClientState(GL_VERTEX_ARRAY)

    def draw_with_norm(self, points, mode, norms):
        """頂点配列を指定した描画モードで描画する.

        それぞれの頂点に法線ベクトルを設定することで陰影をつけて描画する.

        Parameters
        ----------
        points : allay_like
            頂点の座標を保持する2次元配列.
            2次元目の要素数は3である必要がある.
        mode : OpenGL.constant.IntConstant
            描画モードを表す定数.
            pgutil.GL_POINTS, pgutil.GL_LINES, pgutil.GL_LINE_STRIP,
            pgutil.TRIANGLES, pgutil.TRIANGLE_STRIP, pgutil.GL_QUADS
            や他のOpenGLの描画モードが使用可能.
            デフォルト値はpgutil.GL_POINTSで頂点を点として描画する.
        norms : allay_like
            描画に使用する法線ベクトルを格納した配列.
            1番目の次元の要素数はpointsに渡す配列と同じで, 2番目の次元の要素数は3でなければならない.

        """
        if not self.light_initialized: self.init_lighting()

        glEnable(GL_LIGHT0)
        glEnable(GL_NORMALIZE)
        glEnable(GL_LIGHTING)

        glEnableClientState(GL_NORMAL_ARRAY)
        glEnableClientState(GL_VERTEX_ARRAY)

        glNormalPointerf(norms)
        glVertexPointerf(points)
        glDrawArrays(mode, 0, len(points))

        glDisableClientState(GL_NORMAL_ARRAY)
        glDisableClientState(GL_VERTEX_ARRAY)

        glDisable(GL_LIGHT0)
        glDisable(GL_NORMALIZE)
        glDisable(GL_LIGHTING)

    def draw_with_tex(self, points, mode, texcords):
        """頂点配列を指定した描画モードで描画する.

        それぞれの頂点にテクスチャ座標を指定することでポリゴンにテクスチャを貼り付けて描画する.
        テクスチャは事前にバインドしておく必要があり、使用できるテクスチャは2Dテクスチャのみ.

        Parameters
        ----------
        points : allay_like
            頂点の座標を保持する2次元配列.
            2次元目の要素数は3である必要がある.
        mode : OpenGL.constant.IntConstant
            描画モードを表す定数.
            pgutil.GL_POINTS, pgutil.GL_LINES, pgutil.GL_LINE_STRIP,
            pgutil.TRIANGLES, pgutil.TRIANGLE_STRIP, pgutil.GL_QUADS
            や他のOpenGLの描画モードが使用可能.
            デフォルト値はpgutil.GL_POINTSで頂点を点として描画する.
        texcords : allay_like
            描画に使用するテクスチャ座標を格納した配列.
            すべての値は[0-1]の範囲に入っている必要があり、1番目の次元の要素数はpointsに渡す配列と同じで, 2番目の次元の要素数は2でなければならない.

        """

        glEnable(GL_TEXTURE_2D)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)

        glColor4dv((1., 1., 1., 1.))

        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_TEXTURE_COORD_ARRAY)

        glVertexPointerf(points)
        glTexCoordPointerf(texcords)
        glDrawArrays(mode, 0, len(points))

        glDisableClientState(GL_VERTEX_ARRAY)
        glDisableClientState(GL_TEXTURE_COORD_ARRAY)

        glDisable(GL_TEXTURE_2D)

    def draw_with_tex_3d(self, points, mode, texcords):
        """頂点配列を指定した描画モードで描画する.

        それぞれの頂点にテクスチャ座標を指定することでポリゴンにテクスチャを貼り付けて描画する.
        テクスチャは事前にバインドしておく必要があり、使用できるテクスチャは3Dテクスチャのみ.

        Parameters
        ----------
        points : allay_like
            頂点の座標を保持する2次元配列.
            2次元目の要素数は3である必要がある.
        mode : OpenGL.constant.IntConstant
            描画モードを表す定数.
            pgutil.GL_POINTS, pgutil.GL_LINES, pgutil.GL_LINE_STRIP,
            pgutil.TRIANGLES, pgutil.TRIANGLE_STRIP, pgutil.GL_QUADS
            や他のOpenGLの描画モードが使用可能.
            デフォルト値はpgutil.GL_POINTSで頂点を点として描画する.
        texcords : allay_like
            描画に使用するテクスチャ座標を格納した配列.
            すべての値は[0-1]の範囲に入っている必要があり、1番目の次元の要素数はpointsに渡す配列と同じで, 2番目の次元の要素数は3でなければならない.

        """

        glEnable(GL_TEXTURE_3D)
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE)
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE)
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE)

        glColor4dv((1., 1., 1., 1.))

        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_TEXTURE_COORD_ARRAY)

        glVertexPointerf(points)
        glTexCoordPointerf(texcords)
        glDrawArrays(mode, 0, len(points))

        glDisableClientState(GL_VERTEX_ARRAY)
        glDisableClientState(GL_TEXTURE_COORD_ARRAY)

        glDisable(GL_TEXTURE_3D)

    def draw_text(self, text, center, scl=0.3, rotated=False):
        """指定の位置にテキストを描画する.

        Parameters
        ----------
        text : str
            描画するテキスト.
            改行も反映される.
        center : tuple of float
            テキストを描画する座標.
            要素数3のタプルによって指定する.
        scl : float
            表示の際のスケール.
            大きくするほどテキストが大きく表示される.
            デフォルト値は0.3
        rotated : bool
            90°回転するかどうか.
            デフォルトはFalse.
        """

        if text in self.text_cache:
            tex_num, size = self.text_cache[text]
        else:
            tex_num, size = _load_texture(text, self.font)
            self.text_cache[text] = (tex_num, size)

        coords = np.array([[0., 0., 0], [1., 0., 0], [1., 1., 0], [0., 1., 0]]) * scl
        coords[:, 0] *= size[0] / size[1] * 1.16
        texcoords = np.array([[1., 0.], [1., 1.], [0., 1.], [0., 0.]])

        glBindTexture(GL_TEXTURE_2D, tex_num)
        glDepthMask(GL_FALSE)

        glPushMatrix()
        glTranslated(center[0], center[1], center[2])
        glRotated(-self.ay, 0.0, 1.0, 0.0)
        glRotated(-self.ax, 1.0, 0.0, 0.0)
        if rotated: glRotated(90., 0., 0., 1.)
        glTranslated(-size[0] / size[1] * scl * 0.5, -scl * 0.5, 0.)
        self.draw_with_tex(coords, GL_QUADS, texcoords)
        glPopMatrix()

        glDepthMask(GL_TRUE)

    def draw_cube(self, l, color=None):
        """それぞれの辺がx, y, z軸に廃校した立方体を描画する.

        Parameters
        ----------
        l : float
            描画する立方体の1辺の長さの1/2.
            最終的に描画される立方体の1辺の長さは2lとなる.
        color : tuple of float, optional
            描画に使用される色.
            [0-1]の値を持つ要素数3または4のタプルが設定可能.
            指定がない場合はPlayground.colに保存された色が使用される.

        """

        if color is None:
            color = list(map(lambda c: c / 255., self.col))

        self.draw(((-l, -l, -l), (l, -l, -l), (l, l, -l), (-l, l, -l)), GL_LINE_LOOP, color)
        self.draw(((-l, -l, l), (l, -l, l), (l, l, l), (-l, l, l)), GL_LINE_LOOP, color)
        self.draw(((-l, -l, -l), (-l, -l, l), (l, -l, -l), (l, -l, l), (l, l, -l), (l, l, l), (-l, l, -l), (-l, l, l)),
                  GL_LINES,
                  color)

    def draw_axis(self):
        """x, y, z軸に廃校したRGB3色の直線を描画する.

        デバッグ時の座標の確認用.
        """

        d = 1.1
        self.draw_with_color(((-d, -d, -d), (d, -d, -d), (-d, -d, -d), (-d, d, -d), (-d, -d, -d), (-d, -d, d)),
                             GL_LINES,
                             ((1., 0., 0.), (1., 0., 0.), (0., 1., 0.), (0., 1., 0.), (0., 0., 1.), (0., 0., 1.)))

    def gen_voxel_data(self, array, colors=None):
        """3次元配列をボクセルデータとして描画するためのデータを生成する.

        Parameters
        ----------
        array : allay_like
            ボクセルデータが格納されたnp.bool型の3次元配列.
            Trueに対応する場所が不透明であるとして描画される.
        colors : array_like
            描画に使用する色の情報を含んだ配列.
            1次元目の要素数は3または6, 2次元目の要素数は3または4である必要がある.

        Returns
        -------
        verts : allay_like
            生成されたボクセルの頂点の座標が格納された配列.
            GLground.draw_with_color(verts, pgutil.GL_QUADS, cols)のようにして描画に使用できる.
        cols : allay_like
            生成されたボクセルの色の情報が格納された配列.
            GLground.draw_with_color(verts, pgutil.GL_QUADS, cols)のようにして描画に使用できる.

        """

        array = np.swapaxes(array, 0, 1)
        edges = _get_edges(array)
        verts = np.vstack(edges)

        col = colors if colors is not None else self.voxelcol
        cols = np.vstack([np.tile(col[i//2] if len(col)==3 else col[i],(edges[i].shape[0],1)) for i in range(6)])
        return verts, cols


    def draw_voxel(self, array, colors=None):
        """3次元配列をボクセルデータとして描画する.

        Parameters
        ----------
        array : allay_like
            ボクセルデータが格納されたnp.bool型の3次元配列.
            Trueに対応する場所が不透明であるとして描画される.
        colors : array_like
            描画に使用する色の情報を含んだ配列.
            1次元目の要素数は3または6, 2次元目の要素数は3または4である必要がある.

        """

        verts, cols = self.gen_voxel_data(array, colors)
        self.draw_with_color(verts, GL_QUADS, cols)

    def gen_meshdata(self, array, minin=None, maxin=None, scl=1.):
        """2次元配列をメッシュデータとして描画するためのデータを生成する.

        Parameters
        ----------
        array : array_like
            メッシュとして描画する対象となるデータを格納した2次元配列.
            2次元配列でなおかつ1次元目と2次元目の要素数が一致している必要がある.
        minin : float, optional
            表示の際に最小として表示される値.
            これより小さい値は黒として表示される.
            指定がない場合は入力された配列の最小値が使用される.
        maxin : float, optional
            表示の際に最大として表示される値.
            これより大きい値はPlayground.colと同じ色で表示される.
            指定がない場合は入力された配列の最大値が使用される.
        scl : float, optional
            高さ方向のスケールを決定する定数.
            指定がない場合は1が使用され、最大、最小の高さがそれぞれ1, -1になるよう描画される.

        Returns
        -------
        mesh : allay_like
            生成されたメッシュの頂点の座標が格納された配列.
            GLground.draw(mesh, pgutil.GL_LINE_STRIP)のようにして描画に使用できる.

        """

        if minin is not None or maxin is not None:
            array = np.clip(array, minin, maxin)

        ma = maxin if maxin is not None else array.max()
        mi = minin if minin is not None else array.min()

        temp = ((array - mi) / ((ma - mi) / 2)).ravel() - 1.

        if not array.shape in self.mesh_cache:
            self.mesh_cache[array.shape] = _gen_mesh_index(array.shape)

        mesh = self.mesh_cache[array.shape]

        mesh[1][1, :] = temp * scl
        return mesh[1].T[mesh[0]].copy()

    def draw_mesh(self, array, minin=None, maxin=None, scl=1.):
        """2次元配列をメッシュデータとして描画する.

        Parameters
        ----------
        array : array_like
            メッシュとして描画する対象となるデータを格納した2次元配列.
            2次元配列でなおかつ1次元目と2次元目の要素数が一致している必要がある.
        minin : float, optional
            表示の際に最小として表示される値.
            これより小さい値は黒として表示される.
            指定がない場合は入力された配列の最小値が使用される.
        maxin : float, optional
            表示の際に最大として表示される値.
            これより大きい値はPlayground.colと同じ色で表示される.
            指定がない場合は入力された配列の最大値が使用される.
        scl : float, optional
            高さ方向のスケールを決定する定数.
            指定がない場合は1が使用され、最大、最小の高さがそれぞれ1, -1になるよう描画される.

        """

        self.draw(self.gen_meshdata(array, minin, maxin, scl), GL_LINE_STRIP)

    def gen_meshdata_colored(self, array, minin=None, maxin=None, cmap="gist_rainbow_r", scl=1.):
        """2次元配列を色付きのメッシュデータとして描画するためのデータを生成する.

        Parameters
        ----------
        array : array_like
            メッシュとして描画する対象となるデータを格納した2次元配列.
            2次元配列でなおかつ1次元目と2次元目の要素数が一致している必要がある.
        minin : float, optional
            表示の際に最小として表示される値.
            これより小さい値は黒として表示される.
            指定がない場合は入力された配列の最小値が使用される.
        maxin : float, optional
            表示の際に最大として表示される値.
            これより大きい値はPlayground.colと同じ色で表示される.
            指定がない場合は入力された配列の最大値が使用される.
        cmap : str
            描画に使用されるカラーマップ.
            matplotlibで使用可能なカラーマップ全てと"AVS", "AVS_r"が使用可能.
            指定がない場合は"gist_rainbow_r"が選択される.
        scl : float, optional
            高さ方向のスケールを決定する定数.
            指定がない場合は1が使用され、最大、最小の高さがそれぞれ1, -1になるよう描画される.

        Returns
        -------
        mesh : allay_like
            生成されたメッシュの頂点の座標が格納された配列.
            GLground.draw_with_color(mesh, pgutil.GL_LINE_STRIP, color)のようにして描画に使用できる.
        color : allay_like
            生成されたメッシュの色の情報が格納された配列.
            GLground.draw_with_color(mesh, pgutil.GL_LINE_STRIP, color)のようにして描画に使用できる.

        """

        if minin is not None or maxin is not None:
            array = np.clip(array, minin, maxin)

        ma = maxin if maxin is not None else array.max()
        mi = minin if minin is not None else array.min()

        temp = ((array - mi) / ((ma - mi))).ravel()

        color = self._get_cmap(cmap)[(temp * 255).astype(np.uint8)]

        if not array.shape in self.mesh_cache:
            if len(array.shape) != 2 or array.shape[0] != array.shape[1]:
                raise Exception("Illegal array shape.")
            self.mesh_cache[array.shape] = _gen_mesh_index(array.shape)

        mesh = self.mesh_cache[array.shape]

        mesh[1][1, :] = (temp - 0.5) * (2 * scl)
        return mesh[1].T[mesh[0]].copy(), color[mesh[0]]

    def draw_mesh_colored(self, array, minin=None, maxin=None, cmap="gist_rainbow_r", scl=1.):
        """2次元配列を色付きのメッシュデータとして描画する.

        Parameters
        ----------
        array : array_like
            メッシュとして描画する対象となるデータを格納した2次元配列.
            2次元配列でなおかつ1次元目と2次元目の要素数が一致している必要がある.
        minin : float, optional
            表示の際に最小として表示される値.
            これより小さい値は黒として表示される.
            指定がない場合は入力された配列の最小値が使用される.
        maxin : float, optional
            表示の際に最大として表示される値.
            これより大きい値はPlayground.colと同じ色で表示される.
            指定がない場合は入力された配列の最大値が使用される.
        cmap : str
            描画に使用されるカラーマップ.
            matplotlibで使用可能なカラーマップ全てと"AVS", "AVS_r"が使用可能.
            指定がない場合は"gist_rainbow_r"が選択される.
        scl : float, optional
            高さ方向のスケールを決定する定数.
            指定がない場合は1が使用され、最大、最小の高さがそれぞれ1, -1になるよう描画される.

        """

        mesh, color = self.gen_meshdata_colored(array, minin, maxin, cmap, scl)
        self.draw_with_color(mesh, GL_LINE_STRIP, color)

    def gen_surfdata(self, array, minin=None, maxin=None, scl=1.):
        """2次元配列を陰影が処理された曲面として描画するためのデータを生成する.

        Parameters
        ----------
        array : array_like
            曲面として描画する対象となるデータを格納した2次元配列.
            2次元配列でなおかつ1次元目と2次元目の要素数が一致している必要がある.
        minin : float, optional
            表示の際に最小として表示される値.
            これより小さい値は黒として表示される.
            指定がない場合は入力された配列の最小値が使用される.
        maxin : float, optional
            表示の際に最大として表示される値.
            これより大きい値はPlayground.colと同じ色で表示される.
            指定がない場合は入力された配列の最大値が使用される.
        scl : float, optional
            高さ方向のスケールを決定する定数.
            指定がない場合は1が使用され、最大、最小の高さがそれぞれ1, -1になるよう描画される.

        Returns
        -------
        surf : allay_like
            生成された曲面の頂点の座標が格納された配列.
            GLground.draw_with_norm(surf, pgutil.GL_TRIANGLE_STRIP, norm)のようにして描画に使用できる.
        norms : allay_like
            生成された曲面の法線ベクトルが格納された配列.
            GLground.draw_with_norm(surf, pgutil.GL_TRIANGLE_STRIP, norm)のようにして描画に使用できる.

        """

        array = np.swapaxes(array, 0, 1)
        if minin is not None or maxin is not None:
            array = np.clip(array, minin, maxin)

        ma = maxin if maxin is not None else array.max()
        mi = minin if minin is not None else array.min()

        temp = ((array - mi) / ((ma - mi) / 2 / scl) - 1. * scl).ravel()

        if not array.shape in self.surf_cache:
            if len(array.shape) != 2 or array.shape[0] != array.shape[1]:
                raise Exception("Illegal array shape.")
            drawing = np.zeros(((array.shape[0] * (array.shape[1] - 1)) * 2, 3))
            xx, yy = np.meshgrid(np.linspace(-1, 1, array.shape[0]), np.linspace(-1, 1, array.shape[1]))
            drawing[::2, 0] = xx.ravel()[:-array.shape[0]]
            drawing[1::2, 0] = xx.ravel()[array.shape[0]:]
            drawing[::2, 2] = yy.ravel()[:-array.shape[0]]
            drawing[1::2, 2] = yy.ravel()[array.shape[0]:]
            drawing[::array.shape[0] * 2, 0] = np.nan
            norms = np.zeros_like(drawing)
            self.surf_cache[array.shape] = (drawing, norms)

        surf = self.surf_cache[array.shape][0]
        surf[::2, 1] = temp[:-array.shape[0]]
        surf[1::2, 1] = temp[array.shape[0]:]

        norms = self.surf_cache[array.shape][1]
        norms[1:-1, :] = np.cross(surf[0:-2, :] - surf[1:-1, :], surf[1:-1, :] - surf[2:, :])
        norms[::2] *= -1

        return surf.copy(), norms.copy()

    def draw_surface(self, array, minin=None, maxin=None, scl=1.):
        """2次元配列を陰影が処理された曲面として描画する.

        Parameters
        ----------
        array : array_like
            曲面として描画する対象となるデータを格納した2次元配列.
            2次元配列でなおかつ1次元目と2次元目の要素数が一致している必要がある.
        minin : float, optional
            表示の際に最小として表示される値.
            これより小さい値は黒として表示される.
            指定がない場合は入力された配列の最小値が使用される.
        maxin : float, optional
            表示の際に最大として表示される値.
            これより大きい値はPlayground.colと同じ色で表示される.
            指定がない場合は入力された配列の最大値が使用される.
        scl : float, optional
            高さ方向のスケールを決定する定数.
            指定がない場合は1が使用され、最大、最小の高さがそれぞれ1, -1になるよう描画される.

        """

        surf, norm = self.gen_surfdata(array, minin, maxin, scl)
        self.draw_with_norm(surf, GL_TRIANGLE_STRIP, norm)

    def gen_marching_data(selfarray, array, thresh=None, reverse=False):
        """マーチングキューブ法を用いて3次元配列の等値面をを滑らかな立体として描画するためのデータを生成する.

        Parameters
        ----------
        array : array_like
            立体として描画するデータを格納した3次元配列.
        thresh : float
            描画する曲面での配列中での値.
        reverse : bool, optional
            不透明であるとして描画する範囲を反転させるかどうか.
            Falseの場合, threshに指定された値より大きい値を持つ場所が不透明であるとして描画される.
            Trueの場合, threshに指定された値より小さい値を持つ場所が不透明であるとして描画される.
            デフォルト値はFalse.

        Returns
        -------
        surf : allay_like
            生成された曲面の頂点の座標が格納された配列.
            GLground.draw_with_norm(surf, pgutil.GL_TRIANGLE_STRIP, norm)のようにして描画に使用できる.
        norms : allay_like
            生成された曲面の法線ベクトルが格納された配列.
            GLground.draw_with_norm(surf, pgutil.GL_TRIANGLE_STRIP, norm)のようにして描画に使用できる.

        """

        width_half = array.shape[0] / 2
        try:
            verts, faces, normals, values = measure.marching_cubes_lewiner(array, thresh)
        except (ValueError, RuntimeError):
            return None, None

        new_norm = -normals[faces] if reverse else normals[faces]
        return (((verts[faces] - (width_half - 0.5)) / (width_half - 0.5))).reshape(faces.shape[0] * 3, 3), new_norm

    def draw_marching_cubes(self, array, thresh=None, reverse=False):
        """マーチングキューブ法を用いて3次元配列の等値面をを滑らかな立体として描画する.

        Parameters
        ----------
        array : array_like
            立体として描画するデータを格納した3次元配列.
        thresh : float
            描画する曲面での配列中での値.
        reverse : bool, optional
            不透明であるとして描画する範囲を反転させるかどうか.
            Falseの場合, threshに指定された値より大きい値を持つ場所が不透明であるとして描画される.
            Trueの場合, threshに指定された値より小さい値を持つ場所が不透明であるとして描画される.
            デフォルト値はFalse.

        """

        surf, norm = self.gen_marching_data(array, thresh, reverse)

        if surf is not None and norm is not None:
            self.draw_with_norm(surf, GL_TRIANGLES, norm)

    def bind_tex_cube(self, array, minin=None, maxin=None, cmap=None):
        """3Dテクスチャを用いて3次元配列を模様付きの立方体として描画するためのテクスチャをバインドする.

        この関数を1度呼び出してからGLground.draw_tex_cubeをarrayの引数を指定無しで呼び出すことで内容の変化しない3次元配列を光速に描画可能.

        Parameters
        ----------
        array : array_like
            描画するデータを格納した3次元配列もしくは4次元配列.
            3次元配列の場合任意の値を持った配列で, 値に応じた色がカラーマップから選ばれて描画される.
            4次元配列の場合は4次元目がRGBの順番で色の情報を格納した配列であるとみなす. この場合, 配列の値は[0-1]の範囲に入っている必要がある.
        minin : float, optional
            表示の際に最小として表示される値.
            これより小さい値は黒として表示される.
            指定がない場合は入力された配列の最小値が使用される.
            arrayに渡された配列の次元数が4の場合は無視される.
        maxin : float, optional
            表示の際に最大として表示される値.
            これより大きい値はPlayground.colと同じ色で表示される.
            指定がない場合は入力された配列の最大値が使用される.
            arrayに渡された配列の次元数が4の場合は無視される.
        cmap : str, optional
            描画に使用されるカラーマップ.
            matplotlibで使用可能なカラーマップ全てと"AVS", "AVS_r"が使用可能.
            指定がない場合は"vivid"が選択される.
            arrayに渡された配列の次元数が4の場合は無視される.

        Returns
        -------
        color : allay_like
            生成された色情報を格納した4次元配列.

        """

        self.last_bind_func = lambda: self.bind_tex_cube(array.copy(), minin, maxin, cmap)

        glBindTexture(GL_TEXTURE_3D, self.cube_tex)
        array = np.swapaxes(array, 0, 2)
        w, h, d = array.shape[:3]
        if array.ndim == 3:
            if minin is not None or maxin is not None:
                array = np.clip(array, minin, maxin)

            ma = maxin if maxin is not None else array.max()
            mi = minin if minin is not None else array.min()

            temp = ((array.astype(np.float64) - mi) * (255 / ((ma - mi)))).astype(np.uint8)
            color = (255 * self._get_cmap(cmap)).astype(np.uint8)[temp]

            glTexImage3D(GL_TEXTURE_3D, 0, GL_RGB, w, h, d, 0, GL_RGB, GL_UNSIGNED_BYTE, color)
        elif array.ndim == 4:
            color = array
        else:
            raise Exception("Dimension of tex array must be 3 or 4")

        glTexImage3D(GL_TEXTURE_3D, 0, GL_RGB, w, h, d, 0, GL_RGB, GL_UNSIGNED_BYTE, color)

        return color

    def draw_tex_cube(self, array=None, minin=None, maxin=None, cmap=None, mouse=False, box=(-1, 1, -1, 1, -1, 1)):
        """3Dテクスチャを用いて3次元配列を模様付きの立方体として描画する.

        Parameters
        ----------
        array : array_like
            描画するデータを格納した3次元配列もしくは4次元配列.
            3次元配列の場合任意の値を持った配列で, 値に応じた色がカラーマップから選ばれて描画される.
            4次元配列の場合は4次元目がRGBの順番で色の情報を格納した配列であるとみなす. この場合, 配列の値は[0-1]の範囲に入っている必要がある.
            Noneの場合, 事前にバインドしたテクスチャがそのまま描画に使用される.
        minin : float, optional
            表示の際に最小として表示される値.
            これより小さい値は黒として表示される.
            指定がない場合は入力された配列の最小値が使用される.
            arrayに渡された配列の次元数が4の場合は無視される.
        maxin : float, optional
            表示の際に最大として表示される値.
            これより大きい値はPlayground.colと同じ色で表示される.
            指定がない場合は入力された配列の最大値が使用される.
            arrayに渡された配列の次元数が4の場合は無視される.
        cmap : str, optional
            描画に使用されるカラーマップ.
            matplotlibで使用可能なカラーマップ全てと"AVS", "AVS_r"が使用可能.
            指定がない場合は"vivid"が選択される.
            arrayに渡された配列の次元数が4の場合は無視される.
        mouse : bool, optional
            描画する立方体の断面をマウス操作によって移動できるようにするかどうか.
            デフォルト値はFalse.
        box : tuple of float
            描画する立方体の各軸の最小, 最大の指定.
            mouseがTrueの場合は無視される.

        """

        if array is not None:
            self.bind_tex_cube(array, minin, maxin, cmap)

        if mouse:
            verts, tex_cords = _gen_cube(-1, self.mouseX / self.size[0] * 2 - 1, -1, 1 - self.mouseY / self.size[1] * 2,
                                         -1, 1)
        else:
            verts, tex_cords = _gen_cube(*box)

        self.draw_with_tex_3d(verts, GL_QUADS, tex_cords)

    def init_lighting(self):
        """照明システムを初期化する.

        法線ベクトルを使用した描画の前に呼び出される必要があるが, ほとんどの場合は自動で呼び出される.

        """
        black = [0.0, 0.0, 0.0, 1.0]
        white = [1.0, 1.0, 1.0, 1.0]
        ambient = [0.1, 0.1, 0.1, 1.0]
        diffuse = [0.9, 0.9, 0.9, 1.0]
        lightPos = [0.0, 0.0, 1.0, 0.0]

        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient)
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse)
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, black)

        glLightfv(GL_LIGHT0, GL_AMBIENT, white)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, white)
        glLightfv(GL_LIGHT0, GL_SPECULAR, white)
        glLightfv(GL_LIGHT0, GL_POSITION, lightPos)

        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, black)

        glEnable(GL_NORMALIZE)

        if not self.light_initialized: self.screenmodefuncs.append(self.init_lighting)

        self.light_initialized = True

    def clear_tex_cache(self):
        """テクスチャ番号のキャッシュを削除する.

        画面モードが切り替わった際に自動で呼び出される.

        """
        self.text_cache = {}

    def _reload_3d_tex(self):
        if self.last_bind_func is not None:
            self.last_bind_func()


from pygame import *
import pygame.gfxdraw as gfxdraw
