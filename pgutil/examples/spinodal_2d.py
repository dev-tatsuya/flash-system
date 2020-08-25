# coding: UTF-8

import pgutil
import numpy as np

pg = pgutil.Playground()

ND = 64
al = 80.0 * 1.e-9
b1 = al / ND

c0 = 0.5

temp = 1000.
rr = 8.3145
delt = 0.008
rtemp = rr * temp

cmob = 1.0
L = 25000. / rtemp
kappa = 5.0e-15 / b1 ** 2 / rtemp


def laplacian_2d(array):
    result = -4 * array
    result += np.roll(array, 1, 0)
    result += np.roll(array, -1, 0)
    result += np.roll(array, 1, 1)
    result += np.roll(array, -1, 1)
    return result


def chemical(cc):
    return L * (1. - 2. * cc) + (np.log(cc) - np.log(1. - cc))


def dcdt(cc):
    return laplacian_2d(chemical(cc) - 2. * kappa * laplacian_2d(cc))


def step(cc):
    return cc + delt * dcdt(cc)


@pg.keyfunc("F1")
def init():
    global c
    c = c0 + 0.01 * (2 * np.random.rand(ND, ND) - 1.)


init()
while pg.run():
    if pg.ticking:
        for i in range(10):
            c = step(c)

    pg.transform_blit(c, 0, 1)
