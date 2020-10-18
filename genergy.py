# coding: UTF-8
import numpy as np
import matplotlib.pyplot as plt
lalpha = 18500.
RR = 8.31
TT = 1000.
c = np.arange(0., 1., 0.01)
for i in range(24):
    del_lalpha = 500. * i
    galpa = (lalpha + del_lalpha) * c * (1. - c) + RR * TT * (c * np.log(c) + (1. - c) * np.log(1. - c))
    plt.plot(c, galpa,label=f"L = {int(lalpha + del_lalpha)}")
plt.legend()
plt.show()
