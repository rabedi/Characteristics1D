import numpy as np
import pandas as pd
import math
from math import sqrt as sqrt
from math import pi as pi
from math import exp as exp
from scipy.special import erfinv, erf

import matplotlib.pyplot as plt
import matplotlib.lines as mlines


class StatOfVec:
    def __init__(self):
        self.delta = 0.5
        self.alpha = 2.0
        self.lc = 0.01
        #self.n = 512
        self.n = 16384
        #self.n = 128
        # self.n = 10 ** 10 -> For Gumbel test it gives m = 1.95, sigma = 5.49e=6. b = 1.81, delta = 0.5244 FOR delta = 0.5 and alpha = 2 which is pretty good. For other values it's not that good

    
    def Test_Gumbel2Weibull(self):
        # Gumbel has n realization for SN, want to see if the pull back to this distribution works
        n = self.n
        z0 = sqrt(2.0) * erfinv(2.0 / n - 1)
        pz0 = 1.0 / sqrt(2.0 * pi) * exp(-0.5 * z0 * z0)
        Npz0 = n * pz0
        lss10 = np.arange(-7.0, -4.0, 0.01)
        sz = lss10.shape[0]
        lss = np.zeros(sz, dtype=float)
        pdfwt = np.zeros(sz, dtype=float)
        ss = np.zeros(sz, dtype=float)
        wts = np.zeros(sz, dtype=float)
        cdfs = np.zeros(sz, dtype=float)
        sqrt2 = sqrt(2.0)
        for i, ls10 in enumerate(lss10):
            s = 10 ** ls10
            ss[i] = s
            lss[i] = np.log(s)
            arg = -1 + (s / self.delta) ** self.alpha
            z = sqrt2 * erfinv(arg)
            test_s = self.delta * (1.0 + erf(z / sqrt2)) ** (1.0 / self.alpha)
            weibullTest = -Npz0 * (z0 - z)

            cdf_sn_z = 0.5 * (1 + erf(z / sqrt2))
            prob = 1 - (1 - cdf_sn_z) ** n
            weibullTestB = np.log(-np.log(1 - prob))

            wts[i] = weibullTest
            cdfs[i] = 1 - exp(-exp(-weibullTest))
        
        m_w = (wts[-1] - wts[0])/(lss[-1] - lss[0])
        log_sigma = lss[0] - wts[0] / m_w
        sigma_w = exp(log_sigma)
        b = 1.0 / n * sigma_w ** -m_w 
        recovered_delta = 1.0 / (2.0 * b) ** (1.0 / self.alpha)
        print(f'm {m_w}')
        print(f'sigma {sigma_w}')
        print(f'b {b}')
        print(f'recovered_delta {recovered_delta}')

        plt.plot(lss, wts)

        # Adding labels and title
        plt.xlabel('log(s)')
        plt.ylabel('log(-log(1-P))')
        plt.title('Weibull Test')

        # Display the plot
        plt.show()













