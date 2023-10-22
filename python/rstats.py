import numpy as np
import pandas as pd
import math
from math import sqrt as sqrt
from math import pi as pi
from math import exp as exp
from scipy.special import erfinv, erf

import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from rutil import pwl_root_crossing_segmentStats

# Hurst package
import numpy as np
import matplotlib.pyplot as plt
from hurst import compute_Hc, random_walk

# Higuchi Fractal Dimension (HFD)
import HiguchiFractalDimension as hfd

# PSD
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec

# autocorrelation function
from scipy import signal

def Test_HD_Cor():
    cl = -2.5
    fn = '../../InhomogeneousFiles/cl' + str(cl) + '_np16385/initial_values_0.txt'
    corl = 10**cl
    # fn = '../../InhomogeneousFiles/cl-1.5_np16385/initial_values_0.txt'
    # fn = '../../InhomogeneousFiles/cl-0.5_np16385/initial_values_0.txt'
    vals = np.loadtxt(fn)[1:]
    valsOrig = vals
    mn = np.mean(vals)
    sd = np.std(vals)
    print(f"mean = {mn}\n")
    print(f"sd = {sd}\n")
    sz = len(vals)
    vals = (vals - mn) / sd
    dl = 1 / (sz - 1)
    xs = np.arange(0, 1 + dl, dl)
    seg_num, seg_ave, seg_mn, seg_mx, seg_sdiv = pwl_root_crossing_segmentStats(xs, vals)
    print(f"seg_num = {seg_num:.5f}\n")
    print(f"seg_ave = {seg_ave:.5f}\n")
    print(f"seg_mn = {seg_mn:.5f}\n")
    print(f"seg_mx = {seg_mx:.5f}\n")
    print(f"seg_sdiv = {seg_sdiv:.5f}\n")

#        rng = np.random.default_rng()
#        sig = rng.standard_normal(1000)
    autocorr = signal.fftconvolve(vals, vals[::-1], mode='full')
    ac = autocorr[sz - 1:2 * sz]
    ac = ac / ac[0]
    sz2 = len(ac)
    dl = 1 / (sz2 - 1)
    x = np.arange(0, 1 + dl, dl)
    first_negative_index = np.argmax(ac < 0)
    ac2abs = abs(ac)
    ac2abs2 = ac2abs**2
    l1all = np.trapz(ac2abs, x)
    l2all = np.trapz(ac2abs2, x)
    l1allRed = l1all
    l2allRed = l2all
    xCorNeg = np.nan
    if (first_negative_index > 0):
        xCorNeg = x[first_negative_index]
        xRed = x[0:first_negative_index]
        acRed = ac[0:first_negative_index]
        acRed2abs = abs(acRed)
        acRed2abs2 = acRed2abs**2
        l1allRed = np.trapz(acRed2abs, xRed)
        l2allRed = np.trapz(acRed2abs2, xRed)
    print(f"xCorNeg = {xCorNeg:.5f}\n")
    print(f"L1(ac) = {l1all:.5f}\n")
    print(f"L2(ac) = {l2all:.5f}\n")
    print(f"L1(acRed) = {l1allRed:.5f}\n")
    print(f"L2(acRed) = {l2allRed:.5f}\n")
    print(f"L1(acRedPredicted) = {corl * np.sqrt(np.pi / 2):.5f}\n")
    print(f"L2(acRedPredicted) = {corl * np.sqrt(np.pi / 2)/2:.5f}\n")


    autocorr = autocorr / sz
    import matplotlib.pyplot as plt
    fig, (ax_orig, ax_mag) = plt.subplots(2, 1)
    ax_orig.plot(x, valsOrig)
    ax_orig.set_title('random field')
    #sz2 = sz #5000
    xm = min(100 * corl, 1)
    # ax_mag.plot(np.arange(-sz+1,sz)/sz, autocorr)
    ax_mag.plot(x, ac, label="computed")
    acPredicted = np.exp(-(x/corl)**2)
    ax_mag.plot(x, acPredicted, label="original")
    plt.xlim(0, xm)
    plt.legend()
    ax_mag.set_title('Autocorrelation')
    fig.tight_layout()
    fig.show()




    if False:
        step = 1/sz
        x = 0 * vals
        for i in range(1,sz):
            x[i] = i * step 
        plt.subplot(311)
        plt.plot(x, vals)
        plt.subplot(313)
        plt.psd(vals, 512, 1 / step)
        plt.show()

    """
    x = np.random.randn(10000)
    y = np.empty(9900)
    for i in range(x.size-100):
        y[i] = np.sum(x[:(i+100)])
    hfd.hfd(x)
    """

    hfd_v = hfd.hfd(vals)
    print(f"hfd =  {hfd_v:.3f}")

    # Use random_walk() function or generate a random walk series manually:
    # series = random_walk(99999, cumprod=True)
    np.random.seed(42)
    random_changes = 1. + np.random.randn(99999) / 1000.
    #vals = np.cumprod(random_changes)  # create a random walk from random changes
    #vals = np.power(10.0, vals)

    # Evaluate Hurst equation
    # H, c, data = compute_Hc(vals, kind='price', simplified=True)
    H, c, data = compute_Hc(vals, kind='change', simplified=True)

    f, ax = plt.subplots()
    ax.plot(data[0], c*data[0]**H, color="deepskyblue")
    ax.scatter(data[0], data[1], color="purple")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Time interval')
    ax.set_ylabel('R/S ratio')
    ax.grid(True)
    plt.show()

    print("H={:.4f}, c={:.4f}".format(H,c))
    return


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













