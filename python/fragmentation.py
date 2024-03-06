import numpy as np
import pandas as pd
#import pickle
import glob
import math
import os
import statistics
from pathlib import Path

# Hurst package
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

class frgT:
    sz = 4
    names = [""] * sz
    lstyle = [""] * sz
    Zhu6a = 0
    names[Zhu6a] = 'Zhu-a'
    lstyle[Zhu6a] = "solid"
    Zhu6b = 1
    names[Zhu6b] = 'Zhu-b'
    lstyle[Zhu6b] = "dashed"
    Glenn = 2
    names[Glenn] = 'Glenn'
    lstyle[Glenn] = "dashdot"
    Grady = 3
    names[Grady] = 'Grady'
    lstyle[Grady] = "dotted"

# one loading rate results
class Frag1D1la:
    def __init__(self):
        # either one > 0 will be used
        # loading rate
        self.a = -1
        # normalized loading rate
        self.ap = -1
        # length
        self.L = 1
        # G factor (0.5 for Ortiz, e for Xu-Needleman)
        self.energyScale = 0.5
        # fragmentation normalized sizes
        self.s_ps = []
        self.lines_ps = []
        # normalized energy dissipation
        self.e_ps = []
        self.lines_es = []


    def ComputeVals(self):
        zhu6_epsilonDotScale = 1.0 / self.energyScale
        zhu6_sBarScale = self.energyScale
        if (self.ap < 0):
            self.ap = self.a / zhu6_epsilonDotScale
        self.s_ps = np.zeros(frgT.sz)    
        self.s_ps[frgT.Zhu6a] = zhu6_sBarScale * 4.5 / (1.0 + 6.00 * np.power(self.ap, 2.0 / 3.0)) # //(
        self.s_ps[frgT.Zhu6b] = zhu6_sBarScale * 4.5 / (1.0 + 0.77 * pow(self.ap, 0.25) + 5.4 * pow(self.ap, 0.75)) # (29a)
        self.s_ps[frgT.Grady] = zhu6_sBarScale * pow(24.0 / self.ap / self.ap, 1.0 / 3.0) # (1')
        self.s_ps[frgT.Glenn] = zhu6_sBarScale * 4.0 / self.ap * np.sinh(1.0 / 3 * math.asinh(1.5 * self.ap)) # (2')

        self.e_ps = np.zeros(frgT.sz)    
        for i in range(frgT.sz):
            self.e_ps[i] = 1.0 / self.s_ps[i]
 

class Frag1Dlap_curves:
    def __init__(self):
        self.lap_min = -5
        self.lap_max = 4
        self.lap_del = 0.02
        self.laps = []
        self.laps_min = 0
        self.laps_max = 0

        self.aps = []

        self.lsps = []
        self.lines_ps = []
        self.lines_ps_ymin = 0
        self.lines_ps_ymax = 0

        self.leps = []
        self.lines_es = []
        self.lines_es_ymin = 0
        self.lines_es_ymax = 0

    def ReturnLines(self, isLog_phid, isLog_spAve):
        self.Compute(0)
        if (isLog_phid):
            return True, self.lines_es, self.laps_min, self.laps_max, self.lines_es_ymin, self.lines_es_ymax
        if (isLog_spAve):
            return True, self.lines_ps, self.laps_min, self.laps_max, self.lines_ps_ymin, self.lines_ps_ymax 
        return False, [], [], [], [], []

    def Compute(self, doPlots = 0):
        clr = 'darkgray'
        self.laps = np.arange(self.lap_min, self.lap_max, self.lap_del)
        self.laps_min = np.min(self.laps)
        self.laps_max = np.max(self.laps)

        self.aps = np.power(10.0, self.laps)
        n = len(self.aps)
        self.lsps = np.zeros([frgT.sz, n])
        self.leps = np.zeros([frgT.sz, n])
        self.lines_ps_ymin = 1e100
        self.lines_ps_ymax = -1e100
        self.lines_es_ymin = 1e100
        self.lines_es_ymax = -1e100

        for i in range(0, n):
            ap = self.aps[i]
            fd = Frag1D1la()
            fd.ap = ap
            fd.ComputeVals()
            for j in range(0, frgT.sz):
                tmp = np.log10(fd.s_ps[j])
                if (tmp < self.lines_ps_ymin):
                    self.lines_ps_ymin = tmp
                if (tmp > self.lines_ps_ymax):
                    self.lines_ps_ymax = tmp

                self.lsps[j][i] = tmp

                tmp = np.log10(fd.e_ps[j])
                if (tmp < self.lines_es_ymin):
                    self.lines_es_ymin = tmp
                if (tmp > self.lines_es_ymax):
                    self.lines_es_ymax = tmp
                self.leps[j][i] = tmp

        self.lines_ps = []
        self.lines_es = []
        for j in range(0, frgT.sz):
            linep = mlines.Line2D([], [])
            linep.set_xdata(self.laps)
            linep.set_ydata(self.lsps[j])
            linep.set_color(clr)
            linep.set_linestyle(frgT.lstyle[j])
            linep.set_label(frgT.names[j])
            self.lines_ps.append(linep)

            linee = mlines.Line2D([], [])
            linee.set_xdata(self.laps)
            linee.set_ydata(self.leps[j])
            linee.set_color(clr)
            linee.set_linestyle(frgT.lstyle[j])
            linee.set_label(frgT.names[j])
            self.lines_es.append(linee)

        if (doPlots == 0):
            return
        # Create two plots
        plt.figure(figsize=(10, 5))

        # Plot the linear relationship

        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
        for j in range(0, frgT.sz):
            ax1.add_line(self.lines_ps[j])

#        for j in range(0, frgT.sz):
#            plt.plot(self.laps,self.lsps[j], label = frgT.names[j], linestyle = frgT.lstyle[j], color = clr)
        # ax1.set_xlabel('x')
        # ax1.set_ylabel('y')
            
        # recompute the ax.dataLim
        ax1.relim()
        # update ax.viewLim using the new dataLim
        ax1.autoscale_view()            

        ax1.set_title('fragment size')
        ax1.legend()

        # Plot the quadratic relationship
        # fig, ax = plt.subplot(1, 2, 2)
        for j in range(0, frgT.sz):
            ax2.add_line(self.lines_es[j])
#        for j in range(0, frgT.sz):
#            plt.plot(self.laps,self.leps[j], label = frgT.names[j], linestyle = frgT.lstyle[j], color = clr)
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')
        ax2.set_title('energy')
        ax2.legend()
        ax2.relim()
        # update ax.viewLim using the new dataLim
        ax2.autoscale_view()            
        fig.show()

    def ComputeAbsoleteCanBeDeleted(self):
        self.laps = np.arange(self.lap_min, self.lap_max, self.lap_del)
        self.aps = np.power(10.0, self.laps)
        n = len(self.aps)
        n = len(self.aps)
        self.lsps = np.zeros([frgT.sz, n])
        self.leps = np.zeros([frgT.sz, n])
        for i in range(0, n):
            ap = self.aps[i]
            fd = Frag1D1la()
            fd.ap = ap
            fd.ComputeVals()
            for j in range(0, frgT.sz):
                self.lsps[j][i] = np.log10(fd.s_ps[j])
                self.leps[j][i] = np.log10(fd.e_ps[j])
        


        # Create two plots
        plt.figure(figsize=(10, 5))
        clr = 'darkgray'

        # Plot the linear relationship
        plt.subplot(1, 2, 1)
        for j in range(0, frgT.sz):
            plt.plot(self.laps,self.lsps[j], label = frgT.names[j], linestyle = frgT.lstyle[j], color = clr)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('fragment size')
        plt.legend()

        # Plot the quadratic relationship
        plt.subplot(1, 2, 2)
        for j in range(0, frgT.sz):
            plt.plot(self.laps,self.leps[j], label = frgT.names[j], linestyle = frgT.lstyle[j], color = clr)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('energy')
        plt.legend()
        plt.show()

