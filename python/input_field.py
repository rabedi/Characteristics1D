import numpy as np
import pandas as pd
#import pickle
import glob
import math
import os
import statistics

class StatOfVec:
    def __init__(self, valsIn = []):
            if (len(valsIn) > 0):
                self.Compute_Vec_Stat(valsIn)

    def Compute_Vec_Stat(self, valsIn):
        self.num_nan = np.count_nonzero(np.isnan(valsIn))
        self.total_count = len(valsIn)
        self.count = self.total_count - self.num_nan
        if (self.count > 0):
            self.mnV = np.nanmin(valsIn)
            self.mxV = np.nanmax(valsIn)
            self.meanV = np.nanmean(valsIn)
            if (self.count > 2):
                self.std = np.nanstd(valsIn)
                self.cov = 0.0
                if (abs(self.meanV) > 1e-15):
                    self.cov = self.std / self.meanV
                self.span = self.mxV - self.mnV
            else:
                self.std = np.nan
                self.cov = np.nan
                self.span = np.nan
        else:
            self.mnV = np.nan
            self.mxV = np.nan
            self.meanV = np.nan
            self.std = np.nan
            self.cov = np.nan
            self.span = np.nan

    def GetVecVals(self, valsIn = [], retVafter_std = False):
        vals = valsIn
        vals.append(self.meanV)
        vals.append(self.mnV)
        vals.append(self.mxV)
        vals.append(self.cov)
        vals.append(self.std)
        if (retVafter_std):
            return vals
        vals.append(self.span)
        return vals

    @staticmethod
    def Get_Names(pre_name = "inp_sx_", namesIn = [], retVafter_std = False):
        names = namesIn
        names.append(pre_name + "mean")
        names.append(pre_name + "mn")
        names.append(pre_name + "mx")
        names.append(pre_name + "cov")
        names.append(pre_name + "std")
        if (retVafter_std):
            return names
        names.append(pre_name + "span")
        return names

def valReduction(valsIn, meshp2In, valsAtVert, meshp2Out = -1, reduction_sso = 3, isPeriodic = True):
    if (meshp2Out < 0):
        return valsIn
    redP2 = meshp2In - meshp2Out
    if (redP2 == 0):
        return valsIn
    step = (int)(2**redP2)
    resOut = (int)(2**meshp2Out)
    szOut = (int)(resOut) + (int)(valsAtVert)
    bkVal = valsIn[-1]
    valsOut = np.zeros(szOut)
    for i in range(resOut):
        vec = np.zeros(step)
        st = i * step
        for j in range(step):
            vec[j] = valsIn[st + j]
        if (reduction_sso == 3): # minimum val
            valsOut[i] = np.min(vec)
        elif (reduction_sso == 0): # start value
            valsOut[i] = vec[0]
        elif (reduction_sso == 1): # end value
            valsOut[i] = vec[-1]
        elif (reduction_sso == 6): # mean arithmetic
            valsOut[i] = np.mean(vec)
        elif (reduction_sso == 10): # harmonic average
            valsOut[i] = statistics.harmonic_mean(vec)
        elif (reduction_sso == 5): # maximum value
            valsOut[i] = np.max(vec)

    if valsAtVert:
        if isPeriodic:
            valsOut[-1] = valsOut[0]
        else:
            valsOut[-1] = bkVal
    return valsOut

class InpF:
    rootIFFOlder = "../../InhomogeneousFiles/"
    def __init__(self):
        pass
    @staticmethod
    def setInputMeshRootFolder(rootFolderIn = "../../InhomogeneousFiles/"):
        rootIFFOlder = rootFolderIn
    # dd2 is the "span" parameter of pointwise PDF
    # sso_valStart = 0, sso_valEnd = 1, sso_min = 3, sso_max = 5, sso_mean_arithmetic = 6, sso_mean_harmonic = 10
    def Initialize_InpF(self, valsAtVert = True, meshp2 = 14, serNo = 0, llc = -4, dd2 = 0.7, isPeriodic = True, reduction_sso = 3, meshp2_4Simulation = -1, meshp2_4Output = -1, shape = 2):
        nSeg = 2**meshp2
        if (meshp2_4Simulation == -1):
            meshp2_4Simulation = meshp2
        if (meshp2_4Output == -1):
            meshp2_4Output = meshp2
        nVert = nSeg + 1
        nVals = nSeg
        if (valsAtVert):
            nVals += 1
        self.vals = np.ones(nVals)
        if (dd2 > 1e-4):
            str_llc = "z"
            if (llc > -4.1):
                str_llc = "{:.1f}".format(llc)
            if (meshp2 == 14):
                fn = InpF.rootIFFOlder + "/cl" + str_llc + "_np16385/initial_values_" + str(int(serNo)) + ".txt"
            self.vals = np.loadtxt(fn)[1:]
            minV = 1 - dd2
            maxV = 1 + dd2
            inv_sqrt2 = 1.0 / np.sqrt(2.0)
            # periodic domain for ring problem
            if (valsAtVert and isPeriodic):
                vm = 0.5 * (self.vals[0] + self.vals[-1])
                self.vals[0] = vm
                self.vals[-1] = vm
            
            def SN_2_genTri(xstandardNormalValue): # x is a sandard normal value is is going to be mapped to a triangular one with mean 1 and min/max 1 - dd2, 1 + dd2
                p = 0.5 * (1.0 + math.erf(xstandardNormalValue * inv_sqrt2)) # p = sn_cdf
                if (shape == 2):
                    if (p < 0.5):
                        return 1.0 + dd2 * (-1.0 + np.sqrt(2.0 * p)) 
                    return 1 + dd2 * (1.0 - np.sqrt(2.0 * (1.0 - p)))
                if (shape == 1):
                    return 1 + dd2 * (2.0 * p - 1.0)
                inv_shape = 1.0 / shape
                if (p < 0.5):
                    return 1.0 + dd2 * (-1.0 + pow(2.0 * p, inv_shape)) 
                return 1 + dd2 * (1.0 - pow(2.0 * (1.0 - p), inv_shape))

            self.vals = list(map(SN_2_genTri, self.vals))

            # this is the actual input mesh (with potentially reduced resolution)
            # computing stats
            self.vals4Simulation = valReduction(self.vals, meshp2, valsAtVert, meshp2_4Simulation, reduction_sso, isPeriodic)
            self.stats4SimulationFld = StatOfVec(self.vals4Simulation)
            if (meshp2 == meshp2_4Simulation):
                self.stats4RawInputFld = self.stats4SimulationFld
            else:
                self.stats4RawInputFld = StatOfVec(self.vals)
    
            # computing values for output
            # the following is not correct as the final output is derived from the field for stat (that is the actual input mesh)
            # self.vals4Output = valReduction(self.vals, meshp2, valsAtVert, meshp2_4Stat, reduction_sso, meshp2_4Output)
            if ((meshp2_4Output == meshp2) or (meshp2_4Output > meshp2_4Simulation)):
                self.vals4Output = self.vals
            else:
                reduction_sso_output = 1 # start value
                self.vals4Output = valReduction(self.vals4Simulation, meshp2_4Simulation, valsAtVert, meshp2_4Output, reduction_sso_output, isPeriodic)
        else:
            self.stats4SimulationFld = StatOfVec(self.vals)
            self.stats4RawInputFld = self.stats4SimulationFld

            nPts4o = 2**meshp2_4Output
            if (valsAtVert):
                nPts4o += 1
            self.vals4Output = np.ones(nPts4o)

    def GetStatVecVals(self, addRaw = True):
        vals = []
        self.stats4SimulationFld.GetVecVals(vals)
        if (addRaw):
            self.stats4RawInputFld.GetVecVals(vals)
        return vals

    @staticmethod
    def GetStatNames(pre_name = "inp_sx_", addRaw = True):
        names = []
        StatOfVec.Get_Names(pre_name + "sim_", names)
        if (addRaw):
            StatOfVec.Get_Names(pre_name + "raw_", names)
        return names
