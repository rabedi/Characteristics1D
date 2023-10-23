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

from rutil import pwl_root_crossing_segmentStats

class StatOfVec:
    sv_calculate_cor_fun = 1
    sv_calculate_HD = 1
    sv_calculate_crossing = 1
    def __init__(self):
        self.num_nan = 0
        self.total_count = 0
        self.count = 0
        self.mnV = np.nan
        self.mxV = np.nan
        self.meanV = np.nan
        self.std =  np.nan
        self.cov =  np.nan
        self.span = np.nan
        # H & D
        self.hfD = np.nan
        self.HVal = np.nan
        self.Hc = np.nan

        # segment (crossing from 0) stats
        self.seg_num = np.nan
        self.seg_meanV = np.nan
        self.seg_mnV = np.nan
        self.seg_mxV = np.nan
        self.seg_std = np.nan

        # cov function
        self.ac_ys = []
        self.ac_l1intAll = np.nan
        self.ac_l2intAll = np.nan
        self.ac_l1intRed = np.nan
        self.ac_l2intRed = np.nan
        self.ac_xCorNeg = np.nan

    def Compute_Vec_Stat(self, valsIn, compute_fieldStats):
        self.num_nan = np.count_nonzero(np.isnan(valsIn))
        self.total_count = len(valsIn)
        self.count = self.total_count - self.num_nan
        if (self.count > 0):
            self.mnV = np.nanmin(valsIn)
            self.mxV = np.nanmax(valsIn)
            self.meanV = np.nanmean(valsIn)
            if (self.count > 2):
                self.std = np.nanstd(valsIn)
                sz = self.count
                dl = 1 / (sz - 1)
                xs = np.arange(0, 1 + dl, dl)
                if (abs(self.std) > 0):
                    valsN = (valsIn - self.meanV)/ self.std
                else:
                    valsN = (valsIn - self.meanV)

                self.cov = 0.0
                valsEqual = False
                absMean = abs(self.meanV)
                if (absMean > 1e-15):
                    self.cov = self.std / absMean
                    valsEqual = (self.cov < 1e-4)
                else:
                    mx = np.max(np.abs(valsIn))
                    if (mx < 1e-10):
                        valsEqual = True
                    else:
                        valsEqual = (self.std < 1e-4 * mx)
                self.span = self.mxV - self.mnV
                if ((valsEqual == False) and (compute_fieldStats == True)): 
                    if ((StatOfVec.sv_calculate_HD != 0) and (len(valsIn) > 100)): 
                        self.hfD = hfd.hfd(valsIn)
                        # Evaluate Hurst equation
                        self.HVal, self.Hc, data = compute_Hc(valsIn, kind='change', simplified=True)
                    if (StatOfVec.sv_calculate_crossing != 0):
                        self.seg_num, self.seg_meanV, self.seg_mnV, self.seg_mxV, self.seg_std = pwl_root_crossing_segmentStats(xs, valsN)
                    if (StatOfVec.sv_calculate_cor_fun != 0):
                        autocorr_ret = signal.fftconvolve(valsN, valsN[::-1], mode='full')
                        self.ac_ys = autocorr_ret[sz - 1:2 * sz]
                        self.ac_ys = self.ac_ys / self.ac_ys[0]
                        sz2 = len(self.ac_ys)
                        dl = 1 / (sz2 - 1)
                        ac_x = np.arange(0, 1 + dl, dl)
                        first_negative_index = np.argmax(self.ac_ys < 0)
                        ac2abs = abs(self.ac_ys)
                        ac2abs2 = ac2abs**2
                        self.ac_l1intAll = np.trapz(ac2abs, ac_x)
                        self.ac_l2intAll = np.trapz(ac2abs2, ac_x)
                        self.ac_l1intRed = self.ac_l1intAll
                        self.ac_l2intRed = self.ac_l2intAll
                        if (first_negative_index > 0):
                            self.ac_xCorNeg = ac_x[first_negative_index]
                            xRed = ac_x[0:first_negative_index]
                            acRed = self.ac_ys[0:first_negative_index]
                            acRed2abs = abs(acRed)
                            acRed2abs2 = acRed2abs**2
                            self.ac_l1intRed = np.trapz(acRed2abs, xRed)
                            self.ac_l2intRed = np.trapz(acRed2abs2, xRed)

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
        if (StatOfVec.sv_calculate_HD):
            vals.append(self.hfD)
            vals.append(self.HVal)
            vals.append(self.Hc)
        if (StatOfVec.sv_calculate_cor_fun):
            vals.append(self.ac_xCorNeg)
            vals.append(self.ac_l1intAll)
            vals.append(self.ac_l2intAll)
            vals.append(self.ac_l1intRed)
            vals.append(self.ac_l2intRed)
        if (StatOfVec.sv_calculate_crossing):
            vals.append(self.seg_num)
            vals.append(self.seg_meanV)
            vals.append(self.seg_mnV)
            vals.append(self.seg_mxV)
            vals.append(self.seg_std)
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

        if (StatOfVec.sv_calculate_HD):
            names.append(pre_name + "hfD")
            names.append(pre_name + "HVal")
            names.append(pre_name + "Hc")
        if (StatOfVec.sv_calculate_cor_fun):
            names.append(pre_name + "ac_xCorNeg")
            names.append(pre_name + "ac_l1intAll")
            names.append(pre_name + "ac_l2intAll")
            names.append(pre_name + "ac_l1intRed")
            names.append(pre_name + "ac_l2intRed")
        if (StatOfVec.sv_calculate_crossing):
            names.append(pre_name + "seg_num")
            names.append(pre_name + "seg_meanV")
            names.append(pre_name + "seg_mnV")
            names.append(pre_name + "seg_mxV")
            names.append(pre_name + "seg_std")
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
        self.valid = True
    @staticmethod
    def setInputMeshRootFolder(rootFolderIn = "../../InhomogeneousFiles/"):
        rootIFFOlder = rootFolderIn
    # dd2 is the "span" parameter of pointwise PDF
    # sso_valStart = 0, sso_valEnd = 1, sso_min = 3, sso_max = 5, sso_mean_arithmetic = 6, sso_mean_harmonic = 10
    def Initialize_InpF(self, valsAtVert = True, meshp2 = 14, serNo = 0, llc = -3, dd2 = 0.7, isPeriodic = True, reduction_sso = 3, meshp2_4Simulation = -1, meshp2_4Output = -1, shape = 2, useOriginalSN_Mesh4Raw = 0):
        self.valid = True
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

            path = Path(fn)
            if not path.exists():
                self.valid = False
                return
            self.vals = np.loadtxt(fn)[1:]
            minV = 1 - dd2
            maxV = 1 + dd2
            inv_sqrt2 = 1.0 / np.sqrt(2.0)
            # periodic domain for ring problem
            if (valsAtVert and isPeriodic):
                vm = 0.5 * (self.vals[0] + self.vals[-1])
                self.vals[0] = vm
                self.vals[-1] = vm
            self.valsBK = self.vals
            
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

            if (shape >= -0.5): # a transformation from SN to a different field occurs
                self.vals = list(map(SN_2_genTri, self.vals))

                # this is the actual input mesh (with potentially reduced resolution)
                # computing stats
                self.vals4Simulation = valReduction(self.vals, meshp2, valsAtVert, meshp2_4Simulation, reduction_sso, isPeriodic)
                self.stats4SimulationFld = StatOfVec()
                self.stats4SimulationFld.Compute_Vec_Stat(self.vals4Simulation, True)
                if (useOriginalSN_Mesh4Raw == 1):
                    self.stats4RawInputFld = StatOfVec()
                    self.stats4RawInputFld.Compute_Vec_Stat(self.valsBK, True)
                else:
                    if ((meshp2 == meshp2_4Simulation) or (useOriginalSN_Mesh4Raw == 2)):
                        self.stats4RawInputFld = self.stats4SimulationFld
                    else:
                        self.stats4RawInputFld = StatOfVec()
                        self.stats4RawInputFld.Compute_Vec_Stat(self.vals, True)
            else:
                self.stats4SimulationFld = StatOfVec()
                self.stats4SimulationFld.Compute_Vec_Stat(self.vals, True)
                self.stats4RawInputFld = self.stats4SimulationFld
    
            # computing values for output
            # the following is not correct as the final output is derived from the field for stat (that is the actual input mesh)
            # self.vals4Output = valReduction(self.vals, meshp2, valsAtVert, meshp2_4Stat, reduction_sso, meshp2_4Output)
            if ((meshp2_4Output == meshp2) or (meshp2_4Output > meshp2_4Simulation)):
                self.vals4Output = self.vals
            else:
                reduction_sso_output = 1 # start value
                self.vals4Output = valReduction(self.vals4Simulation, meshp2_4Simulation, valsAtVert, meshp2_4Output, reduction_sso_output, isPeriodic)
        else:
            self.valsBK = 0 * self.vals
            self.stats4SimulationFld = StatOfVec()
            self.stats4SimulationFld.Compute_Vec_Stat(self.vals, True)

            if (useOriginalSN_Mesh4Raw):
                self.stats4RawInputFld = StatOfVec()
                self.stats4RawInputFld.Compute_Vec_Stat(0 * self.vals, True)
            else:
                self.stats4RawInputFld = self.stats4SimulationFld

            nPts4o = 2**meshp2_4Output
            if (valsAtVert):
                nPts4o += 1
            self.vals4Output = np.ones(nPts4o)

    def GetStatVecVals(self, addRaw = True, addMapped = True):
        vals = []
        if (addMapped):
            self.stats4SimulationFld.GetVecVals(vals)
        if (addRaw):
            self.stats4RawInputFld.GetVecVals(vals)
        return vals

    @staticmethod
    def GetStatNames(pre_nameMapped = "inp_sx_sim_", pre_nameRaw = "inp_sx_raw_", addMapped = True, addRaw = True): 
        names = []
        if (addMapped):
            StatOfVec.Get_Names(pre_nameMapped, names)
        if (addRaw):
            StatOfVec.Get_Names(pre_nameRaw, names)
        return names


# class for outputting all data for one shape, one Delta, on corr length, ALL serial numbers
class InpFsOuput:
    rootIFFOlder = "../../InhomogeneousFiles/"
    def __init__(self):
        pass
    @staticmethod
    def setInputMeshRootFolder(rootFolderIn = "../../InhomogeneousFiles/"):
        rootIFFOlder = rootFolderIn
        InpF.setInputMeshRootFolder(rootFolderIn)

    # dd2 is the "span" parameter of pointwise PDF
    # sso_valStart = 0, sso_valEnd = 1, sso_min = 3, sso_max = 5, sso_mean_arithmetic = 6, sso_mean_harmonic = 10
    # step4Cov, how many points to skip for the function y values, 1: all taken, 2: every other, 3: one in every three ...
    def Print(self, serNos = range(0,10), shapes = [1,2, 4, 10], llcs = [-3], dd2s = [0.0], fileNameBase = "fnb", includeStatsOfSN = True, includeStatsOfMappedFld = True, printCOVFunction = True, step4Cov = 1, printFieldVals = True, stepFieldVals = 1, isPeriodic = True, valsAtVert = True):
        meshp2 = 14
        useOriginalSN_Mesh4Raw = 1
        reduction_sso = 3
        meshp2_4Simulation = -1
        meshp2_4Output = -1

        pre_nameMapped = "map_"
        pre_nameRaw = "sn_"
        addRaw = True
        addMapped = True

        if any(x < 0 for x in shapes):
            includeStatsOfMappedFld = False

        if (includeStatsOfSN == False):
            addRaw = False
            useOriginalSN_Mesh4Raw = 2 # this inactivates the call to calculate the stat of the original fields
            if (includeStatsOfMappedFld == False):
                print("This is an invalid option as original and mapped field are both turned off\n")
        else:
            if (includeStatsOfMappedFld == False):
                # a trick to turn off mapped field is to use negative shape
                shapes = [-1]
                llcs = [-4]
                dd2s = [0.7] 
                addMapped = False


        fnStat = fileNameBase + "_stat.csv"
        fstat = open(fnStat, 'w')
        fnCovSN = fileNameBase + "_cov_SN.csv"
        fnCovMap = fileNameBase + "_cov_mapped.csv"
        cov_includeStatsOfSN = printCOVFunction and includeStatsOfSN
        if (cov_includeStatsOfSN):
            fCovSN = open(fnCovSN, 'w')
            fCovSN.write('shape,llc,dd2,serNo,sz,vals\n')
        cov_includeStatsOfMappedFld = printCOVFunction and includeStatsOfMappedFld
        if (cov_includeStatsOfMappedFld):
            fCovMap = open(fnCovMap, 'w')
            fCovMap.write('shape,llc,dd2,serNo,sz,vals\n')

        fnfvSN = fileNameBase + "_fieldVal_SN.csv"
        fnfvMap = fileNameBase + "_fieldVal_mapped.csv"
        fv_includeStatsOfSN = printFieldVals and includeStatsOfSN
        if (fv_includeStatsOfSN):
            ffvSN = open(fnfvSN, 'w')
            ffvSN.write('shape,llc,dd2,serNo,sz,vals\n')
        fv_includeStatsOfMappedFld = printFieldVals and includeStatsOfMappedFld
        if (fv_includeStatsOfMappedFld):
            ffvMap = open(fnfvMap, 'w')
            ffvMap.write('shape,llc,dd2,serNo,sz,vals\n')

        names = InpF.GetStatNames(pre_nameMapped, pre_nameRaw, addMapped, addRaw) 
        sz_stat_fields = len(names)
        comma_separated_data = ",".join(names)
        fstat.write('shape,llc,dd2,serNo,')
        fstat.write(comma_separated_data)
        fstat.write('\n')

        for shape in shapes:
            for llc in llcs:
                for dd2 in dd2s:
                    for serNo in serNos:
                        pinp = InpF()
                        pinp.Initialize_InpF(valsAtVert, meshp2, serNo, llc, dd2, \
                                            isPeriodic, reduction_sso, \
                                            meshp2_4Simulation, meshp2_4Output, shape, useOriginalSN_Mesh4Raw)
                        if (self.valid == False):
                            continue
                        if (serNo % 100 == 0):
                            print(f"serNo = {serNo}")
                        statVals = pinp.GetStatVecVals(addRaw, addMapped)
                        str_lab = str(shape) + ',' + str(llc) + ',' + str(dd2) + ',' + str(serNo)
                        fstat.write(str_lab)
                        for vl in statVals:
                            fstat.write(f",{vl}")
#                        comma_separated_data = ",".join(map(str, statVals))
#                        fstat.write(comma_separated_data)
                        fstat.write('\n')

                        if (cov_includeStatsOfSN):
                            fCovSN.write(str_lab)
                            vc = pinp.stats4SimulationFld.ac_ys
                            if (step4Cov > 1):
                                vc = vc[::step4Cov]
                            sz = len(vc)
                            fCovSN.write(f",{sz}")
                            for vl in vc:
                                fCovSN.write(f",{vl}")
                            fCovSN.write('\n')
                        if (cov_includeStatsOfMappedFld):
                            fCovMap.write(str_lab)
                            vc = pinp.stats4SimulationFld.ac_ys
                            if (step4Cov > 1):
                                vc = vc[::step4Cov]
                            sz = len(vc)
                            fCovMap.write(f",{sz}")
                            for vl in vc:
                                fCovMap.write(f",{vl}")
                            fCovMap.write('\n')

                        if (fv_includeStatsOfSN):
                            ffvSN.write(str_lab)
                            vc = pinp.valsBK
                            if (stepFieldVals > 1):
                                vc = vc[::stepFieldVals]
                            sz = len(vc)
                            ffvSN.write(f",{sz}")
                            for vl in vc:
                                ffvSN.write(f",{vl}")
                            ffvSN.write('\n')
                        if (fv_includeStatsOfMappedFld):
                            ffvMap.write(str_lab)
                            vc = pinp.vals
                            if (stepFieldVals > 1):
                                vc = vc[::stepFieldVals]
                            sz = len(vc)
                            ffvMap.write(f",{sz}")
                            for vl in vc:
                                ffvMap.write(f",{vl}")
                            ffvMap.write('\n')

        fstat.close()
        if (includeStatsOfSN):
            fCovSN.close()
        if (includeStatsOfMappedFld):
            fCovMap.close()

        if (fv_includeStatsOfSN):
            ffvSN.close()
        if (fv_includeStatsOfMappedFld):
            ffvMap.close()

    def Print_AllSerials_OnePara(self, serNoMax = 100, shape = 2, llc = -3, dd2 = 0.7, includeStatsOfSN = True, includeStatsOfMappedFld = True, printCOVFunction = True, step4Cov = 1, printFieldVals = True, stepFieldVals = 1, isPeriodic = True, valsAtVert = True):
        serNos = range(0,serNoMax + 1)
        shapes = [shape]
        llcs = [llc]
        dd2s = [dd2]
        fileNameBase = "fnb_shape_" + str(shape) + "_llc_" + str(llc) + "_dd2_" + str(dd2)
        self.Print(serNos, shapes, llcs, dd2s, fileNameBase, includeStatsOfSN, includeStatsOfMappedFld, printCOVFunction, step4Cov, printFieldVals, stepFieldVals, isPeriodic, valsAtVert)


    def Print_AllSerials_AllPara(self, serNoMax = 1300, shapes = [1, 1.5, 2, 3, 4], llcs = [-2], dd2s = [0.5], includeStatsOfSN = True, includeStatsOfMappedFld = True, printCOVFunction = True, step4Cov = 1, printFieldVals = True, stepFieldVals = 1, isPeriodic = True, valsAtVert = True):

        for shape in shapes:
            print(f"shape = {shape}\n")
            for llc in llcs:
                print(f"llc = {llc}\n")
                for dd2 in dd2s:
                    print(f"dd2 = {dd2}\n")
                    self.Print_AllSerials_OnePara(serNoMax, shape, llc, dd2, includeStatsOfSN, includeStatsOfMappedFld, printCOVFunction, step4Cov, printFieldVals, stepFieldVals, isPeriodic, valsAtVert)
