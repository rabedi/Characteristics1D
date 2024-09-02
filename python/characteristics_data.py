import numpy as np
import pandas as pd
import glob
import os
import shutil
import rutil as rutil
import math
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pickle
import rmulti_index as rmi 
from rmulti_index import SMIndex as SMIndex
from input_field import StatOfVec as StatOfVec
from input_field import InpF as InpF
import time

# mainOptions
# 0, 1: effect of llc, dd2 ldelc, la
#       0: coarse mesh for la = -3, ...
#       1: 2^14 for all meshes

# 2: resolution_x for fracture (SVE homogenization - Andy Tonge discussion)
# 3 similar but with coarsening columns
mo_frac_rate_cf = 0
mo_frac_rate_f  = 1
mo_frac_rate_f_wShape = 2
mo_frac_res_x = 3 # coarsening data for which the energies don't match for high loading rate - see 3 below
mo_frac_res_x_w_delc_fact = 4 # new data 7/23 that includes factoring deltaC to get the correct energy

mainOption = mo_frac_rate_f_wShape # mo_frac_rate_f_wShape #mo_frac_rate_f
rateStudy = ((mainOption == mo_frac_rate_cf) or (mainOption == mo_frac_rate_f) or (mainOption == mo_frac_rate_f_wShape))


# Change
if (rateStudy):
    if (mainOption == mo_frac_rate_f_wShape):
        sortingFields = ["shape", "llc", "dd2", "ldelc", "la"]
    else:
        sortingFields = ["llc", "dd2", "ldelc", "la"]
elif (mainOption == mo_frac_res_x):
    sortingFields = ["dt_resap", "ssoFS", "la", "llc", "resolutionFactor"]
elif (mainOption == mo_frac_res_x_w_delc_fact):
    sortingFields = ["delc_Tf_fact", "dt_resap", "ssoFS", "la", "llc", "resolutionFactor"]

# spatial field data
valsAtVert = True
meshp2 = 14
isPeriodic = True
reduction_sso = 3
meshp2_4Output = -1
meshp2_4Simulation = -1
addRawStat = True

# take the log of these fields and adds them to all
dist_logs = {}

dist_logs["time_F"] = 10
dist_logs["time_M"] = 10
dist_logs["time_sigma0"] = 10
dist_logs["time_norm_tauc_F"] = 10
dist_logs["time_norm_tauc_M"] = 10
dist_logs["time_sigma0_norm_tauc"] = 10
dist_logs["time_norm_taul_F"] = 10
dist_logs["time_norm_taul_M"] = 10
dist_logs["time_sigma0_norm_taul"] = 10
dist_logs["phi_d_tFin"] = 10
dist_logs["phi_d_norm_phi0_tFin"] = 10
dist_logs["phi_d_norm_phi_input_F"] = 10
dist_logs["psi_f"] = 10
dist_logs["psi_f_norm_psiC"] = 10
dist_logs["eps_f"] = 10
dist_logs["eps_f_norm"] = 10
dist_logs["lbar_F"] = 10
dist_logs["lbar_tF"] = 10
dist_logs["lbar_DelU_tF"] = 10
dist_logs["lbar_M"] = 10
dist_logs["lbar_DelU_M"] = 10

dist_logs["lbar_F_orig"] = 10
dist_logs["lbar_tF_orig"] = 10
dist_logs["lbar_DelU_tF_orig"] = 10
dist_logs["lbar_M_orig"] = 10
dist_logs["lbar_DelU_M_orig"] = 10


def read_updata_csv(filename = "../StochasticPostprocessor/data/2023_03_20/pps3Out_ub1.csv"):
    if (rateStudy):
        hasShape = (mainOption == mo_frac_rate_f_wShape)
        pd = read_updata_csv_frac_rate(filename, hasShape)
    elif ((mainOption == mo_frac_res_x) or (mainOption == mo_frac_res_x_w_delc_fact)):
        has_deltaCFactor = (mainOption == mo_frac_res_x_w_delc_fact)
        pd = read_updata_csv_frac_res_x(filename, has_deltaCFactor)
    nondimensionLFactor = 2
    fragmentSizeNorm_wrt_lcoh = 1
    pd = Update4NondimensionalVals(pd, nondimensionLFactor, fragmentSizeNorm_wrt_lcoh)
    return pd

def Update4NondimensionalVals(db, nondimensionLFactor = 1, fragmentSizeNorm_wrt_lcoh = 0):
    nondimensionLFactorInv = 1.0 / nondimensionLFactor
    lf = np.log(nondimensionLFactor) / np.log(10)
    c_la_v = "inp_s_la"
    i_la_v = db.columns.get_loc(c_la_v)
    c_lap_v = "inp_s_lap"
    i_lap_v = -1
    has_lap = 0
    if (c_lap_v in db.columns):
        i_lap_v = db.columns.get_loc(c_lap_v)
        has_lap = 1
    else:
        i_lap_v = i_la_v + 1
        db.insert(i_lap_v,c_lap_v, 0)

    c_llc_v = "inp_s_llc"
    i_llc_v = db.columns.get_loc(c_llc_v)
    c_llcp_v = "inp_s_llcp"
    i_llcp_v = -1
    has_llcp = 0
    if (c_llcp_v in db.columns):
        i_llcp_v = db.columns.get_loc(c_llcp_v)
        has_llcp = 1
    else:
        i_llcp_v = i_llc_v + 1
        db.insert(i_llcp_v,c_llcp_v, 0)

    colsNormalized_t_tauc = []
    if (nondimensionLFactor != 1):
        for j, coln in enumerate(db.columns):
            if '_norm_tauc' in coln:
                colsNormalized_t_tauc.append(j)

    # figuring out which columns are related to fragments
    cfrag_N = []
    cfrag_len = []
    cfrag_lbar = []
    if (fragmentSizeNorm_wrt_lcoh == 1):
        for j, coln in enumerate(db.columns):
            if '_Nl_' in coln:
                cfrag_N.append(j)
                nmNew = coln + '_orig'
                db[nmNew] = db[coln]
            elif ('_lbar_' in coln): 
                cfrag_len.append(j)
                cfrag_lbar.append(j)
                nmNew = coln + '_orig'
                db[nmNew] = db[coln]
            elif (('_lmin_' in coln) or ('_lmax_' in coln) or ('_sdiv_l_' in coln)):
                cfrag_len.append(j)

    c_lldelc_v = "inp_s_ldelc"
    i_lldelc_v = db.columns.get_loc(c_lldelc_v)

    i = 0
    for row in db.iterrows():
        la = db.iloc[i, i_la_v]
        ldelc = db.iloc[i, i_lldelc_v]
        delc = np.power(10.0, ldelc)
        lTilde = nondimensionLFactorInv * delc
        ldelcf = ldelc - lf
        if not has_lap:
            lap = la + ldelcf
            db.iloc[i, i_lap_v] = lap

        if not has_llcp:
            llc = db.iloc[i, i_llc_v]
            llcp = llc - ldelcf
            db.iloc[i, i_llcp_v] = llcp

        for j in colsNormalized_t_tauc:
            db.iloc[i, j] *= nondimensionLFactor

        if (fragmentSizeNorm_wrt_lcoh == 1):
            lTildeInv = 1.0 / lTilde
            for j in cfrag_N:
                db.iloc[i, j] *= lTilde
            for j in cfrag_len:
                db.iloc[i, j] *= lTildeInv
        i += 1
    return db

def read_updata_csv_frac_rate(filename = "../StochasticPostprocessor/data/2023_03_20/pps3Out_ub1.csv", hasShape = 0):
    db = pd.read_csv(filename)
    db = db.apply(pd.to_numeric, errors='coerce')
    #change
    c_ver0 = "verNo"
    c_ver1 = "verWOffsetNo"
    c_lldelc_ind = "indldelc"
    c_lldelc_v = "inp_s_ldelc"
    c_la_ind = "indla"
    c_la_v = "inp_s_la"
    c_dd2_ind = "inddd2"
    c_dd2_v = "inp_s_dd2"
    c_llc_ind = "indllc"
    c_llc_v = "inp_s_llc"
#    c_shape_ind = "indshape"
    c_shape_v = "inp_s_shape"

#    i_shape_ind = -1
    i_shape_v = -1

#    if (not hasShape):
#        index_position = db.columns.get_loc("inddd2") + 2
#        column_name = c_shape_ind
#        db.insert(index_position, column_name, 2)
#        index_position = index_position + 1
#        column_name = c_shape_v
#        db.insert(index_position, column_name, 2.0)
        #		1.0 #		1.5 #	2.0 #		3.0 #		4.0

    if (hasShape):
#        i_shape_ind = db.columns.get_loc(c_llc_ind)
        i_shape_v = db.columns.get_loc(c_shape_v)

    i_ver0 = db.columns.get_loc(c_ver0)
    i_ver1 = db.columns.get_loc(c_ver1)
    i_ldelc_ind = db.columns.get_loc(c_lldelc_ind)
    i_ldelc_v = db.columns.get_loc(c_lldelc_v)
    i_la_ind = db.columns.get_loc(c_la_ind)
    i_la_v = db.columns.get_loc(c_la_v)
    i_dd2_ind = db.columns.get_loc(c_dd2_ind)
    i_dd2_v = db.columns.get_loc(c_dd2_v)
    i_llc_ind = db.columns.get_loc(c_llc_ind)
    i_llc_v = db.columns.get_loc(c_llc_v)
    column_name = 'out_s_time_sigma0'
    i_col_ts0 = -1
    if column_name in db.columns:
        i_col_ts0 = db.columns.get_loc(column_name)    

    i_col_phidtf = db.columns.get_loc("out_s_phi_d_tFin")
    log10Inv = 1.0 / np.log(10)

    nRows = db.shape[0]
    ii = 0
    for i, row in db.iterrows():
        ver0 = int(row[i_ver0])
        ver1 = int(row[i_ver1])
        la_ind = int(row[i_la_ind])
        la_v = row[i_la_v]
        ldelc_v = row[i_ldelc_v]
        dd2_ind = int(row[i_dd2_ind])
        dd2_v = row[i_dd2_v]
        llc_v = row[i_llc_v]
 
        phidtf = row[i_col_phidtf]
        if (phidtf < 0):
            db = db.drop(i, axis=0)
            continue
        timeSigma0 = 1.0
        if (i_col_ts0 > 0):
            timeSigma0 = row[i_col_ts0]
        if((math.isnan(timeSigma0)) or (timeSigma0 <= 0)):
            db = db.drop(i, axis=0)
            continue
        logphid = np.log(phidtf) * log10Inv
        if (la_v > -1.1):
            logphidMin = -0.45 + la_v * 0.6125
        else:
            logphidMin = -1.5
        logphidMin = logphidMin + (ldelc_v + 1)
        if (la_v < -1.1):
            logphidMin -= (dd2_v - 0.1) 

        # decrementing by 1 so to drop bad data
        logphidMin = logphidMin - 2
        if (logphid < logphidMin):
            db = db.drop(i, axis=0)
            continue

        ver0 = min(ver0 % 2000, ver1 % 2000)
        ver1 = ver0

        if (hasShape):
            dd2_ind = 0
        else:
            if abs(dd2_v) < 0.001:
                ver0 += 5000
                dd2_ind = 0
            elif (abs(dd2_v - 0.1) < 0.0001):
                dd2_ind = 1
            elif (abs(dd2_v - 0.3) < 0.0001):
                dd2_ind = 2
            elif (abs(dd2_v - 0.5) < 0.0001):
                dd2_ind = 3
            elif (abs(dd2_v - 0.7) < 0.0001):
                dd2_ind = 4
            elif (abs(dd2_v - 0.9) < 0.0001):
                dd2_ind = 5

        if (la_v < 0):
            ver0 += 10000
        la_ind = int(round(2 * la_v)) + 6         
        if (hasShape):
            shape_v = row[i_shape_v]
            ver0 = ver0 + math.floor(100000 * shape_v)
        ver1 = ver0

        db.iloc[ii, i_ver0] = ver0 
        db.iloc[ii, i_ver1] = ver1 

        db.iloc[ii, i_la_ind] = la_ind
        db.iloc[ii, i_la_v] = la_v
         
        db.iloc[ii, i_dd2_ind] = dd2_ind
        db.iloc[ii, i_dd2_v] = dd2_v

        # llc_ind = int(-2 * llc_v) - 1
        # 2023/12/17: I am going to add llc_v = 0 to represent homogeneous field - llc_v = 0 (fake value for homogeneous field) turns to 0 now
        llc_ind = int(-2 * llc_v)
        db.iloc[ii, i_llc_ind] = llc_ind
        ii += 1
    nRows = db.shape[0]
    return db

def read_updata_csv_frac_res_x(filename = "../../data/resolution_x_fracture_scalars/pps3Out_resF_la2p0_004_0_9.csv", has_deltaCFactor = 0):
    goodRun = "nresF"
    zeroRes = "_001_"
    badHarmonic = not ((goodRun in filename) or (zeroRes in filename))
    alreadyHasDeltaCFactor = ("mres" in filename)

    db = pd.read_csv(filename)

    if ((not alreadyHasDeltaCFactor) and (has_deltaCFactor)):
        c_ver1 = "verWOffsetNo"
        i_ver1 = db.columns.get_loc(c_ver1)
        index_position = i_ver1 + 1
        column_name = 'inddelc_Tf_fact'
        db.insert(index_position, column_name, 0)

        index_position = i_ver1 + 2
        column_name = 'inp_s_delc_Tf_fact'
        db.insert(index_position, column_name, 1)

    colnames = db.columns
    numCols = len(colnames)
    nRows = db.shape[0]

    c_ssoFS_ind = "indssoFS"
    c_ssoFS_v = "inp_s_ssoFS"
    i_ssoFS_ind = db.columns.get_loc(c_ssoFS_ind)
    i_ssoFS_v = db.columns.get_loc(c_ssoFS_v)

    pdOut = pd.DataFrame(columns=colnames)

    db = db.apply(pd.to_numeric, errors='coerce')

    ii = 0
    for i, row in db.iterrows():
        sso_v = row[i_ssoFS_v]
#        sso_i = row[i_ssoFS_ind]
        if (sso_v == 3): #"min"):
            sso_i = 0
        elif (sso_v == 6): #"mean_arithmetic"):
            sso_i = 1
        elif (sso_v == 10): #"mean_harmonic"):
            sso_i = 2
        elif (sso_v == 0): #"valStart"):
            sso_i = 3

        db.iloc[i, i_ssoFS_ind] = sso_i 
        # db.iloc[i, i_ssoFS_v] = sso_i 

 
    c_ver0 = "verNo"
    c_ver1 = "verWOffsetNo"

    c_dt_resap_ind = "inddt_resap"
    c_dt_resap_v = "inp_s_dt_resap"

    c_la_ind = "indla"
    c_la_v = "inp_s_la"
    c_llc_ind = "indllc"
    c_llc_v = "inp_s_llc"

    c_resolutionFactor_ind = "indresolutionFactor"
    c_resolutionFactor_v = "inp_s_resolutionFactor"

    i_ver0 = db.columns.get_loc(c_ver0)
    i_ver1 = db.columns.get_loc(c_ver1)

    i_dt_resap_ind = db.columns.get_loc(c_dt_resap_ind)
    i_dt_resap_v = db.columns.get_loc(c_dt_resap_v)

    i_la_ind = db.columns.get_loc(c_la_ind)
    i_la_v = db.columns.get_loc(c_la_v)
    i_llc_ind = db.columns.get_loc(c_llc_ind)
    i_llc_v = db.columns.get_loc(c_llc_v)

    i_resolutionFactor_ind = db.columns.get_loc(c_resolutionFactor_ind)
    i_resolutionFactor_v = db.columns.get_loc(c_resolutionFactor_v)

    i_col_phidtf = db.columns.get_loc("out_s_phi_d_tFin")
    log10Inv = 1.0 / np.log(10)

    n_resaps = 2
    n_ssoFS = 4
    n_llc = 9
    n_la = 15
    n_resolutionFactor = 10

    nn_resolutionFactor = 1
    fact = n_resolutionFactor
    nn_llc = fact
    fact *= n_llc
    nn_la = fact
    fact *= n_la
    nn_ssoFS = fact
    fact *= n_ssoFS
    nn_resap = fact

    nRows = db.shape[0]
    #ii = 0
    for ii, row in db.iterrows():
        # ver0 = int(row[i_ver0])
        # ver1 = int(row[i_ver1])
        la_ind = int(row[i_la_ind])
        la_v = row[i_la_v]

        ldelc_v = -1
        dd2_ind = 0
        dd2_v = 0.9

        llc_v = row[i_llc_v]

        dt_resap_ind = int(row[i_dt_resap_ind])
        dt_resap_v = row[i_dt_resap_v]
        dt_resap_ind = dt_resap_v
        db.iloc[ii, i_dt_resap_ind] = dt_resap_ind 
        db.iloc[ii, i_dt_resap_v] = dt_resap_v 

        dt_ssoFS_ind = int(row[i_ssoFS_ind])
        dt_ssoFS_v = int(row[i_ssoFS_v])
        if ((dt_ssoFS_v == 10) and badHarmonic):
            continue 

        resolutionFactor_ind = int(row[i_resolutionFactor_ind])
        resolutionFactor_v = row[i_resolutionFactor_v]
        resolutionFactor_ind = round(np.log(resolutionFactor_v) / np.log(2))
        db.iloc[ii, i_resolutionFactor_v] = resolutionFactor_ind 
        db.iloc[ii, i_resolutionFactor_ind] = resolutionFactor_ind 


        phidtf = row[i_col_phidtf]
        if (phidtf < 0):
            # db = db.drop(i, axis=0)
            continue
        logphid = np.log(phidtf) * log10Inv
        if (la_v > -1.1):
            logphidMin = -0.45 + la_v * 0.6125
        else:
            logphidMin = -1.5
        logphidMin = logphidMin + (ldelc_v + 1)
        if (la_v < -1.1):
            logphidMin -= (dd2_v - 0.1) 

        # decrementing by 1 so to drop bad data
        logphidMin = logphidMin - 2
        if (logphid < logphidMin):
            # db = db.drop(i, axis=0)
            continue

        la_ind = int(round(2 * la_v)) + 6

        # db.iloc[ii, i_ver0] = ver0 
        # db.iloc[ii, i_ver1] = ver1 

        db.iloc[ii, i_la_ind] = la_ind
        db.iloc[ii, i_la_v] = la_v
         
        llc_ind = int(-2 * llc_v) - 1
        db.iloc[ii, i_llc_ind] = llc_ind

        row = db.loc[ii]
        if (resolutionFactor_ind == 0):
            dt_resaps = [0, 1]
            dt_ssoFS = [0, 1, 2, 3]
            dt_ssoFSV = [3, 6, 10, 0]
        else:
            dt_resaps = [dt_resap_ind]
            dt_ssoFS = [dt_ssoFS_ind]
            dt_ssoFSV = [dt_ssoFS_v]
        sz_dt_resaps = len(dt_resaps)
        sz_dt_ssoFS = len(dt_ssoFS)

        for resaps in dt_resaps:
            for j, ssoFS in enumerate(dt_ssoFS):
                rowNew = row.copy()
                ssoV = dt_ssoFSV[j]
                rowNew[i_ssoFS_ind] = ssoFS 
                rowNew[i_ssoFS_v] = ssoV 
                rowNew[i_dt_resap_ind] = resaps 
                rowNew[i_dt_resap_v] = resaps

                verN = nn_resolutionFactor * resolutionFactor_ind + \
                nn_llc * llc_ind + \
                nn_la * la_ind + \
                nn_ssoFS * ssoFS + \
                nn_resap + resaps
                
                rowNew[i_ver0] = verN
                # rowNew[i_ver1] = verN
#                pdOut = pdOut.append(rowNew, ignore_index=True)
                new_df = pd.DataFrame([rowNew])
                pdOut = pd.concat([pdOut, new_df], axis=0, join='outer', ignore_index=True)
#        ii += 1
    return pdOut

    #db.to_csv(filename, index=False,header=True)
    #db.to_csv("data/x.csv", index=False,header=True)

#def updata_csvs(root = "data/2023_03_22"):
#        sources = root + "/*.{}"
#        file_names = glob.glob(sources.format('csv') )
#        for fn in file_names:
#            updata_csv(fn)

# this includes all random field realization data for 1 combination of input scalar values (e.g. dd2, ...)
class OneScalarSetData:
    rawName = "raw"
    def __init__(self):
        # start and end (non-inclusive) range of rows of all.csv file that correspond to this combination of values
        self.st = -1
        self.en = -1
        # indices and values for scalar inputs corresponding to this data set
        self.scalar_input_indices = []
        self.scalar_input_vals = []


class Characteristics_data:

    lineColors = ["black", "blue", "red", "darkgreen", "peru", "darkorchid", "gray", "pink", "violet", \
                  "cyan", "firebrick", "royalblue", "lime", "darkslategrey", "orange"]
    markers = ["o", "s", "D", "p", "^", "<", ">", "P", "d", "X", "+", "x", "*", ".", "1", "2", "3", "4", "8"]

    def __init__(self):
        self.print_spatial_field_csv = True #True #change
        self.add_spatial_field_stat = True

        self.pd_data = pd.DataFrame()
        self.nm_ind_sortingFields = [] # names of indices of scalar inputs (ind)
        self.nm_inp_sortingFields = [] # names of values of scalar inputs (inp_s_)

        self.pos_ind_sortingFields = [] # positions of indices of scalar inputs (ind)
        self.pos_inp_sortingFields = [] # positions of values of scalar inputs (inp_s_)

        # for sorting and going from multiindex to single index and back
        self.num_ind_sorting = 0    # number of scalar values for sorting (e.g. 4 when llc, dd2, ldelc, and la are used)
        self.sz_ind_sortings = []   # size (e.g. maximum index + 1) for each index
        self.s2mi = SMIndex()
        # set of data for each random field set realization
        # self.runNoSets
        self.runNoSets_validPos = [] 

        # initial field related data
        self.pd_data_w_iniField = []
        self.folderSource = ''
        self.folderDest = ''

    def Main_Characteristics_data(self, folderSource = "data/2023_03_22_source", folderDest = "data/Characteristics_data"):
        self.folderSource = folderSource 
        self.folderDest = folderDest
        fn = self.folderDest + "/allData/mainMembers.txt"
        if os.path.isfile(fn):
            return    
        self.GenerateOneSortedFile()
        self.FormMultiIndexMatrix_AfterSortedMatrix()
        self.WriteTextData()
        self.Add_SpatialFieldData()
        self.CalculateStatesOverSameInputParaSet()

    def Main_Plot_Scatter(self, folderSource = "../data/2023_03_22_source", folderDest = "../data/Characteristics_data"):
        self.folderSource = folderSource
        self.folderDest = folderDest
        self.folderSource = folderSource 
        self.folderDest = folderDest
        fn = self.folderDest + "/allData/mainMembers.txt"
        if not os.path.isfile(fn):
            self.GenerateOneSortedFile()
            self.FormMultiIndexMatrix_AfterSortedMatrix()
            self.WriteTextData()
            self.print_spatial_field_csv = True # change
            self.add_spatial_field_stat = True
            self.Add_SpatialFieldData()
        else:
            fdestAll = self.folderDest + "/allData"
            self.GenerateOneSortedFile()
            self.FormMultiIndexMatrix_AfterSortedMatrix()
            allFile_wxs = fdestAll + "/all_wxs.csv"
            if (os.path.isfile(allFile_wxs)):
                self.pd_data_w_iniField = pd.read_csv(allFile_wxs)
            else:
                print(f"Cannot open file {allFile_wxs}")

        self.PlotScatter()

    def WriteTextData(self):
        print("start WriteTextData ...\n")
        fn = self.folderDest + "/allData/mainMembers.txt"
        with open(fn, "w") as fl:
            print("runNoSets_validPos\n",file=fl)
            np.savetxt(fl, self.runNoSets_validPos, fmt='%d', delimiter=",")
            print("s2mi\n",file=fl)
            self.s2mi.savetxt(fl)

    def GenerateOneSortedFile(self):
        print("start GenerateOneSortedFile ...\n")
        # sortingOrder = [True, True, True, True] # Don't change this as the rest of multi-index depends on all True sorting
        self.nm_ind_sortingFields = []
        self.nm_inp_sortingFields = []
        for v in sortingFields:
            self.nm_ind_sortingFields.append("ind" + v)
            self.nm_inp_sortingFields.append("inp_s_" + v)

        sf = self.nm_ind_sortingFields.copy()
        #so = sortingOrder
        sf.append("runNo")
        #so.append("True")
        fdestAll = self.folderDest + "/allData"
        allFile = fdestAll + "/all_raw.csv"
        if os.path.isfile(allFile):
            self.pd_data = pd.read_csv(allFile)
        else:
            rutil.rmkdir(self.folderDest)
            rutil.rmkdir(fdestAll)

            # generate One sorted file
            file_names = glob.glob(self.folderSource + "/*.{}".format('csv') )
            # print('2')
            # time.sleep(2)
            # print('file_names = ', file_names)
            # print('self.folderSource = ', self.folderSource)
            # print(self.folderSource)
            # time.sleep(10)

            cnt = 0
            for fn in file_names:
                print(fn)
                if cnt==0:
                    self.pd_data = read_updata_csv(fn)
                else:
                    db = read_updata_csv(fn)
                    self.pd_data = pd.concat([self.pd_data, db], ignore_index=True)
                print("cnt: %d" % cnt)
                cnt = cnt +1

            self.pd_data = self.pd_data.apply(pd.to_numeric, errors='coerce')

            trues = [True for _ in range(len(sf))]
            self.pd_data = self.pd_data.sort_values(by=sf, ascending=trues)
#            if ((mainOption == mo_frac_rate_cf) or (mainOption == mo_frac_rate_f)):
#                    self.pd_data = self.pd_data.sort_values(by=sf, ascending=[True, True, True, True, True])
#            elif (mainOption == mo_frac_res_x):            
#                    self.pd_data = self.pd_data.sort_values(by=sf, ascending=[True, True, True, True, True, True])

            for key, logbase in dist_logs.items():
                if (logbase <= 0):
                    continue
                clmn = "out_s_" + key
                try:
                    loc = self.pd_data.columns.get_loc(clmn)
                    clmn_new = "out_s_log_" + key
                    self.pd_data.insert(loc + 1, clmn_new, np.log(self.pd_data[clmn]) / np.log(logbase))
                except KeyError:
                    print(f"Column {clmn} does not exist in the DataFrame")                    

            self.pd_data.to_csv(fdestAll + "/all_raw.csv", index=False,header=True)

    def FormMultiIndexMatrix_AfterSortedMatrix(self):
        print("start FormMultiIndexMatrix_AfterSortedMatrix ...\n")
        self.pos_ind_sortingFields = []
        for v in self.nm_ind_sortingFields:
            self.pos_ind_sortingFields.append(self.pd_data.columns.get_loc(v))
        self.pos_inp_sortingFields = []
        for v in self.nm_inp_sortingFields:
            self.pos_inp_sortingFields.append(self.pd_data.columns.get_loc(v))
        self.num_ind_sorting = len(self.nm_inp_sortingFields)
        self.sz_ind_sortings = np.zeros(self.num_ind_sorting, dtype=int)
        for i, cn in enumerate(self.nm_ind_sortingFields):
            unique = self.pd_data[cn].unique()
            unique.sort()
            sz = len(unique)
            self.sz_ind_sortings[i] = (int)(unique[sz - 1] + 1)
        self.s2mi.Initialize_MultiIndex(self.sz_ind_sortings)
        totSize = self.s2mi.totalSz
        self.runNoSets = np.empty(totSize, dtype=object)
        self.runNoSets_validPos = [] 
        
        cur_linPos = -1
        numRows = self.pd_data.shape[0]
        numCols = self.pd_data.shape[1]
        for i in range(numRows):
            mind = []
            for jc in self.pos_ind_sortingFields:
                mind.append((int)(self.pd_data.iloc[i, jc]))
            linPos = int(self.s2mi.MI2SI(mind))
            if (linPos != cur_linPos):
                self.runNoSets[linPos] = OneScalarSetData()
                self.runNoSets[linPos].st = i
                self.runNoSets[linPos].scalar_input_indices = mind
                self.runNoSets[linPos].scalar_input_vals = []
                for jc in self.pos_inp_sortingFields:
                    self.runNoSets[linPos].scalar_input_vals.append(self.pd_data.iloc[i, jc])
                cur_linPos = linPos
                self.runNoSets_validPos.append(cur_linPos)
            self.runNoSets[linPos].en = i + 1

    def Add_SpatialFieldData(self):
        print("start Add_SpatialFieldData ...\n")
        fdestAll = self.folderDest + "/allData"
        allFile = fdestAll + "/all_raw.csv"
        allFile_wxs = fdestAll + "/all_wxs.csv"
        fieldFile = fdestAll + "/all_fx.csv"
        if ((self.print_spatial_field_csv == True) and os.path.isfile(fieldFile)):
            self.print_spatial_field_csv = False
        if ((self.add_spatial_field_stat) and os.path.isfile(allFile_wxs)):
            self.pd_data_w_iniField = pd.read_csv(allFile_wxs)
            self.pd_data_w_iniField = self.pd_data_w_iniField.apply(pd.to_numeric, errors='coerce')
            self.add_spatial_field_stat = False
        elif (self.add_spatial_field_stat == False):
            self.pd_data_w_iniField = self.pd_data
        if ((self.print_spatial_field_csv == False) and (self.add_spatial_field_stat == False)):
            return
        if (self.add_spatial_field_stat == True):
            self.pd_data_w_iniField = self.pd_data.copy()
        nRows = self.pd_data.shape[0]
        
        # now adding the spatial field and a separate file for input scalar fields

        nScalarFld, env = 0, 0
        if (self.print_spatial_field_csv == True):
            fieldFile = fdestAll + "/all_fx.csv"
            eol = "\n" #os.linesep
            f = open(fieldFile, 'w')
            nScalarFld = (int)(2**meshp2)
            if (valsAtVert):
                nScalarFld+= 1
            env = nScalarFld - 1
            for j in range(env):
                f.write(f"inp_fx_iniSc{(int)(j)},")
            f.write(f"inp_fx_iniSc{(int)(env)}{eol}")


        self.scalar_out_pos = [i for i, x in enumerate(self.pd_data.columns) if "out_s_" in x]
        self.scalar_out_field = [i for i, x in enumerate(self.pd_data.columns) if "out_f_" in x]

        if (self.add_spatial_field_stat == True):
            pre_nameMapped = "inp_sx_sim_"
            pre_nameRaw = "inp_sx_raw_"
            addMapped = True
            addRaw = True
            names = InpF.GetStatNames(pre_nameMapped, pre_nameRaw, addMapped, addRaw) 
            sz_stat_fields = len(names)
            first_out_pos = self.scalar_out_pos[0]
            for i in range(len(self.scalar_out_pos)):
                self.scalar_out_pos[i] += sz_stat_fields
            for i in range(len(self.scalar_out_field)):
                self.scalar_out_field[i] += sz_stat_fields
            first_out_pos_old = first_out_pos
            first_out_pos = self.scalar_out_pos[0]
            fx_col_pos = range(first_out_pos_old, first_out_pos)
            for i, fxsi in enumerate(fx_col_pos):
                self.pd_data_w_iniField.insert(fxsi, names[i], np.nan)
        
        # now looping over rows get the random field realization and process it:
        hasDeltFactor = (mainOption == mo_frac_res_x_w_delc_fact)
        iCol_delC_fact = -1
        iCol_shape = -1
        if (rateStudy):
            iCol_llc = self.pd_data.columns.get_loc("inp_s_llc") #self.pos_inp_sortingFields[0]
            iCol_dd2 = self.pd_data.columns.get_loc("inp_s_dd2") #self.pos_inp_sortingFields[1]
            iCol_la = self.pd_data.columns.get_loc("inp_s_la") #self.pos_inp_sortingFields[3]
            iCol_resF = -1
            if (mainOption == mo_frac_rate_f_wShape):
                iCol_shape = self.pd_data.columns.get_loc("inp_s_shape")
        elif ((mainOption == mo_frac_res_x) or hasDeltFactor):
            if (hasDeltFactor):
                iCol_delC_fact = self.pd_data.columns.get_loc("inp_s_delc_Tf_fact")
            iCol_llc = self.pd_data.columns.get_loc("inp_s_llc")
            iCol_dd2 = self.pd_data.columns.get_loc("inp_s_dd2")
            iCol_la = self.pd_data.columns.get_loc("inp_s_la")
            iCol_resF = self.pd_data.columns.get_loc("inp_s_resolutionFactor")
            iCol_sso = self.pd_data.columns.get_loc("inp_s_ssoFS")
            
        for i in range(nRows):
            if (i % 1000 == 0):
                print(f"i {i} / {nRows}\n")
            llc = self.pd_data.iloc[i, iCol_llc]
            dd2 = self.pd_data.iloc[i, iCol_dd2]
            la = self.pd_data.iloc[i, iCol_la]
            meshp2_4Simulation = -1
            # change:
            sso = reduction_sso
            shape = 2
            useOriginalSN_Mesh4Raw = 1
            OneScalarSetData.rawName = "origSN"
            if (mainOption == mo_frac_rate_cf):
                if (la < 0.0): 
                    if (la < -2.4):
                        # -3.0 -> 16
                        # -2.5 -> 16
                        meshp2_4Simulation = 10
                    elif (la < -1.9):
                        # -2.0 -> 8
                        meshp2_4Simulation = 11
                    elif (la < -1.4):
                        # -1.5 -> 4
                        meshp2_4Simulation = 12
                    else:
                        # -1.0 -> 2
                        # -0.5 -> 2
                        meshp2_4Simulation = 13
                useOriginalSN_Mesh4Raw = 0 # 1 uses original mesh for "raw" mesh (i.e. instead of using a possibly a lower resolution mesh); 2 does not calculate raw/original stat altogether 
                OneScalarSetData.rawName = "raw"
            elif ((mainOption == mo_frac_res_x) or hasDeltFactor):
                resF = self.pd_data.iloc[i, iCol_resF]
                meshp2_4Simulation = 14 - resF
                sso = self.pd_data.iloc[i, iCol_sso]
                useOriginalSN_Mesh4Raw = 0 # 1 uses original mesh for "raw" mesh (i.e. instead of using a possibly a lower resolution mesh); 2 does not calculate raw/original stat altogether 
                OneScalarSetData.rawName = "raw"
            if (mainOption == mo_frac_rate_f_wShape):
                shape = self.pd_data.iloc[i, iCol_shape]
    
            serNo = self.pd_data.iloc[i, 0]
            pinp = InpF()
            pinp.Initialize_InpF(valsAtVert, meshp2, serNo, llc, dd2, \
                                isPeriodic, sso, \
                                meshp2_4Simulation, meshp2_4Output, shape, useOriginalSN_Mesh4Raw)
            if (self.add_spatial_field_stat == True):
                stat_vals = pinp.GetStatVecVals(addRawStat)
                for ci, fxsi in enumerate(fx_col_pos):
                    self.pd_data_w_iniField.iloc[i, fxsi] = stat_vals[ci]

            if (self.print_spatial_field_csv == True):
                for j in range(env):
                    f.write(f"{pinp.vals4Output[j]},")
                f.write(f"{pinp.vals4Output[env]}{eol}")

        if (self.add_spatial_field_stat == True):
            self.pd_data_w_iniField.to_csv(allFile_wxs, index=False,header=True)
        else:
            shutil.copy(allFile, allFile_wxs)
        # Close the file
        if (self.print_spatial_field_csv == True):
            f.close()

    # calculating mean, max, ... of dataset
    def CalculateStatesOverSameInputParaSet(self):
        print("start CalculateStatesOverSameInputParaSet ...\n")
        self.pd_data_w_iniField
        # create multiple files for
        # 0 -> mean        # 1 -> min        # 2 -> max        # 3 -> cov        # 4 -> sdiv 
        # uses StatOfVec.GetVecVals()
        retVafter_std = True
        namesStat = []
        fdestAll = self.folderDest + "/allData"
        namesStat = StatOfVec.Get_Names("_", namesStat, retVafter_std)
        if (os.path.isfile(fdestAll + "/stat" + namesStat[0] + ".csv")):
            return
        szStat = len(namesStat)
        colnames = self.pd_data_w_iniField.columns
        self.pdStats = np.empty(szStat, dtype=object)
        numColsStat = len(colnames)
        numRowsStat = len(self.runNoSets_validPos)
        self.pdStat_vals = np.empty([numRowsStat, numColsStat], dtype=object)

        self.scalar_out_pos_wfx = [i for i, x in enumerate(self.pd_data_w_iniField.columns) if "out_s_" in x]
        stat_calStart = self.scalar_out_pos_wfx[0]
        sfx_cols = [i for i, x in enumerate(self.pd_data_w_iniField.columns) if "inp_sfx_" in x]
        if (len(sfx_cols) > 0):
            stat_calStart = sfx_cols[0]

        for si in range(szStat):
            pdd = pd.DataFrame(columns=colnames)
            for cnm in pdd.columns:
                pdd[cnm] = np.full(numRowsStat, np.nan)
            self.pdStats[si] = pdd
        
        # now generating data
        # dv_ri stands for "distinct valid row index" referring to the compbinations of scalar input that have valid input in the data set
        for cntr, dv_ri in enumerate(self.runNoSets_validPos):
            # start and end row number in allRaw file:
            st_row = self.runNoSets[dv_ri].st
            en_row = self.runNoSets[dv_ri].en
            sz_row = en_row - st_row
            # now taking care of other stats
            for si in range(szStat):
                for cj in range(stat_calStart):
                    self.pdStats[si].iloc[cntr, cj] = self.pd_data_w_iniField.iloc[st_row, cj]
                self.pdStats[si].iloc[cntr, 0] = sz_row # count how many values are in this set
            # now actually calculating the stats for other columns
            for cj in range(stat_calStart):
                vl = self.pd_data_w_iniField.iloc[st_row, cj]
                vc = np.full(sz_row, vl).tolist()
                self.pdStat_vals[cntr, cj] = vc

            for cj in range(stat_calStart, numColsStat):
                vecVals = np.full(sz_row, np.nan)
                for ri in range(st_row, en_row):
                    vecVals[ri - st_row] = self.pd_data_w_iniField.iloc[ri, cj]
                self.pdStat_vals[cntr, cj] = vecVals
                sov = StatOfVec()
                sov.Compute_Vec_Stat(vecVals, False)
                statVals = sov.GetVecVals([], retVafter_std)
                for si in range(szStat):
                    self.pdStats[si].iloc[cntr, cj] = statVals[si]
            
        # now printing the stat files
        for si in range(szStat):
            self.pdStats[si].apply(pd.to_numeric, errors='coerce')
            self.pdStats[si].to_csv(fdestAll + "/stat" + namesStat[si] + ".csv", index=False,header=True)
        with open(fdestAll + "/stat_vals.pkl", 'wb') as f:
            pickle.dump(self.pdStat_vals, f)            

    def PlotScatter(self):
        posMean = self.pd_data_w_iniField.columns.get_loc("inp_sx_sim_mean")
        posMin = self.pd_data_w_iniField.columns.get_loc("inp_sx_sim_mn")
        poscov = self.pd_data_w_iniField.columns.get_loc("inp_sx_sim_cov")
        posstd = self.pd_data_w_iniField.columns.get_loc("inp_sx_sim_std")
        posD = self.pd_data_w_iniField.columns.get_loc("inp_sx_sim_hfD")
        pos_xCorNeg = self.pd_data_w_iniField.columns.get_loc("inp_sx_sim_ac_xCorNeg")
        pos_N = self.pd_data_w_iniField.columns.get_loc("inp_sx_sim_seg_num")
        pos_seg_std = self.pd_data_w_iniField.columns.get_loc("inp_sx_sim_seg_std")
        clmns = [posMin, posMean, poscov, posstd, posD, pos_xCorNeg, pos_N, pos_seg_std]
        nms = ["min", "mean", "cov", "std", "D", "xCor0", "N", "seg_len_std"]
        # col_out_indices = [i for i, name in enumerate(self.pd_data_w_iniField.columns) if (('out_' in name) or ("inp_sfx_" in name))]
        col_out_indices = [i for i, name in enumerate(self.pd_data_w_iniField.columns) if ('out_' in name)]
        # sz = min(len(nms), 4)
        sz = len(nms)
        for ci in range(sz):
            self.PlotScatterAux(ci, clmns[ci], nms[ci], col_out_indices)


    def PlotScatterAux(self, ci, clmn, nm, col_out_indices):
        # sortingFields = ["llc", "dd2", "ldelc", "la"]
        plotRoot = "scatters_" + str(ci) + nm 
        rutil.rmkdir(plotRoot)
        plotRoot += "/"
        xlbl = nm + "(s(x))"

        totalIndSize = len(self.sz_ind_sortings)
        fldNames = []
        for cyi, cy in enumerate(col_out_indices):
            nmTmp = "fld_" + str(cyi) + "_" + self.pd_data_w_iniField.columns[cy]
            fldNames.append(nmTmp)  
            rutil.rmkdir(plotRoot + nmTmp)

        outer_indices = [2, 3, 1]
        inner_indices = [0]
        sz_outer = len(outer_indices)
        directionalSizes_outer = []
        for ind in outer_indices:
            directionalSizes_outer.append(self.sz_ind_sortings[ind])
        smi_outer = SMIndex()
        smi_outer.Initialize_MultiIndex(directionalSizes_outer)

        sz_inner = len(inner_indices)
        directionalSizes_inner = []
        for ind in inner_indices:
            directionalSizes_inner.append(self.sz_ind_sortings[ind])
        smi_inner = SMIndex()
        smi_inner.Initialize_MultiIndex(directionalSizes_inner)

        for olin in range(smi_outer.totalSz):
            omulti = smi_outer.SI2MI(olin)

            outName = ""
            for [ocntr, oind] in enumerate(outer_indices):
                str_ind = str(omulti[ocntr])
                strTmp = "I_" + sortingFields[oind] + "=" + str_ind
                if (ocntr > 0):
                    outName += "_"
                outName += strTmp

            numRuns = 0
            sts = []
            ens = []
            imultis = []
            ilins = []
            linIndexAlls = []
            for ilin in range(smi_inner.totalSz):
                imulti = smi_inner.SI2MI(ilin)
                # getting full index of this run
                indTotal = np.zeros(totalIndSize)
                for jj, po in enumerate(outer_indices):
                    indTotal[po] = omulti[jj]
                for jj, pi in enumerate(inner_indices):
                    indTotal[pi] = imulti[jj]

                # now go from total index to lin index over all runs
                linIndexAll = self.s2mi.MI2SI(indTotal)
                # see if there is data for this
                if (self.runNoSets[linIndexAll] == None):
                    continue
                # so there is data in it
                numRuns += 1
                st = self.runNoSets[linIndexAll].st
                en = self.runNoSets[linIndexAll].en
                sts.append(st)
                ens.append(en)
                imultis.append(imulti)
                ilins.append(ilin)
                linIndexAlls.append(linIndexAll)

            if (numRuns == 0):
                continue

            rutil.rmkdir(plotRoot + outName)

            for cyi, cy in enumerate(col_out_indices):
                fig, ax = plt.subplots()
                ax.set_adjustable('box')
                of = plotRoot + outName 
                rutil.rmkdir(of)
                fldName = fldNames[cyi]
                nm_out_fld = of + "/" + outName + fldName
                nm_fld_out = plotRoot + "/" + fldName + "/" + fldName + outName

                lc = Characteristics_data.lineColors[0]
                mkr = Characteristics_data.markers[0]

                for rn in range(numRuns):
                    st = sts[rn]
                    en = ens[rn]
                    imulti = imultis[rn]
                    ilin = ilins[rn]
                    linIndexAll = linIndexAlls[rn]
                    xVals = self.pd_data_w_iniField.iloc[st:en, clmn]
                    if (sz_inner > 0):
                        lc = Characteristics_data.lineColors[imulti[0]]
                        if (sz_inner > 1):
                            mkr = Characteristics_data.markers[imulti[1]]
                    lbl = ""
                    if (sz_inner == 1):
                        lbl = str(self.runNoSets[linIndexAll].scalar_input_vals[0])
                    elif (sz_inner > 1):
                        for [icntr, iind] in enumerate(inner_indices):
                            str_vl = str(self.runNoSets[linIndexAll].scalar_input_vals[iind])
                            strTmp = sortingFields[iind] + "=" + str_vl
                            if (icntr > 0):
                                lbl += ","
                            lbl += strTmp

                    yVals = self.pd_data_w_iniField.iloc[st:en, cy]
                    plt.scatter(xVals, yVals, marker=mkr, color=lc, label=lbl)

                if sz_inner > 0:
                    plt.legend()
                plt.xlabel(xlbl)
                plt.ylabel(fldName)
                # fig.show()

                plotFullNameWOExt = nm_out_fld #
                pltName = plotFullNameWOExt + ".svg"
                plt.savefig(pltName, format='svg')
                pltName = plotFullNameWOExt + ".png"
                plt.savefig(pltName, format='png')
                fig.clf()
                plt.close()

                shutil.copy(plotFullNameWOExt + ".png", nm_fld_out + ".png")
                shutil.copy(plotFullNameWOExt + ".svg", nm_fld_out + ".svg")
