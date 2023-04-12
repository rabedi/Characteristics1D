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

# Change
sortingFields = ["llc", "dd2", "ldelc", "la"]

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
dist_logs["phi_d_tFin"] = 10
dist_logs["phi_d_tFin"] = 10
dist_logs["phi_d_norm_phi0_tFin"] = 10
dist_logs["phi_d_norm_phi_input_F"] = 10
dist_logs["psi_f"] = 10
dist_logs["psi_f_norm_psiC"] = 10
dist_logs["eps_f"] = 10
dist_logs["eps_f_norm"] = 10
dist_logs["lbar_F"] = 10

def read_updata_csv(filename = "../StochasticPostprocessor/data/2023_03_20/pps3Out_ub1.csv"):
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
        ver1 = ver0

        db.iloc[ii, i_ver0] = ver0 
        db.iloc[ii, i_ver1] = ver1 

        db.iloc[ii, i_la_ind] = la_ind
        db.iloc[ii, i_la_v] = la_v
         
        db.iloc[ii, i_dd2_ind] = dd2_ind
        db.iloc[ii, i_dd2_v] = dd2_v

        llc_ind = int(-2 * llc_v) - 1
        db.iloc[ii, i_llc_ind] = llc_ind
        ii += 1
    nRows = db.shape[0]
    return db
    #db.to_csv(filename, index=False,header=True)
    #db.to_csv("data/x.csv", index=False,header=True)

#def updata_csvs(root = "data/2023_03_22"):
#        sources = root + "/*.{}"
#        file_names = glob.glob(sources.format('csv') )
#        for fn in file_names:
#            updata_csv(fn)

# this includes all random field realization data for 1 combination of input scalar values (e.g. dd2, ...)
class OneScalarSetData:
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
        self.print_spatial_field_csv = False #True
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
        self.GenerateOneSortedFile()
        self.FormMultiIndexMatrix_AfterSortedMatrix()
        self.WriteTextData()
        self.print_spatial_field_csv = False
        self.add_spatial_field_stat = True
        self.Add_SpatialFieldData()
        self.PlotScatter()

    def WriteTextData(self):
        fn = self.folderDest + "/allData/mainMembers.txt"
        with open(fn, "w") as fl:
            print("runNoSets_validPos\n",file=fl)
            np.savetxt(fl, self.runNoSets_validPos, fmt='%d', delimiter=",")
            print("s2mi\n",file=fl)
            self.s2mi.savetxt(fl)

    def GenerateOneSortedFile(self):
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
            self.pd_data = self.pd_data.sort_values(by=sf, ascending=[True, True, True, True, True])

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
            pre_name = "inp_sfx_"
            addRawStat = True
            names = InpF.GetStatNames(pre_name, addRawStat)
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
        iCol_llc = self.pos_inp_sortingFields[0]
        iCol_dd2 = self.pos_inp_sortingFields[1]
        iCol_la  = self.pos_inp_sortingFields[3]
        for i in range(nRows):
            llc = self.pd_data.iloc[i, iCol_llc]
            dd2 = self.pd_data.iloc[i, iCol_dd2]
            la = self.pd_data.iloc[i, iCol_la]
            meshp2_4Simulation = -1
            # change:
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

            serNo = self.pd_data.iloc[i, 0]
            pinp = InpF()
            pinp.Initialize_InpF(valsAtVert, meshp2, serNo, llc, dd2, \
                                isPeriodic, reduction_sso, \
                                meshp2_4Simulation, meshp2_4Output)
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
                sov.Compute_Vec_Stat(vecVals)
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
        posMean = self.pd_data_w_iniField.columns.get_loc("inp_sfx_sim_mean")
        posMin = self.pd_data_w_iniField.columns.get_loc("inp_sfx_sim_mn")
        poscov = self.pd_data_w_iniField.columns.get_loc("inp_sfx_sim_cov")
        posstd = self.pd_data_w_iniField.columns.get_loc("inp_sfx_sim_std")
        clmns = [posMin, posMean, poscov, posstd]
        nms = ["min", "mean", "cov", "std"]
        # col_out_indices = [i for i, name in enumerate(self.pd_data_w_iniField.columns) if (('out_' in name) or ("inp_sfx_" in name))]
        col_out_indices = [i for i, name in enumerate(self.pd_data_w_iniField.columns) if ('out_' in name)]
        sz = min(len(nms), 4)
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
