import numpy as np
import pandas as pd
import pickle
import glob
import math
import os
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import shutil
import time

import rmulti_index as rmi 
from rmulti_index import SMIndex as SMIndex
import input_field
from input_field import InpF as InpF
#from enum import Enum
import rutil as rutil
from characteristics_data import Characteristics_data as Characteristics_data 
import characteristics_data as charDat 

from input_field import StatOfVec as StatOfVec
from input_field import InpFsOuput as InpFsOuput

from fragmentation import Frag1D1la
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

from rstats import  Test_HD_Cor

from fragmentation import Frag1D1la
from fragmentation import Frag1Dlap_curves

class plParameters:
    lineStyles = ["solid", "dashed", "dashdot", "dotted", (0, (1, 10)), (0, (1, 1)), \
                 (5, (10, 3)), (0, (5, 10)), (0, (5, 1)), (0, (3, 10, 1, 10)), \
                 (0, (3, 1, 1, 1)), (0, (3, 5, 1, 5, 1, 5)), (0, (3, 10, 1, 10, 1, 10)), \
                    (0, (3, 1, 1, 1, 1, 1))]
    #lineStyles = ["solid", "dashed", "dashdot", "dotted", "loosely dotted", "densely dotted", \
    #             "long dash with offset", "loosely dashed", "densely dashed", "loosely dashdotted", \
    #             "densely dashdotted", "dashdotdotted", "loosely dashdotdotted", "densely dashdotdotted", "dashdotted"]
    lineColors = ["black", "blue", "red", "darkgreen", "peru", "darkorchid", "gray", "pink", "violet", \
                  "cyan", "firebrick", "royalblue", "lime", "darkslategrey", "orange", "purple", "brown", \
                  "yellow", "green", "magenta",  \
                  "aliceblue","antiquewhite", "aqua", "aquamarine", "azure", \
                  "beige", "bisque", "blanchedalmond", \
                  "blueviolet", "burlywood","cadetblue", \
                  "chartreuse", "chocolate", "coral", "cornflowerblue", \
                  "cornsilk", "crimson", "darkblue", "darkcyan", \
                  "darkgoldenrod", "darkgray", \
                  "darkkhaki", "darkmagenta", "darkolivegreen", "darkorange", \
                  "darkorchid", "darkred", "darksalmon", "darkseagreen", \
                  "darkslateblue", "darkslategray", "darkslategrey", \
                  "darkturquoise", "darkviolet", "deeppink", "deepskyblue", \
                  "dimgray", "dimgrey", "dodgerblue", "firebrick", \
                  "floralwhite", "forestgreen", "fuchsia", "gainsboro", \
                  "ghostwhite", "gold", "goldenrod", \
                  "greenyellow", "honeydew", "hotpink", "indianred", "indigo", \
                  "ivory", "khaki", "lavender", "lavenderblush", "lawngreen", \
                  "lemonchiffon", "lightblue", "lightcoral", "lightcyan", \
                  "lightgoldenrodyellow", "lightgray", "lightgrey", \
                  "lightgreen", "lightpink", "lightsalmon", "lightseagreen", \
                  "lightskyblue", "lightslategray", "lightslategrey", \
                  "lightsteelblue", "lightyellow", "lime", "limegreen", \
                  "linen", "magenta", "maroon", "mediumaquamarine", \
                  "mediumblue", "mediumorchid", "mediumpurple", \
                  "mediumseagreen", "mediumslateblue", "mediumspringgreen", \
                  "mediumturquoise", "mediumvioletred", "midnightblue", \
                  "mintcream", "mistyrose", "moccasin", "navajowhite", "navy", \
                  "oldlace", "olive", "olivedrab", "orangered", \
                  "orchid", "palegoldenrod", "palegreen", "paleturquoise", \
                  "palevioletred", "papayawhip", "peachpuff", \
                  "plum", "powderblue", "rosybrown", \
                  "royalblue", "saddlebrown", "salmon", "sandybrown", \
                  "seagreen", "seashell", "sienna", "silver", "skyblue", \
                  "slateblue", "slategray", "slategrey", "snow", "springgreen", \
                  "steelblue", "tan", "teal", "thistle", "tomato", "turquoise", \
                  "wheat", "white", "whitesmoke", "yellowgreen"]
    markers = ["o", "s", "D", "p", "^", "<", ">", "P", "d", "X", "+", "x", "*", ".", "1", "2", "3", "4", "8"]
    filltyles = ["full", "top", "bottom", "left", "right", "full", "top", "bottom", "left", "right"] #"none"
    linewidths = [1, 2, 3, 4, 5, 1, 2, 3, 1, 2, 3]
    markersizes = [3, 5, 6, 7, 8, 9, 10, 11, 12]
    #markerfacecolors = ["none", "auto", ]

    @staticmethod
    # lineIn and lineOut are mlines.Line2D
    def UpdateLineProp(lineIn, propName, propIndex):
#        lineOut = mlines.Line2D([], [])
        lineOut = lineIn
        if (propName == "ls"):
            lineOut.set_linestyle(plParameters.lineStyles[propIndex])
        elif (propName == "c"):
            lineOut.set_color(plParameters.lineColors[propIndex])
        elif (propName == "marker"):
            lineOut.set_marker(plParameters.markers[propIndex])
        elif (propName == "fillstyle"):
            lineOut.set_fillstyle(plParameters.filltyles[propIndex])
        elif (propName == "lw"):
            lineOut.set_linewidth(plParameters.linewidths[propIndex])
        elif (propName == "ms"):
            lineOut.set_markersize(plParameters.markersizes[propIndex])
        elif (propName == "mfc"):
            if (propIndex == 0):
                lineOut.set_fillstyle("none")
                lineOut.set_markerfacecolor("none")
            elif (propIndex == 1):
                lineOut.set_markerfacecolor(lineOut.get_color())
            else:
                lineOut.set_markerfacecolor(plParameters.lineColors[(propIndex + 10) % 15])
        return lineOut
    @staticmethod
    def GetLineProperty_from_dict(dict):
        lineOut = mlines.Line2D([], [])
        try:
            value = dict["c"]
            lineOut = plParameters.UpdateLineProp(lineOut, "c", value)
        except KeyError:
            lineOut.set_color("black")

        for key, value in dict.items():
            if (key == "c"):
                continue
            lineOut = plParameters.UpdateLineProp(lineOut, key, value)
        return lineOut

# stores the range of indices and values of a particular column in a database generated by the characteristics code (or others following the same naming convention)
class DB_InputColumn:
    # nameBase is the base name of a data base (inp)
    # nameBase is like "la"
    def Initialize_DB_InputColumn(self, nameBase = "la", lineprop = "none"):
        self.nameBase = nameBase
        self.lineprop = lineprop
        if (self.nameBase != "runNo"):
            self.nameIndex = "ind" + self.nameBase
            self.nameValue = "inp_s_" + self.nameBase
        else:
            self.nameIndex = self.nameBase
            self.nameValue = self.nameBase

        #self.nameValue = "inp_" + self.nameBase

        # computed from input data 
        # 
        self.dbColumnIndex = 0
        self.dbColumnValue = 0
        # it is assumed indices go from zero to index Sz - 1
        self.indexSz = 0
        self.vals = []

# this is a scheme that in conjunction with PlotOMISplit can generate groups of 2D plots
# outputs -> refers to quantities that are changed outside the plot
# inner -> quantity that is varied across versionNos (combinations of input scalars) that creates a curve
#       >= 0 -> column for inner input column
#       <  0 -> data does not have this inner loop, e.g. plot of PDF for each run
# middle -> anything that is left. If emp
class PlotOMISplit_Instruction:
    def __init__(self):
        # the base names for split are given -> from which the indices below are computed
        # inputs:
        self.out_col_nameBases = []
        self.mid_col_nameBases = []
        self.in_col_nameBase = "none"
        splitNo = 0
        self.name = "splt_" + str(splitNo)
        # field name becomes before outer plot name
        # True (e.g. phid_llc=-2) where phid is field name and llc=-2 is the outer name
        # False will create llc=-2_phid
        self.fbeforeo = True
        # follow up to above, if mkfolder_fo == True, the name will be phid/llc=-2 and False, will be phid_llc=-2
        self.mkfolder_fo = True

        # computed
        # these indices are relative to scalar input parameters provided (e.g. from 0 to 6)
        self.out_col_inds = []
        self.mid_col_inds = []
        self.in_col_ind = -1 # -1 means it each version (combination of integrer inputs) give a curve (e.g. PDF) rather than a point

        # these indices are the absolute column positions of the scalar indices (not values) in the db files
        self.out_db_col_inds = []
        self.mid_db_col_inds = []
        self.in_db_col_ind = -1
        self.in_db_col_val = -1

        self.out_sz = 0
        self.mid_sz = 0
        self.has_in = False

# this is one curve info
class PlotOMISplit_Curve:
    def __init__(self):
        # this includes the legend entry and line properties BUT not the actual x and y
        self.oneCurve_indices = []
        self.oneCurve_vals = []
        self.oneCurve_dict_line_prop_2_ind = {}
        self.llabel = "" 
        self.row_numbers4curve = []

# this has a list of curves (PlotOMISplit_Curve) in it 
class PlotOMISplit_1Plot:
    def __init__(self):
        self.onePlot_indices = []
        self.onePlot_vals = []
        self.onePlot_dict_line_prop_2_ind = {}
        # names that can be used in nameing files, first one uses index values
        self.plot_nameBaseFromIndices = ""
        # second one uses values
        self.plot_nameBaseFromVals = ""
        self.plot_curves = [] # vector of curves
        self.plot_title_base = "" 


# this is the actual plot instructions that uses PlotOMISplit_Instruction and one DataBaseSplits to form plot informations
class PlotOMISplit_AllPlots:
    def __init__(self):
        self.colIndices = []
        self.colNames = []
        self.colLinepropKey = []

        # this is an mx2 matrix where m is the number of independent scalar value combinations
        #   1st column is plotNumber for this row
        #   2nd column is curveNumber
        # the opposite direction (i.e. from plot/curve to what rows it includes) is stored in PlotOMISplit_Curve.row_numbers4curve
        self.db_row2_plot_curveNos = []


        self.all_plots = [] # a vector of 1Plots
        self.smi_out = SMIndex()
        self.smi_mid = SMIndex()

class DataBaseSplits:
    # inpColDict is a dictionary of input columns
    def Initialize_DataBaseSplits(self, db, inpColDict, dist_logs = {}, folderDest = "data/Characteristics_data", plotFillMode = 0):
        self.plotFillMode = plotFillMode
        self.col_out_indices = [i for i, name in enumerate(db.columns) if (('out_' in name) or ("inp_sfx_" in name))]
        self.num_out = len(self.col_out_indices)

        self.db = db
        if (self.plotFillMode == 0):
            summaryName_min = folderDest  + "/allData/stat_mn.csv"
            summaryName_max = folderDest  + "/allData/stat_mx.csv"
            try:
                self.db_min = pd.read_csv(summaryName_min)
                try:
                    self.db_max = pd.read_csv(summaryName_max)
                    self.db_has_db_mM = 1
                except IOError:
                    self.db_has_db_mM = 0
            except IOError:
                self.db_has_db_mM = 0
        else:
            summary_std = folderDest  + "/allData/stat_std.csv"
            try:
                self.db_std = pd.read_csv(summary_std)
                self.db_max = self.db + self.db_std
                self.db_min = self.db - self.db_std
            except IOError:
                self.db_has_db_mM = 0

        self.db_inp_cols = []
        self.db_colnames = []
        self.db_logs = []
        repParts = ["out_s_", "out_", "inp_s_", "inp_"]
        for col in db.columns:
            nname = rutil.changeName(col, repParts)
            self.db_colnames.append(nname)
            try:
                lv = dist_logs[nname]
            except KeyError:
                lv = -1
            self.db_logs.append(lv)            

        for key, value in inpColDict.items():
            inp_col = DB_InputColumn()
            inp_col.Initialize_DB_InputColumn(key, value)
            try:
                inp_col.dbColumnIndex = db.columns.get_loc(inp_col.nameIndex)
                inp_col.dbColumnValue = db.columns.get_loc(inp_col.nameValue)
                unique = db[inp_col.nameIndex].unique()
                unique.sort()
                sz = len(unique)
                inp_col.indexSz = (int)(unique[sz - 1] + 1)
                unique_names = range(inp_col.indexSz)
                unique_row_indices = []
                for name in unique_names:
                    lst = db.index[db[inp_col.nameIndex] == name].tolist()
                    if (len(lst) > 0):
                        unique_row_indices.append(lst[0])
                    else:
                        unique_row_indices.append(-1)
                inp_col.vals = np.zeros(inp_col.indexSz, dtype=float)
                for j in range(inp_col.indexSz):
                    if (unique_row_indices[j] >= 0):
                        inp_col.vals[j] = db.iloc[unique_row_indices[j]][inp_col.dbColumnValue]
                    else:
                        inp_col.vals[j] = np.nan
#                inp_col.indexSz = len(unique)
#                unique_names = range(inp_col.indexSz)
#                unique_row_indices = [db.index[db[inp_col.nameIndex] == name].tolist()[0] for name in unique_names]
#                inp_col.vals = np.zeros(inp_col.indexSz, dtype=float)
#                for j in range(inp_col.indexSz):
#                    inp_col.vals[j] = db.iloc[unique_row_indices[j]][inp_col.dbColumnValue]
            except KeyError:
                print(f"Column {inp_col.nameIndex} cannot be found")
            self.db_inp_cols.append(inp_col)

        # splitInstructions is PlotOMISplit_Instruction
        # it returns updated splitInstructions
    def  Generate_AllPlots_fromSplitInstructions(self, splitInstructions):
        # updating the instructions based on db input columns
        totNumInputCols = len(self.db_inp_cols)
        colFound = np.zeros(totNumInputCols, dtype=bool)
        colFound.fill(False)

        splitInstructions.out_sz = len(splitInstructions.out_col_nameBases)
        splitInstructions.out_col_inds = np.zeros(splitInstructions.out_sz, dtype=int)
        splitInstructions.out_db_col_inds = np.zeros(splitInstructions.out_sz, dtype=int)
        for i in range(splitInstructions.out_sz):
            colBaseName = splitInstructions.out_col_nameBases[i]
            splitInstructions.out_col_inds[i] = -1
            for j in range(totNumInputCols):
                if (self.db_inp_cols[j].nameBase == colBaseName):
                    colFound[j] = True
                    splitInstructions.out_col_inds[i] = j
                    splitInstructions.out_db_col_inds[i] = self.db_inp_cols[j].dbColumnIndex
                    break

        sz_middle = totNumInputCols - splitInstructions.out_sz
        # inside
        splitInstructions.has_in = False
        splitInstructions.in_col_ind = -1
        splitInstructions.in_db_col_ind = -1
        splitInstructions.in_db_col_val = -1
        if (splitInstructions.in_col_nameBase != "none"):
            for j in range(totNumInputCols):
                if (self.db_inp_cols[j].nameBase == splitInstructions.in_col_nameBase):
                    colFound[j] = True
                    splitInstructions.in_col_ind = j
                    splitInstructions.in_db_col_ind = self.db_inp_cols[j].dbColumnIndex
                    if (self.db_inp_cols[j].nameBase != "runNo"):
                        splitInstructions.in_db_col_val = splitInstructions.in_db_col_ind + 1
                    else:
                        splitInstructions.in_db_col_val = splitInstructions.in_db_col_ind
                    splitInstructions.has_in = True
                    sz_middle -= 1
                    break

        # middle
        splitInstructions.mid_sz = len(splitInstructions.mid_col_nameBases)
        if (sz_middle == splitInstructions.mid_sz): # means middle base names are provided
            splitInstructions.mid_col_inds = np.zeros(splitInstructions.mid_sz, dtype=int)
            for i in range(splitInstructions.mid_sz):
                colBaseName = splitInstructions.mid_col_nameBases[i]
                splitInstructions.mid_col_inds[i] = -1
                for j in range(totNumInputCols):
                    if (self.db_inp_cols[j].nameBase == colBaseName):
                        colFound[j] = True
                        splitInstructions.mid_col_inds[i] = j
                        break
        else:
            splitInstructions.mid_col_nameBases = []
            splitInstructions.mid_col_inds = np.where(~colFound)[0]
            for i in splitInstructions.mid_col_inds:
                splitInstructions.mid_col_nameBases.append(self.db_inp_cols[i].nameBase)
            splitInstructions.mid_sz = len(splitInstructions.mid_col_inds)

        splitInstructions.mid_db_col_inds = np.zeros(splitInstructions.mid_sz, dtype=int)
        for jj in range(splitInstructions.mid_sz):
            pos = splitInstructions.mid_col_inds[jj] 
            splitInstructions.mid_db_col_inds[jj] = self.db_inp_cols[pos].dbColumnIndex
 
        allPlots = PlotOMISplit_AllPlots()
        szs = []
        for i in splitInstructions.out_col_inds:
            szs.append(self.db_inp_cols[i].indexSz)
        allPlots.smi_out.Initialize_MultiIndex(szs)

        szs = []
        for i in splitInstructions.mid_col_inds:
            szs.append(self.db_inp_cols[i].indexSz)
        allPlots.smi_mid.Initialize_MultiIndex(szs)

        sz_out = allPlots.smi_out.totalSz
        sz_mid = allPlots.smi_mid.totalSz

        allPlots.all_plots = np.empty(sz_out, dtype=object)
        allPlots.colIndices = splitInstructions.out_col_inds
        allPlots.colNames = [self.db_inp_cols[i].nameBase for i in allPlots.colIndices]
        allPlots.colLinepropKey = [self.db_inp_cols[i].lineprop for i in allPlots.colIndices]

        for i in range(sz_out):
            onePlot = PlotOMISplit_1Plot()
            onePlot.onePlot_indices = allPlots.smi_out.si2mi[i]
            onePlot.onePlot_dict_line_prop_2_ind = {}
            onePlot.onePlot_vals = []
            onePlot.plot_nameBaseFromIndices = ""
            onePlot.plot_nameBaseFromVals = ""

            for ii in range(splitInstructions.out_sz):
                ind = splitInstructions.out_col_inds[ii]
                indI = onePlot.onePlot_indices[ii]
                val = self.db_inp_cols[ind].vals[indI]
                tmpIndName = "I" + allPlots.colNames[ii] + "=" + str(indI)
                onePlot.plot_nameBaseFromIndices += tmpIndName
                if val.is_integer():
                    str_num = f"{val:.0f}"  # No decimal places for integers
                else:
                    str_num = f"{val:.1f}"                 
                tmpValName = "V" + allPlots.colNames[ii] + str_num
                onePlot.plot_nameBaseFromVals += tmpValName
                onePlot.onePlot_vals.append(val)
                onePlot.onePlot_dict_line_prop_2_ind[self.db_inp_cols[ind].lineprop] = onePlot.onePlot_indices[ii]

            onePlot.plot_curves = np.empty(sz_mid, dtype=object)
            hasMid = (sz_mid > 0)
            for j in range(sz_mid):
                oneCurve = PlotOMISplit_Curve()

                oneCurve.oneCurve_indices = allPlots.smi_mid.si2mi[j]
                oneCurve.oneCurve_dict_line_prop_2_ind = onePlot.onePlot_dict_line_prop_2_ind.copy()
                oneCurve.oneCurve_vals = []
                legentry = ""
                for jj in range(splitInstructions.mid_sz):
                    ind = splitInstructions.mid_col_inds[jj]
                    indj = oneCurve.oneCurve_indices[jj]
                    val = self.db_inp_cols[ind].vals[indj]
                    tmplegentry = self.db_inp_cols[ind].nameBase + "=" + str(val)
                    if (jj > 0):
                        legentry += ","
                    legentry += tmplegentry 
                    oneCurve.oneCurve_vals.append(val)
                    oneCurve.oneCurve_dict_line_prop_2_ind[self.db_inp_cols[ind].lineprop] = \
                        oneCurve.oneCurve_indices[jj]
                if (hasMid):
                    oneCurve.llabel = legentry

                onePlot.plot_curves[j] = oneCurve
            allPlots.all_plots[i] = onePlot

        # going over the rows of the database and storing the row in the corresponding plot number / curve number and vice versa
        num_rows = self.db.shape[0]
        allPlots.db_row2_plot_curveNos = np.zeros((num_rows, 2), dtype=int)
        for i in range(num_rows):
            # calculating plot position
            ind_out = np.zeros(splitInstructions.out_sz, dtype=int)
            for j in range(splitInstructions.out_sz):
                colPos = splitInstructions.out_db_col_inds[j]
                index = self.db.iloc[i, colPos]
                ind_out[j] = index
            plotPos = allPlots.smi_out.MI2SI(ind_out)
            
            # calculating curve position
            ind_mid = np.zeros(splitInstructions.mid_sz, dtype=int)
            for j in range(splitInstructions.mid_sz):
                colPos = splitInstructions.mid_db_col_inds[j]
                index = self.db.iloc[i, colPos]
                ind_mid[j] = index
            curvePos = allPlots.smi_mid.MI2SI(ind_mid)

            # db row number -> plot / curve
            allPlots.db_row2_plot_curveNos[i, 0] = plotPos
            allPlots.db_row2_plot_curveNos[i, 1] = curvePos

            # plot / curve -> db row number
            allPlots.all_plots[plotPos].plot_curves[curvePos].row_numbers4curve.append(i)

        return allPlots, splitInstructions

    # this is for plots that each row of the pandas DataFrame will provide only one value 
    # to generate a curve multiple lines are combined through looping over in attribute of splitInstructions
    # if fi = -1, x axis is the inner scalar parameters (e.g. la log of loading rate)
    # if >= 0 it refers to a column index of the data  set
    def  Plot_Scalar_Per_Row_Fi_vs_Fj(self, splitInstructions = PlotOMISplit_Instruction(), allPlots = PlotOMISplit_AllPlots(), fjs = 0, fis = -1, plotFill = 0, addExact = False, propName4MultiField = "c"):
        if (type(fjs) == int):
            fjs = [fjs]
        if (type(fis) == int):
            fis = [fis]
        if (fjs[0] == -1):
            for fj in self.col_out_indices:
                #if (fj != 140): log normalized fragment size
                #    continue # 
                self.Plot_Scalar_Per_Row_Fi_vs_Fj(splitInstructions, allPlots, [fj], fis, plotFill, addExact, propName4MultiField)
 
        # number of fill plot related
        num_fpr = 1
        if (plotFill):
            num_fpr = 3

        schemeFolder = splitInstructions.name
        rutil.rmkdir(schemeFolder)        
        fi = fis[0]
        fj = fjs[0]
        fjName = "fld" + str(fj) + self.db_colnames[fj]
        fiName = splitInstructions.in_col_nameBase
        isLog_phid = False
        isLog_spAve = False
        if (addExact): # (fiName == 'lap'):
            isLog_phid = (self.db_colnames[fj] == 'log_phi_d_tFin')
            isLog_spAve = (self.db_colnames[fj] == 'log_lbar_F')
        use_exactSln = (isLog_phid or isLog_spAve)

        sz_fj = len(fjs)
        moreThan1F = (sz_fj > 1)
        if (moreThan1F and (len(fis) == 1)):
            fis = np.zeros(sz_fj, dtype=int)
            fis.fill(fi)
        fiField = (fi >= 0)
        fifjName = fjName
        if (fiField):
            fiName = "fld" + str(fi) + self.db_colnames[fi]
            fifjName = fiName + "_vs_" + fjName

        for api, ap in enumerate(allPlots.all_plots):
            plName = ap.plot_nameBaseFromIndices
            tmpn = "plt" + str(api)
            plotFullNameWOExt = ""
            fnpart = ""
            if (splitInstructions.fbeforeo):
                if (splitInstructions.mkfolder_fo):
                    plotFolder = schemeFolder + "/" + tmpn + "_" + plName
                    rutil.rmkdir(plotFolder)        
                else:
                    plotFolder = schemeFolder
                fnpart =  tmpn + "_" + fifjName + "_" + plName
            else:
                if (splitInstructions.mkfolder_fo):
                    plotFolder = schemeFolder +  "/_field_" + fifjName
                    rutil.rmkdir(plotFolder)        
                else:
                    plotFolder = schemeFolder
                fnpart =  tmpn + "_" + plName + "_" + fifjName
            plotFullNameWOExt = plotFolder + "/" + fnpart
            
            # making the plot
            fig, ax = plt.subplots()
            ax.set_adjustable('box')

            xm = 1e10
            xM = -1e10
            ym = 1e10
            yM = -1e10

            if (use_exactSln):
                f1dc = Frag1Dlap_curves()
                hasLines, lns, txmn, txmx, tymn, tymx = f1dc.ReturnLines(isLog_phid, isLog_spAve)
                if (hasLines):
                    xm = txmn
                    xM = txmx
                    ym = tymn
                    yM = tymx
                    for ln in lns:
                        ax.add_line(ln)


            # ax.set_box_aspect(None)
            #ax.aspect('auto')
            pad = 0.05
            cntr_curves = 0
            yaxis_log = False
            xaxis_log = False
            if (splitInstructions.has_in): # inside gives point data
                for find in range(sz_fj):
                    fi = fis[find]
                    fj = fjs[find]
                    for curvei, curve in enumerate(ap.plot_curves):
                        sz_ptrs = len(curve.row_numbers4curve)
                        if (sz_ptrs == 0):
                            continue
                        if (plotFill):
                            xms = np.zeros(sz_ptrs, dtype=float)
                            yms = np.zeros(sz_ptrs, dtype=float)
                            xMs = np.zeros(sz_ptrs, dtype=float)
                            yMs = np.zeros(sz_ptrs, dtype=float)
                        for fpri in range(num_fpr):                        
                            x = np.zeros(sz_ptrs, dtype=float)
                            y = np.zeros(sz_ptrs, dtype=float)
                            colx = fi
                            if (not fiField):
                                colx = splitInstructions.in_db_col_val
                            for indri, ri in enumerate(curve.row_numbers4curve):
                                x[indri] = self.db.iloc[ri, colx]
                            if (fpri == 0): # average - or no fill
                                for indri, ri in enumerate(curve.row_numbers4curve):
                                    y[indri] = self.db.iloc[ri, fj] # + np.random.rand()
                            elif (fpri == 1):  # minimum value
                                for indri, ri in enumerate(curve.row_numbers4curve):
                                    y[indri] = self.db_min.iloc[ri, fj] # + np.random.rand()
                                if (fiField):
                                    for indri, ri in enumerate(curve.row_numbers4curve):
                                        x[indri] = self.db_min.iloc[ri, colx]
                            elif (fpri == 2):  # maximum value
                                for indri, ri in enumerate(curve.row_numbers4curve):
                                    y[indri] = self.db_max.iloc[ri, fj] # + np.random.rand()
                                if (fiField):
                                    for indri, ri in enumerate(curve.row_numbers4curve):
                                        x[indri] = self.db_max.iloc[ri, colx]

                            dict = curve.oneCurve_dict_line_prop_2_ind
                            if (moreThan1F):
                                dict[propName4MultiField] = find 
                            line = plParameters.GetLineProperty_from_dict(curve.oneCurve_dict_line_prop_2_ind)
                            if (plotFill):
                                if (fpri == 1):
                                    line.set_linestyle("dashdot")
                                elif (fpri == 0): 
                                    line.set_linestyle("solid")
                                elif (fpri == 2): 
                                    line.set_linestyle("dashed")

                            line.set_xdata(x)
                            if (self.db_logs[fj] > 0):
                                y = [math.log(abs(yv), self.db_logs[fj]) for yv in y]
                                yaxis_log = True
                            if (fiField and (self.db_logs[colx] > 0)):
                                x = [math.log(abs(xv), self.db_logs[colx]) for xv in x]
                                xaxis_log = True
                            line.set_ydata(y)

                            if (plotFill):
                                if (fpri == 1):  # minimum value
                                    xms = x.copy()
                                    yms = y.copy()
                                elif (fpri == 2):  # maximum value
                                    xMs = x.copy()
                                    yMs = y.copy()

                            if ((len(curve.llabel) > 0) and (fpri == 0)):
                                if (not moreThan1F):
                                    line.set_label(curve.llabel)
                                else:
                                    lab = self.db_colnames[fj] + "," + curve.llabel
                                    line.set_label(lab)

                            ax.add_line(line)
                            for xi in x:
                                if (np.isnan(xi) or np.isinf(xi)):
                                    continue
                                if (xi < xm):
                                    xm = xi
                                if (xi > xM):
                                    xM = xi
                            for yi in y:
                                if (np.isnan(yi) or np.isinf(yi)):
                                    continue
                                if (yi < ym):
                                    ym = yi
                                if (yi > yM):
                                    yM = yi
                        if (plotFill):
                            ax.fill_between(xms, yms, yMs, alpha=0.2)
                        ++cntr_curves

            #fig.show()

            # fig.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9)
            if ((xM >= xm) and (not math.isnan(xM)) and (not math.isnan(xm))):
                xmM = max(np.abs(xm), np.abs(xM))
                txmM = max(1e-3 * xmM, 1e-11)
                delx = xM - xm
                if (delx < txmM):
                    delx = txmM
                dx = pad * delx
                xmdx = xm - dx
                xMdx = xM + dx
                plt.xlim([xmdx, xMdx])

            if ((yM >= ym) and (not math.isnan(yM)) and (not math.isnan(ym))):
                ymM = max(np.abs(ym), np.abs(yM))
                tymM = max(1e-3 * ymM, 1e-11)
                dely = yM - ym
                if (dely < tymM):
                    dely = tymM
                dy = pad * dely
                ymdy = ym - dy
                yMdy = yM + dy
                plt.ylim([ymdy, yMdy])

            if (math.isnan(xm) or math.isnan(xM) or math.isnan(ym) or math.isnan(yM)):
                continue

            if (xaxis_log):
                xl = "log(" + fiName + ")"
            else:
                xl = fiName
            plt.xlabel(xl)
            if (yaxis_log):
                yl = "log(" + fjName + ")"
            else:
                yl = fjName
            plt.ylabel(yl)
            if splitInstructions.mid_sz > 0:
                plt.legend()
            # fig.show()

            pltName = plotFullNameWOExt + ".svg"
            plt.savefig(pltName, format='svg')
            pltName = plotFullNameWOExt + ".png"
            plt.savefig(pltName, format='png')
            # pltName = plotFullNameWOExt + ".pdf"
            # plt.savefig(pltName, format='pdf')
#            pltName = plotFullNameWOExt + ".fig.pickle"
#            pickle.dump(fig, open(pltName, 'wb')) 
            fig.clf()
            plt.close()

            cpFolder = schemeFolder + "/" +  "_" + fifjName
            rutil.rmkdir(cpFolder)        
            cpFWOExt = cpFolder + "/" + fnpart
            shutil.copy(plotFullNameWOExt + ".png", cpFWOExt + ".png")
            shutil.copy(plotFullNameWOExt + ".svg", cpFWOExt + ".svg")

def read_csv(root = "data/Characteristics_data", readMainLineMode = 0):
    if (readMainLineMode == -1): # all raw data
        summaryName = root + "/allData/all_wxs.csv"
    elif (readMainLineMode == 0):
        summaryName = root + "/allData/stat_mean.csv"
    elif (readMainLineMode == 1):
        summaryName = root + "/allData/stat_cov.csv"
    elif (readMainLineMode == 2):
        summaryName = root + "/allData/stat_std.csv"
    try:
        pd_data = pd.read_csv(summaryName)
    except IOError:
        print(f"Cannot open file {summaryName}\n")
    return pd_data

def main_function():

    # fd = Frag1D1la()
    # fd.ap = 1e-3
    # fd.ComputeVals()

    #fdc = Frag1Dlap_curves()
    #fdc.Compute(doPlots = 1)
    #return
    # from scipy.interpolate import interp1d
    # from scipy.optimize import root_scalar

    if False:
        ifo = InpFsOuput()
        # ifo.Print()
        # For Ali 
        ifo.Print_AllSerials_AllPara()
        return

    if False:
        a1 = InpF()
        a1.Initialize_InpF(valsAtVert = True, meshp2 = 14, serNo = 0, llc = -3.0, dd2 = 0.2, isPeriodic = True, reduction_sso = 3, meshp2_4Simulation = -1, meshp2_4Output = -1, shape = 2, useOriginalSN_Mesh4Raw = 1)
        # a2 = InpF()
        # a2.Initialize_InpF(valsAtVert = True, meshp2 = 14, serNo = 0, llc = -1.5, dd2 = 0.6, isPeriodic = True, reduction_sso = 3, meshp2_4Simulation = -1, meshp2_4Output = -1, shape = 2, useOriginalSN_Mesh4Raw = 1)
        with open("test1.txt", "w") as fl:
            print(a1.vals,file=fl)
        with open("test1_stat.txt", "w") as fl:
            nms = a1.stats4SimulationFld.Get_Names("inp_sx_")
            print(nms,file=fl)
            stats = a1.GetStatVecVals()
            print(stats,file=fl)
        #with open("test2.txt", "w") as fl:
        #    print(a2.vals,file=fl)
        return

    if False:
        cl = -2.5
        fn = '../../InhomogeneousFiles/cl' + str(cl) + '_np16385/initial_values_0.txt'
        vals = np.loadtxt(fn)[1:]
        sov = StatOfVec()
        sov.Compute_Vec_Stat(vals)
        # Test_HD_Cor()
        return

    # sov = StatOfVec()
    # sov.Test_Gumbel2Weibull()
    # return

    mo_frac_rate_cf = 0
    mo_frac_rate_f  = 1
    mo_frac_rate_f_wShape = 2
    mo_frac_res_x = 3 # coarsening data for which the energies don't match for high loading rate - see 3 below
    mo_frac_res_x_w_delc_fact = 4 # new data 7/23 that includes factoring deltaC to get the correct energy

    mainOption = mo_frac_rate_f
    rateStudy = ((mainOption == mo_frac_rate_cf) or (mainOption == mo_frac_rate_f) or (mainOption == mo_frac_rate_f_wShape))


    #README:
    # 1. set the path to "InhomogeneousFiles" folder
    InpF.setInputMeshRootFolder("../InhomogeneousFiles")
    # 2. set path to "2023_03_24". I put it inside a "data" folder for myself. You can adjust this path
    folderSource = ''
    if (mainOption == mo_frac_rate_cf):
        folderSource = "../../data/2023_04_11"
    if (mainOption == mo_frac_rate_f):    
        folderSource = "../../data/axt_orig"
    if (mainOption == mo_frac_rate_f_wShape):    
        folderSource = "../../data/axt_shape"
    if (mainOption == mo_frac_res_x):    
        folderSource = "../../data/resolution_x_fracture_scalars"
    if (mainOption == mo_frac_res_x_w_delc_fact):    
        folderSource = "../../data/resolution_x_fracture_scalars_w_coarsening"
    # 3. If needed, adjust output path were the files are generated
    folderDest = "../../data/Characteristics_data"
    generatePlots = True
    cd = Characteristics_data()

    # scatter plots
    if (False):
        cd.Main_Plot_Scatter(folderSource, folderDest)
        return
 
    cd.Main_Characteristics_data(folderSource, folderDest)
    if (not generatePlots):
        return

    # root = "data/2023_03_13_x_resolution_F/_PPS3/"
    # root = "data/2023_03_20/_PPS3/"
    # option 0 # -> default
    # option 10 # la axis, delc in plot lines
    option = 10 #10 #0 #100 # -> la, ldelc out, dd2 mid, 1 lc axis, dd2 middle, 2 lc axis, la midle
    # options >= 100 are going to be used for plotting rawData (runNo is used)

    readMainLineMode = 0  # 0 -> mean, 1 -> cov 2 -> std        | -1 -> rawData rather than stats
    plotFillMode = 0 # 0 -> min, max, 1 -> mean -/+ std
    plotFill = True #sFalse
    if (plotFill):
        readMainLineMode = 0

    inpColDict = {}
    lasym = "la"
    #lasym = "lap" don't uncomment
    lasymXaxis = lasym
    lasymXaxis = "lap"
    if (rateStudy):
        # sortingFields = ["llc", "dd2", "ldelc", lasym]
        if (not plotFill):
            if (mainOption == mo_frac_rate_f_wShape):
                inpColDict["shape"] = "none"
            if (option == 0):
                inpColDict[lasym] = "fillstyle" #"c" # color
                inpColDict["llc"] = "c" # line style
                inpColDict["ldelc"] = "ls" #"fillstyle" # marker style
            #    inpColDict["ssoFS"] = "mfc" #"mfc" # marker face color (not-filled, filled, other colors)
                inpColDict["dd2"] = "marker" # marker fill color
            elif (option == 10):
                inpColDict[lasym] = "fillstyle" #"c" # color
                inpColDict["ldelc"] = "c" # line style
                inpColDict["llc"] = "ls" #"fillstyle" # marker style
            #    inpColDict["ssoFS"] = "mfc" #"mfc" # marker face color (not-filled, filled, other colors)
                inpColDict["dd2"] = "marker" # marker fill color
            elif (option == 1):
                inpColDict[lasym] = "marker" #"c" # color
                inpColDict["llc"] = "fillstyle" # line style
                inpColDict["ldelc"] = "ls" #"fillstyle" # marker style
                inpColDict["dd2"] = "c" # marker fill color
            elif (option == 2):
                inpColDict[lasym] = "c" #"c" # color
                inpColDict["llc"] = "fillstyle" # line style
                inpColDict["ldelc"] = "ls" #"fillstyle" # marker style
                inpColDict["dd2"] = "marker" # marker fill color
            elif (option == 3): # for shape, with shape inside
                if (mainOption == mo_frac_rate_f_wShape):
                    inpColDict["shape"] = "c"
                    inpColDict[lasym] = "fillstyle" #"c" # color
                    inpColDict["llc"] = "c" # line style
                    inpColDict["ldelc"] = "ls" #"fillstyle" # marker style
                    inpColDict["dd2"] = "marker" # marker fill color
        else:
            inpColDict[lasym] = "ls" #"c" # color
            inpColDict["llc"] = "fillstyle" # line style
            inpColDict["ldelc"] = "c" #"fillstyle" # marker style
        #    inpColDict["ssoFS"] = "mfc" #"mfc" # marker face color (not-filled, filled, other colors)
            inpColDict["dd2"] = "marker" # marker fill color
    elif ((mainOption == mo_frac_res_x) or (mainOption == mo_frac_res_x_w_delc_fact)):
        # sortingFields = ["dt_resap", "ssoFS", lasym, "llc", "resolutionFactor"]
        if (option < 100): # stat files
            if (option == 0):
                inpColDict["dt_resap"] = "fillstyle"
                inpColDict["resolutionFactor"] = "marker"
                inpColDict["ssoFS"] = "c"
                inpColDict["llc"] = "ls"
                inpColDict[lasym] = "marker"
                if (mainOption == mo_frac_res_x_w_delc_fact):
                    inpColDict["delc_Tf_fact"] = "c"
                    inpColDict["dt_resap"] = "fillstyle"
                    inpColDict["resolutionFactor"] = "none"
                    inpColDict["ssoFS"] = "ls"
                    inpColDict["llc"] = "none"
                    inpColDict[lasym] = "marker"
            elif (option == 1):
                inpColDict["dt_resap"] = "fillstyle"
                inpColDict["resolutionFactor"] = "c"
                inpColDict["ssoFS"] = "ls"
                inpColDict["llc"] = "marker" # marker fill color
                inpColDict[lasym] = "none"
                if (mainOption == mo_frac_res_x_w_delc_fact):
                    inpColDict["delc_Tf_fact"] = "marker"
                    inpColDict["llc"] = "none"
                    inpColDict["ssoFS"] = "none"
        else: # raw data
            inpColDict["runNo"] = "c"
            inpColDict["dt_resap"] = "fillstyle"
            inpColDict["resolutionFactor"] = "fillstyle"
            inpColDict["ssoFS"] = "ls"
            inpColDict["llc"] = "lw" # marker fill color
            inpColDict[lasym] = "marker"
            if (mainOption == mo_frac_res_x_w_delc_fact):
                inpColDict["delc_Tf_fact"] = "marker"
                inpColDict[lasym] = "none"
                inpColDict["ssoFS"] = "none"

    # trying to see if "runNo" is one of the entries   -> dealing with raw data
    if "runNo" in inpColDict:
        readMainLineMode = -1
        plotFill = False

    pd_data = read_csv(folderDest, readMainLineMode)

    dist_logs = {}
    
    """ Not needed anymore as these fields are added to the dataset
    dist_logs["time_F"] = 10
    dist_logs["phi_d_tFin"] = 10
    dist_logs["phi_d_tFin"] = 10
    dist_logs["phi_d_norm_phi0_tFin"] = 10
    dist_logs["phi_d_norm_phi_input_F"] = 10
    dist_logs["psi_f"] = 10
    dist_logs["psi_f_norm_psiC"] = 10
    dist_logs["lbar_F"] = 10
    """
    dbs = DataBaseSplits()
    dbs.Initialize_DataBaseSplits(pd_data, inpColDict, dist_logs, folderDest, plotFillMode)

    splitInstructions = PlotOMISplit_Instruction()
#    splitInstructions.out_col_nameBases.append("ssoFS")
    if (rateStudy):
        # sortingFields = ["llc", "dd2", "ldelc", lasym]
        if (not plotFill):
            # default option
            if ((mainOption == mo_frac_rate_f_wShape) and (option < 3)):
                splitInstructions.out_col_nameBases.append("shape")
            if (option == 0):
                splitInstructions.out_col_nameBases.append("dd2")
                splitInstructions.out_col_nameBases.append("ldelc")
                #   splitInstructions.out_col_nameBases.append(lasym)
                #   splitInstructions.out_col_nameBases.append("llc")
                splitInstructions.in_col_nameBase = lasym
            elif (option == 10):
                splitInstructions.out_col_nameBases.append("dd2")
                splitInstructions.out_col_nameBases.append("llc")
                #   splitInstructions.out_col_nameBases.append(lasym)
                #   splitInstructions.out_col_nameBases.append("ldelc")
                splitInstructions.in_col_nameBase = lasym
            elif (option == 1):
                splitInstructions.out_col_nameBases.append(lasym)
                splitInstructions.out_col_nameBases.append("ldelc")
                splitInstructions.in_col_nameBase = "llc"
            elif (option == 2):
                splitInstructions.out_col_nameBases.append("dd2")
                splitInstructions.out_col_nameBases.append("ldelc")
                splitInstructions.in_col_nameBase = "llc"
            elif (option == 3):
                splitInstructions.out_col_nameBases.append("dd2")
                splitInstructions.out_col_nameBases.append("ldelc")
                splitInstructions.out_col_nameBases.append("llc")
                splitInstructions.in_col_nameBase = lasym
        else:
            splitInstructions.out_col_nameBases.append("dd2")
            #   splitInstructions.out_col_nameBases.append("ldelc")
            #   splitInstructions.out_col_nameBases.append(lasym)
            splitInstructions.out_col_nameBases.append("llc")
            splitInstructions.in_col_nameBase = lasym
    elif ((mainOption == mo_frac_res_x) or (mainOption == mo_frac_res_x_w_delc_fact)):
        # sortingFields = ["dt_resap", "ssoFS", lasym, "llc", "resolutionFactor"]
        if (option < 100): # stat files
            if (option == 0):
                splitInstructions.out_col_nameBases.append("dt_resap")
                if (mainOption == mo_frac_res_x_w_delc_fact):
                    splitInstructions.out_col_nameBases.append("ssoFS")
                splitInstructions.out_col_nameBases.append("llc")
                splitInstructions.out_col_nameBases.append(lasym)
                splitInstructions.in_col_nameBase = "resolutionFactor"
            elif (option == 1):
                if (mainOption == mo_frac_res_x_w_delc_fact):
                    splitInstructions.out_col_nameBases.append("resolutionFactor")
                splitInstructions.out_col_nameBases.append("ssoFS")
                splitInstructions.out_col_nameBases.append("dt_resap")
                splitInstructions.out_col_nameBases.append("llc")
                splitInstructions.in_col_nameBase = lasym
        else: # raw data
            if (mainOption == mo_frac_res_x_w_delc_fact):
                splitInstructions.out_col_nameBases.append("delc_Tf_fact")
            splitInstructions.out_col_nameBases.append("ssoFS")
            splitInstructions.out_col_nameBases.append("dt_resap") # can comment this out
            splitInstructions.out_col_nameBases.append("llc")
            splitInstructions.out_col_nameBases.append(lasym)
            splitInstructions.in_col_nameBase = "resolutionFactor" #"runNo"

#    splitInstructions.in_col_nameBase = "resolutionFactor"

    # want to make sure it can figure it out automatically
    # splitInstructions.mid_col_nameBases.append("dd2")

    allPlots, splitInstructions = dbs.Generate_AllPlots_fromSplitInstructions(splitInstructions)

    # plott all values versus inner parameter
    if True:
        fj = -1 # all output fields
        fi = -1 # inner loop
        addExact = False
        if ((splitInstructions.in_col_nameBase == "la") and (lasymXaxis == "lap")):
            fi = dbs.db.columns.get_loc("inp_s_lap")
            addExact = True
        dbs.Plot_Scalar_Per_Row_Fi_vs_Fj(splitInstructions, allPlots, fj, fi, plotFill, addExact)

    # plotting mean(min(strenth)) (or mean(strength), cov(strength)), etc. where the second operator is over space
    # versus all fields
    # 22 -> inp_sfx_sim_mean
    # 23 -> inp_sfx_sim_mn
    # 25 -> inp_sfx_sim_cov
    if False:
        fj = -1 # all output fields
        fi = 22 # inner loop
        dbs.Plot_Scalar_Per_Row_Fi_vs_Fj(splitInstructions, allPlots, fj, fi, plotFill)

    # example of plotting one field version the other: 
    if False:
        fi = 40 #phid
        fj = 56 # sigmaM
        dbs.Plot_Scalar_Per_Row_Fi_vs_Fj(splitInstructions, allPlots, fj, fi)

    # example of running multiple fields as y axis in a plot
    # last argument is how to differentiate lines for different fields
    if False:
        fi = [-1]
        fj = [21, 22] # time of max and 0
        dbs.Plot_Scalar_Per_Row_Fi_vs_Fj(splitInstructions, allPlots, fj, fi, "ls")

    if False:
        pinp = input_field.InpF()
        pinp.Initialize_InpF(valsAtVert = True, meshp2 = 14, serNo = 0, llc = -4, dd2 = 0.7, \
                            isPeriodic = True, reduction_sso = 3, \
                            meshp2_4Simulation = 12, meshp2_4Output = 10)
    #                            meshp2_4Stat = -1, meshp2_4Output = -1)
        pre_name = "inp_s_"
        addRaw = True
        names = input_field.InpF.GetStatNames(pre_name, addRaw)
        vals = pinp.GetStatVecVals(addRaw)
        print(names)
        print(vals)

    # sample random mesh realizations
    if (False):
        a1 = InpF()
        a1.Initialize_InpF(valsAtVert = True, meshp2 = 14, serNo = 0, llc = -3.0, dd2 = 0.2, isPeriodic = True, reduction_sso = 3, meshp2_4Simulation = -1, meshp2_4Output = -1)
        a2 = InpF()
        a2.Initialize_InpF(valsAtVert = True, meshp2 = 14, serNo = 0, llc = -1.5, dd2 = 0.6, isPeriodic = True, reduction_sso = 3, meshp2_4Simulation = -1, meshp2_4Output = -1)
        with open("test1.txt", "w") as fl:
            print(a1.vals,file=fl)
        with open("test2.txt", "w") as fl:
            print(a2.vals,file=fl)
        return


def main():
    #test()
    main_function()
if __name__=="__main__":
    main()