import os
import numpy as np

def rrep(str, repPart):
    if (repPart in str):
        strOut = str.replace(repPart, "")
        return strOut, True
    strOut = str
    return strOut, False

def changeName(str, repParts):
    for repPart in repParts:
        strNew, b = rrep(str, repPart)
        if b:
            return strNew

    strNew = str
    return strNew
        

def rmkdir(folderName):
    if os.path.exists(folderName):
        return
    try:
        os.mkdir(folderName)
        #print(f"Directory '{folderName}' created successfully!")
    except OSError as error:
        print(f"Directory '{folderName}' creation failed. Error: {error}")


# ############################### math functions
def pwl_roots(x, y, relTol = 1e-60):
    roots = []
    sz = len(x)
    if (sz == 0):
        return roots
    mxy = max(abs(y))
    absTol = mxy * relTol
    
    xp = x[0]
    yp = y[0]
    sp = 0
    if (yp < -absTol):
        sp = -1
    elif (yp > absTol): 
        sp = 1
    else:
        roots.append(xp)
    for i in range(1, sz):
        xn = x[i]
        yn = y[i]
        sn = 0
        if (yn < -absTol):
            sn = -1
            if (sp == 1):
                xr = (xp * yn - xn * yp) / (yn - yp)
                roots.append(xr)
        elif (yn > absTol): 
            sn = 1
            if (sp == -1):
                xr = (xp * yn - xn * yp) / (yn - yp)
                roots.append(xr)
        else:
            roots.append(xn)
        xp = xn
        yp = yn
        sp = sn
    return np.array(roots)

def pwl_root_crossing_segmentStats(x, y, relTol = 1e-60):
    roots = pwl_roots(x, y, relTol)
    sz = len(roots)

    szx = len(x)
    if ((sz == 0) or (szx == 0)):
        return sz, np.nan, np.nan, np.nan, np.nan
    roots = roots - x[0]
    roots = np.array(roots)
    lengths = np.concatenate(([roots[0]], roots[1:sz] - roots[0:sz - 1], [x[szx - 1] - roots[sz - 1]]))
    sm = np.sum(lengths)
    #    lengths = x[0]
    # lengths.append(roots[2:sz] - roots[1:sz - 1])
    ave = np.mean(lengths)
    mn = np.min(lengths)
    mx = np.max(lengths)
    sdiv = np.std(lengths)
    return len(lengths), ave, mn, mx, sdiv
