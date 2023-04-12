import numpy as np

def outer_product_matrix(directionalSizes):
    order = len(directionalSizes)
    if (order == 0):
        retVal = []
        retVal.append([])
        return retVal, 1, [1]
    if (order == 1):
        sz0 = directionalSizes[0]
        retVal = []
        for i  in range(sz0):
            tmp = [i]
            retVal.append(tmp)
        return retVal, sz0, [sz0]
    if (order == 2):
        sz0 = directionalSizes[0]
        sz1 = directionalSizes[1]
        retVal = []
        for i0  in range(sz0):
            for i1  in range(sz1):
                tmp = [i0, i1]
                retVal.append(tmp)
        return retVal, sz0 * sz1, [sz0, sz1]
    if (order == 3):
        sz0 = directionalSizes[0]
        sz1 = directionalSizes[1]
        sz2 = directionalSizes[2]
        retVal = []
        for i0  in range(sz0):
            for i1  in range(sz1):
                for i2  in range(sz2):
                    tmp = [i0, i1, i2]
                    retVal.append(tmp)
        return retVal, sz0 * sz1 * sz2, [sz0, sz1, sz2]
    if (order == 4):
        sz0 = directionalSizes[0]
        sz1 = directionalSizes[1]
        sz2 = directionalSizes[2]
        sz3 = directionalSizes[3]
        retVal = []
        for i0  in range(sz0):
            for i1  in range(sz1):
                for i2  in range(sz2):
                    for i3  in range(sz3):
                        tmp = [i0, i1, i2, i3]
                        retVal.append(tmp)
        return retVal, sz0 * sz1 * sz2 * sz3, [sz0, sz1, sz2, sz3]
    if (order == 5):
        sz0 = directionalSizes[0]
        sz1 = directionalSizes[1]
        sz2 = directionalSizes[2]
        sz3 = directionalSizes[3]
        sz4 = directionalSizes[4]
        retVal = []
        for i0  in range(sz0):
            for i1  in range(sz1):
                for i2  in range(sz2):
                    for i3  in range(sz3):
                        for i4  in range(sz4):
                            tmp = [i0, i1, i2, i3, i4]
                            retVal.append(tmp)
        return retVal, sz0 * sz1 * sz2 * sz3 * sz4, [sz0, sz1, sz2, sz3, sz4]
    if (order == 6):
        sz0 = directionalSizes[0]
        sz1 = directionalSizes[1]
        sz2 = directionalSizes[2]
        sz3 = directionalSizes[3]
        sz4 = directionalSizes[4]
        sz5 = directionalSizes[5]
        retVal = []
        for i0  in range(sz0):
            for i1  in range(sz1):
                for i2  in range(sz2):
                    for i3  in range(sz3):
                        for i4  in range(sz4):
                            for i5  in range(sz5):
                                tmp = [i0, i1, i2, i3, i4, i5]
                                retVal.append(tmp)
        return retVal, sz0 * sz1 * sz2 * sz3 * sz4 * sz5, [sz0, sz1, sz2, sz3, sz4, sz5]
    print(order)
    raise TypeError("Too high of an order")

# this is a class that enables going back and force between a multiIndex of size "order" and one index, using outer_product_matrix
class SMIndex:
    def __init__(self, directionalSizesIn = []):
        self.directionalSizes = directionalSizesIn
        self.Initialize_MultiIndex(directionalSizesIn)

    def savetxt(self, fl):
        print(f"ordr {(int)(self.order)}\n",file=fl)
        print(f"totalSz {(int)(self.totalSz)}\n",file=fl)
        print("directionalSizes\n",file=fl)
        np.savetxt(fl, self.directionalSizes, delimiter=",")
        print("si2mi\n",file=fl)
        np.savetxt(fl, self.si2mi, fmt='%d', delimiter=",")

    def Initialize_MultiIndex(self, directionalSizesIn):
        self.directionalSizes = directionalSizesIn
        self.order = len(self.directionalSizes)
        self.si2mi, self.totalSz, directionalSzs = outer_product_matrix(self.directionalSizes)
        self.totalDirectionalFactors = np.zeros(self.order, dtype=int)
        if (self.order > 0):
            self.totalDirectionalFactors[self.order - 1] = 1
            for i in range(self.order - 1):
                ind = self.order - i - 2
                self.totalDirectionalFactors[ind] = self.directionalSizes[ind + 1] * self.totalDirectionalFactors[ind + 1]

    # multi index to single index
    def MI2SI(self, multiIndex):
        if (self.order == 0):
            return 0
        ind = 0
        for i in range(self.order):
            ind += self.totalDirectionalFactors[i] * multiIndex[i]
        return (int)(ind)

    # single to multi index
    def SI2MI(self, singleIndex):
        if (self.order == 0):
            return []
        return self.si2mi[singleIndex]
