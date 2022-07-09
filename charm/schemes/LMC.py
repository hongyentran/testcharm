''' LMC scheme
'''

from functools import reduce

from charm.toolbox.pairinggroup import PairingGroup, ZR, G1, G2, pair
import numpy as np
import time
import sys
from math import log

debug = False


def bytes_needed(n):
    if n == 0:
        return 1
    return int(log(n, 256)) + 1


class LMC:
    def __init__(self, groupObj):
        global group
        group = groupObj

    def setup(self, q, l, secparam=None):
        self.q = q
        self.l = l
        alpha = group.random(ZR)
        g = group.random(G1)

        z = [0] * (q+1)
        for i in range(1, q+1):
            z[i] = group.random(ZR)

        G = [g ** 0] * (l+1)
        for j in range(1, l+1):
            G[j] = g ** (alpha ** j)

        H = [[g ** 0] * (2 * l + 1)] * (q + 1)
        for i in range(1, q+1):
            for j in range(1, 2*l+1):
                H[i][j] = g ** (z[i] * (alpha ** j))

        self.sk = {"alpha": alpha, "z": z}
        self.pp = {'g': g, 'G': G, 'H': H}

    def com(self, x):
        C = 1
        for j in range(1, self.l+1):
            C = C * self.pp['G'][j] ** x[j]
        return C

    def open(self, F, y, aux):
        witness = self.pp['g'] ** 0
        for i in range(1, self.q+1):
            for j in range(1, self.l+1):
                for k in range(1,self.l+1):
                    if k != j:
                        witness = witness * (self.pp['H'][i][self.l+1+j-k] ** (F[i][k] * aux[j]))
        return witness

    def openfast(self, F, y, aux):
        witness = self.pp['g'] ** 0
        s = [[0] * (2 * self.l)] * (q+1)

        for i in range(1, self.q+1):
            for delta in range(1-self.l, self.l):
                s[i][delta + self.l] = 0
                for j in range(1, self.l+1):
                    if 0 < j - delta < self.l+1:
                        #try:
                        s[i][delta + self.l] += F[i][j - delta] * aux[j]
                        #except IndexError:
                        #    print('i, j, delta=', i, j, delta)

        for i in range(1, self.q+1):
            for delta in range(1-self.l, self.l):
                if delta != 0:
                    witness *= self.pp['H'][i][self.l + 1 + delta] ** s[i][delta + self.l]

        return witness

    def verify(self, C, F, y, witness):
        u = self.pp['g'] ** 0
        for i in range(1, self.q+1):
            for j in range(1, self.l+1):
                u = u * (self.pp['H'][i][self.l+1-j] ** F[i][j])
        v = self.pp['g'] ** 0
        for i in range(self.q):
            v = v * (self.pp['H'][i][self.l] ** y[i])

        if pair(C, u) == pair(self.pp['G'][1], v) * pair(witness, self.pp['g']):
            print("Fx = y")
        else:
            print("Fx != y")


if __name__ == "__main__":
    debug = True
    q = 1
    l = 128

    print('l = ', l)

    group = PairingGroup('SS512')
    lmc = LMC(group)
    print('Byte size of group order: ', bytes_needed(group.order()))

    start = time.time()
    lmc.setup(q, l)
    print('Set up time:', time.time() - start)
    x = [0] * (l+1)
    for i in range(1, l+1):
        x[i] = group.random(ZR)
    F = [[0] * (l+1)] * (q+1)
    # print(x)
    # a1 = 141521963257168023000607995328529515396172492561
    # a2 = 575220683234108754582923420788674310618729045392
    #
    #
    # t = time.time()
    # for i in range(100000):
    #     a = x[0] ** x[1]
    # print("Time GMP: ", time.time() - t)

    for i in range(1, q+1):
        for j in range(1, l+1):
            F[i][j] = group.random(ZR)

    y = np.dot(np.array(F), np.array(x))

    start = time.time()
    C = lmc.com(x)
    print('Compute commitment time:', time.time() - start)

    start = time.time()
    witness = lmc.open(F, y, x)
    print('Generate proof time:', time.time() - start)

    start = time.time()
    witness = lmc.openfast(F, y, x)
    print('Fast Generate proof time:', time.time() - start)

    start = time.time()
    lmc.verify(C, F, y, witness)
    print('Verify time:', time.time() - start)

    y = y + group.random(ZR)
    lmc.verify(C, F, y, witness)




