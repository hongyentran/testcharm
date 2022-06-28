''' LMC scheme
'''

from functools import reduce

from charm.toolbox.pairinggroup import PairingGroup, ZR, G1, G2, pair
import numpy as np


debug = False


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
        for i in range(q):
            z[i+1] = group.random(ZR)

        G = [g ** 0] * (l + 1)
        for j in range(l):
            G[j+1] = g ** (alpha ** (j+1))

        H = [[g ** 0] * (2*l+1)] * (q+1)
        for i in range(q):
            for j in range(2*l):
                H[i+1][j+1] = g ** (z[i+1] * (alpha ** (j+1)))

        self.sk = {"alpha": alpha, "z": z}
        self.pp = {'g': g, 'G': G, 'H': H}

    def com(self, x):
        C = 1
        for j in range(self.l):
            C = C * self.pp['G'][j+1] ** x[j]
        return C

    def open(self, F, y, aux):
        evidence = self.pp['g'] ** 0
        for i in range(self.q):
            for j in range(self.l):
                for k in range(self.l):
                    if k != j:
                        evidence = evidence * (self.pp['H'][i+1][self.l + 1 + j - k] ** (F[i][k] * aux[j]))
        return evidence

    def verify(self, C, F, y, evidence):
        u = self.pp['g'] ** 0
        for i in range(self.q):
            for j in range(1, self.l):
                u = u * (self.pp['H'][i+1][self.l - j] ** F[i][j])
        v = self.pp['g'] ** 0
        for i in range(self.q):
            v = v * (self.pp['H'][i+1][self.l] ** y[i])

        if pair(C, u) == pair(self.pp['G'][1], v) * pair(evidence, self.pp['g']):
            print("Fx = y")
        else:
            print("Fx != y")


if __name__ == "__main__":
    debug = True
    q = 1
    l = 3

    group = PairingGroup('SS512')
    lmc = LMC(group)

    lmc.setup(q, l)

    x = [group.random(ZR) for i in range(l)]
    F = [[0] * l] * q
    for i in range(q):
        for j in range(l):
            F[i][j] = group.random(ZR)
    y = np.dot(np.array(F), np.array(x))

    C = lmc.com(x)
    evidence = lmc.open(F, y, x)
    lmc.verify(C, F, y, evidence)
    y = y + 1
    lmc.verify(C, F, y, evidence)




