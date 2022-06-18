# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

from charm.schemes.commit.commit_gs08 import Commitment_GS08
from charm.toolbox.pairinggroup import PairingGroup, G1, G2, pair, ZR # ZR?
from charm.core.engine.util import objectToBytes
import gmpy2
import random
debug = True


def testCommitment_GS08():
    groupObj = PairingGroup('SS512')
    cm = Commitment_GS08(groupObj)

    pk = cm.setup()
    if debug:
        print("Public parameters...")
        print("pk =>", pk)

    m = groupObj.random(G1)
    if debug: print("Committing to =>", m)
    (c, d) = cm.commit(pk, m)

    assert cm.decommit(pk, c, d, m), "FAILED to decommit"
    if debug: print("Successful and Verified decommitment!!!")

def testPairing():
    # Test BLS signature schemes from pairings
    # Generate (pk, sk)
    group = PairingGroup('MNT224')  # 'MNT224' represents an asymmetric EC with 224-bit base field
    sk = group.random(ZR)   # generate a random element in ZR, r is group order (r = group.order())
    g = group.random(G2)    # generate a random element in G2
    pk = g ** sk    # cannot use gmpy2.powmod(g, sk, q)

    # Sign a message using sk
    m = 123
    M = objectToBytes(m, group)     # encode m
    sig = group.hash(M, G1) ** sk

    # Verify a message using pk
    m2 = 123
    M2 = objectToBytes(m2, group)
    h = group.hash(M2, G1)
    if pair(g, sig) == pair(h, pk):
        print("Accept!")
    else:
        print("Not accept!")




# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    testPairing()
    # testCommitment_GS08()


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
