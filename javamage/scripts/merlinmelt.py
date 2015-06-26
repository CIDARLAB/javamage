#-------------------------------------------------------------------------------
# Name:        merlinmelt.py
# Purpose:  return the melting temeperature of the input DNA sequence.
#   Adapted from the Biopython API docs
#   http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html
# Usage:    python melt.py [SEQUENCE1, SEQUENCE2...] check=True strict=True
#           c_seq=None shift=0 nn_table=DNA_NN3 tmm_table=DNA_TMM1
#           imm_table=DNA_IMM1 de_table=DNA_DE1 dnac1=25 dnac2=25
#           selfcomp=False Na=50 K=0 Tris=0 Mg=0 dNTPs=0 saltcorr=5
#
# Author:      mquintin
#
# Created:     19/06/2015
# Copyright:   (c) mquintin 2015
# Licence:     BSD 3-Clause License (http://opensource.org/licenses/BSD-3-Clause)
#-------------------------------------------------------------------------------
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import sys

DNA_IMM1=mt.DNA_IMM1
DNA_NN1=mt.DNA_NN1
DNA_NN2=mt.DNA_NN2
DNA_NN3=mt.DNA_NN3
DNA_NN4=mt.DNA_NN4
DNA_TMM1=mt.DNA_TMM1
DNA_DE1=mt.DNA_DE1

seqarr = []

#default values according to BioPython API
check=True
strict=True
c_seq=None
shift=0
nn_table=DNA_NN3
tmm_table=DNA_TMM1
imm_table=DNA_IMM1
de_table=DNA_DE1
dnac1=25 #uM
dnac2=25
selfcomp=False
Na=50 #mM
K=0
Tris=0
Mg=0
dNTPs=0
saltcorr=5

def main():
    for seq in seqarr:
        s = Seq(seq)
        res = mt.Tm_NN(s,check,strict,c_seq,shift,nn_table,
        tmm_table,imm_table,de_table,dnac1,dnac2,selfcomp,Na,
        K,Tris,Mg,dNTPs,saltcorr)
        print('%0.2f' % res)

if __name__ == '__main__':
    list=sys.argv[1:]
    str=' '.join(list)
    str = str.replace('= ','=')
    str = str.replace(' =','=')
    list=str.split(' ')
    for i in list:
        try:
            exec(i)
        except NameError:
            seqarr.append(i)
            pass
    main()
