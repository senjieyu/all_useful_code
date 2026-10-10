import numpy as np
from ase.io import read,write,iread
from ase.build import make_supercell
import os
import sys

def additional_rand(atom,std):
    a=atom.arrays['positions']
    s=atom.get_chemical_symbols()
    rng = np.random.RandomState(1234)
    for i in range(len(s)):
        if s[i] in ["Li",'Na']:
            a[i]=a[i]+rng.normal(scale=std, size=3)
    atom.set_positions(a)

if __name__ == '__main__':
    in_files = sys.argv[1:-1]
    out_file = sys.argv[-1]
    #os.system(f"rm {out_file}")
    sc = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

    print(f'in_files:{in_files} out_file:{out_file}')

    total = []

    for f in in_files:
        aa=read(f)
        a=make_supercell(aa,sc)
        spos=a.get_scaled_positions()
        cell=a.get_cell()
        strain=[1.05,1.03,1.0,0.97,0.95]
        stru=[]
        for s in strain:
            _a=a.copy()
            _a.set_cell(cell*s)
            _a.set_scaled_positions(spos)
            stru.append(_a)
        for a in stru:
            total.append(a)
            #write(out_file,a,append=True,format="extxyz")
            for i in range(8):
                a_=a.copy()
                a_.rattle((i+1)*0.02,rng=np.random)
                #additional_rand(a_,5*(i+1)*0.01)
                total.append(a_)
                #write(out_file,a_,append=True,format="extxyz")
    write(out_file, total, format="extxyz")
    print(len(list(iread(out_file))))
    

