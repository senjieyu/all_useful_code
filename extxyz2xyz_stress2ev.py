from ase.io import iread,write,read
import numpy as np

def stress2ev(data,force_threshold):
    atoms = iread(data)
    temp = []
    count = 0
    for a in atoms:
        f = a.get_forces()
        max_force = (f ** 2).sum(axis=1).max() ** (1 / 2)
        if max_force < force_threshold:
            s = a.get_stress()
            six2nine = np.array([s[0], s[5], s[4], s[5], s[1], s[3], s[4], s[3], s[2]])
            virial = six2nine * -1 * a.get_volume()
            a.info['virial'] = virial
            temp.append(a)
        else:
            count +=1
    print(count)
    return temp

def merge_two(data_1,data_2,name):
    threshold = 100
    new_1 = stress2ev(data_1,threshold)
    new_2 = stress2ev(data_2,threshold)
    total = new_1 + new_2
    write(name,total,format='extxyz')

if __name__ == '__main__':
    # data_2 = 'wt1912.xyz'
    # data_1 ='saltwtNaCl1203.xyz'
    # name = 'NaCl1203_H2O.xyz'

    # data_1 ='saltwtKCl1812.xyz'
    # name = 'KCl1812_H2O.xyz'
    import sys
    input = sys.argv[1]
    output = sys.argv[2]
    #merge_two(data_1, data_2, name)
    total = stress2ev(data=input, force_threshold=100)

    write(output, total, format='extxyz')
