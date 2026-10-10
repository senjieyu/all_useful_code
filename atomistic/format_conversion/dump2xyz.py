from ase.io import iread,write,read
import glob
import sys

def dump2xyz(dump,out_xyz,ele):
    data = list(iread(dump))
    map_dic = {}

    for i in range(len(ele)):
        map_dic.update({i + 1: ele[i]})
    temp = []
    for a in data:
        an = a.get_atomic_numbers()
        new_ele = [map_dic[i] for i in an]
        a.set_chemical_symbols(new_ele)
        temp.append(a)
    write(out_xyz, temp, format='extxyz')
    print(f'stru_num: {len(temp)}')

def dump2last_atom(dump,ele):
    a = list(iread(dump))[-1]
    map_dic = {}

    for i in range(len(ele)):
        map_dic.update({i + 1: ele[i]})

    an = a.get_atomic_numbers()
    new_ele = [map_dic[i] for i in an]
    a.set_chemical_symbols(new_ele)
    return a


if __name__ == '__main__':
    ele = ['Cu']
    #dump = 'nvt.dump'
    #out_xyz = 'mlp_md.xyz'
    dump = sys.argv[1]
    out_xyz = sys.argv[2]
    print(f'dump_file:{dump}, out_xyz:{out_xyz}')
    dump2xyz(dump, out_xyz,ele)

    # dump_list = glob.glob('dump*')
    # total_list = []
    # for dump in dump_list:
    #     total_list.append(read(dump))
    #     #print(read(dump))
    # write('lmp_neb.xyz',total_list,format='extxyz')
    # atoms = list(iread('dump.neb.1'))
    # write ('cal.xyz',atoms,format='extxyz')
