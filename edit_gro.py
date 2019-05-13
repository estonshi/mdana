import numpy as np
import sys
import os
import argparse
import re

def parse_coor(line):
    re = [None] * 10
    re[0] = int(line[:5])
    re[1] = line[5:10].strip()
    re[2] = line[10:15].strip()
    re[3] = int(line[15:20])
    tmp = re.findall(r"\-?\d+\.?\d*", re[20:-24])
    re[4] = float(tmp[0])
    re[5] = float(tmp[1])
    re[6] = float(tmp[2])
    re[7] = float(line[-24:-16])
    re[8] = float(line[-16:-8])
    re[9] = float(line[-8:])
    return re

def read_gro(grofile, ions=['cl']):
    dataset = {}
    info_protein = []
    info_water = []
    info_ion = []
    pos_water = []
    with open(grofile, 'r') as fp:
        lines = fp.readlines()
    for ind, line in enumerate(lines):
        line = line.strip('\n')
        if ind == 0:
            dataset['title'] = line.strip()
            continue
        if ind == 1:
            dataset['atom_num'] = int(line.strip())
            continue
        if ind == 2 + dataset['atom_num']:
            dataset['box'] = line
            continue
        re = parse_coor(line)
        if re[1].strip().lower() != "sol":
            if re[1].strip().lower() in ions:
                info_ion.append(re)
            else:
                info_protein.append(re)
        else:
            info_water.append(re)
            if re[2].strip().lower() == "ow":
                pos_water.append(re[4:7])
    pos_water = np.array(pos_water)
    n_water = len(pos_water)
    # calculate center
    center = np.zeros(3)
    ca_num = 0
    for atom in info_protein:
        if atom[2].strip().lower() in ['ca', 'c']:
            center += np.array(atom[4:7])
            ca_num += 1
    center /= ca_num
    # 
    if n_water != len(info_water)//3:
        return None
    # return dataset
    dataset['center'] = center
    dataset['info_protein'] = info_protein
    dataset['info_water'] = info_water
    dataset['info_ion'] = info_ion
    dataset['pos_water'] = pos_water
    return dataset
    
def write_gro(filename, dataset):
    fmt = "{:5d}{:5}{:>5}{:5d}{:8.3f}{:8.3f}{:8.3f}{:8.4f}{:8.4f}{:8.4f}"
    lines = []
    lines.append(dataset['title']+'\n')
    lines.append(str(dataset['atom_num'])+'\n')
    for re in dataset['info_protein']:
        l = fmt.format(re[0],re[1].strip(),re[2].strip(),re[3],re[4],re[5],re[6],re[7],re[8],re[9])
        lines.append(l + '\n')
    for re in dataset['info_water']:
        l = fmt.format(re[0],re[1].strip(),re[2].strip(),re[3],re[4],re[5],re[6],re[7],re[8],re[9])
        lines.append(l + '\n')
    for re in dataset['info_ion']:
        l = fmt.format(re[0],re[1].strip(),re[2].strip(),re[3],re[4],re[5],re[6],re[7],re[8],re[9])
        lines.append(l + '\n')
    lines.append(dataset['box'])
    with open(filename, 'w') as fp:
        fp.writelines(lines)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gro", help="path of input .gro file", type=str)
    parser.add_argument("--tag", help="tag that appended to the name of output files, default is '0'", type=str, default="0")
    parser.add_argument("--mode", help="edit mode, 'DROPEVA'/...", default="DROPEVA")
    parser.add_argument("--top", help="topol file, default is none", type=str, default="none")
    args = parser.parse_args()

    grofile = args.gro
    if not os.path.isfile(grofile):
        raise RuntimeError("gro file is invalid")
    tag = args.tag
    savefile = os.path.splitext(grofile)[0][0:-3] + "%s" % tag + ".gro"
    mode = args.mode
    if mode not in ['DROPEVA']:
        raise RuntimeError("Don't support mode :%s" % mode)

    D = read_gro(grofile)
    if D is None:
        raise RuntimeError("Fail to read gro file")
    n_water = len(D['pos_water'])

    if mode == "DROPEVA":
        new_water_info = []
        new_ion_info = []
        # calculate eva
        R = np.linalg.norm(D['pos_water']-D['center'], axis=1)
        R_sort = np.sort(R)
        R_sort_index = np.argsort(R)
        R_sort_grad = R_sort[1:] - R_sort[:-1]
        R_grad_max = np.max(R_sort_grad[0:int(n_water*0.2)])*2
        cutoff = np.where(R_sort_grad > R_grad_max)[0]
        if len(cutoff) == 0:
            candidate = np.array([])
        else:
            candidate = R_sort_index[cutoff[0]+1:]
        ######
        print("Dropped SOL mocelule : %d" % len(candidate))
        ######
        # newlines
        mol_idx = D['info_protein'][-1][0]
        atom_idx = D['info_protein'][-1][3]
        for ind, atom in enumerate(D['info_water']):
            if ind//3 in candidate.tolist():
                continue
            atom_idx += 1
            if ind % 3 == 0:
                mol_idx += 1
            atom[0] = mol_idx
            atom[3] = atom_idx
            new_water_info.append(atom)
        D['info_water'] = new_water_info
        for ind, atom in enumerate(D['info_ion']):
            atom_idx += 1
            mol_idx += 1
            atom[0] = mol_idx
            atom[3] = atom_idx
            new_ion_info.append(atom)
        D['info_ion'] = new_ion_info
        D['atom_num'] -= len(candidate)*3
        # write
        print("Edit %s to %s." % (grofile, savefile))
        write_gro(savefile, D)
        # edit topol?
        if os.path.isfile(args.top):
            newtop = os.path.splitext(args.top)[0][0:-3] + '%s' % tag + '.top'
            print("Edit %s to %s." % (args.top, newtop))
            with open(args.top, 'r') as fp:
                lines = fp.readlines()
            lines[-2] = "SOL%13d\n" % (n_water - len(candidate))
            with open(newtop, 'w') as fp:
                fp.writelines(lines)

