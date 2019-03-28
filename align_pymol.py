#! /home/ycshi/software/pymol/bin/python
from pymol import cmd
import sys
import glob
import os
import argparse
    

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_folders", type=str, nargs=2, help="2 folders containing pdb files to compare, use blank to seperate them.\
    Corresponding pdb files in two folders should share the same name.")
    parser.add_argument("--save", type=str, help="save newly generated pdb files to ? e.g. '~/Downloads'")
    args = parser.parse_args()

    d1 = args.pdb_folders[0].rstrip("/")
    d2 = args.pdb_folders[1].rstrip("/")
    dnew = args.save
    save_folder = os.path.join(dnew, os.path.split(d1)[-1]+"-aligned")
    if not os.path.isdir(save_folder):
        os.mkdir(save_folder)

    files1 = glob.glob(os.path.join(d1, "*.pdb"))
    print("There are %d files." % len(files1))
    for ii,f1 in enumerate(files1):
        if not os.path.isfile(d2):
            f2 = os.path.join(d2, os.path.split(f1)[-1])
        else:
            f2 = d2
        cmd.load(f1, "d1")
        cmd.load(f2, "d2")
        cmd.remove("name h* and d1")
        cmd.remove("name h* and d2")
        re = cmd.align("d1 and name c+ca", "d2 and name c+ca")
        cmd.save(os.path.join(save_folder, os.path.split(f1)[-1]), "d1")
        cmd.delete("d1 and d2")
        sys.stdout.write("Processing %d/%d files\r" % (ii, len(files1)))
        sys.stdout.flush()
