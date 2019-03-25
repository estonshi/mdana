import numpy as np
import sys
import os
import argparse
from spipy.image import io



def cal_rms(f1, f2, atoms=['c','ca','n','o']):
    pdb1 = io.readpdb_full(f1)
    pdb2 = io.readpdb_full(f2)

    if len(pdb1.keys()) != len(pdb2.keys()):
        raise RuntimeError("Input 2 pdb files should be the same protein !")

    res_rms = {}
    for res_i in pdb1.keys():
        r1 = pdb1[res_i][1]
        r2 = pdb2[res_i][1]
        for aname, ainfo1 in r1.items():
            if aname.lower() not in atoms:
                continue
            ainfo2 = r2[aname]
            if not res_rms.has_key(res_i):
                res_rms[res_i] = [0,0]
            res_rms[res_i][0] += 1
            res_rms[res_i][1] += np.sum((np.array(ainfo2[2:5])-np.array(ainfo1[2:5]))**2)

    rms = np.zeros(len(res_rms.keys()))
    tot = 0
    n_atoms = 0
    for kn, kv in res_rms.items():
        rms[kn-1] = np.sqrt(kv[1]/kv[0])
        tot += kv[1]
        n_atoms += kv[0]

    return rms, np.sqrt(tot/n_atoms)



if __name__ == "__main__":
    import glob
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_folders", type=str, nargs=2, help="2 folders containing pdb files to compare, use blank to seperate them.\
    Corresponding pdb files in two folders should share the same name.")
    parser.add_argument("--frames", type=int, help="number of frames", default=-1)
    parser.add_argument("--dt", type=int, help="time step between two pdbs, picosecond", default=100)
    parser.add_argument("--rmsfile", type=str, help="save/load rms values", default="none")
    parser.add_argument("--plot_only", action="store_true", default=False, help="only plot rms distribution loaded from --rmsfile")
    args = parser.parse_args()

    if args.rmsfile.lower() != "none" and not os.path.isdir(os.path.split(args.rmsfile)[0]):
        raise RuntimeError("Given saving file path is in valid !")
    if args.rmsfile.lower() == "none":
        print("No rms files given.")
        if args.plot_only:
            raise RuntimeError("Please give rms file to plot.")

    if not args.plot_only:
        d1 = args.pdb_folders[0]
        d2 = args.pdb_folders[1]
        if args.frames <=0 :
            files1 = len(glob.glob(os.path.join(d1,"*.pdb")))
            files2 = len(glob.glob(os.path.join(d1,"*.pdb")))
            args.frames = min(files1, files2)
        print("Total frames to process : %d" % args.frames)

        rms = [None] * args.frames
        tot = [0] * args.frames
        for i in np.arange(args.frames):
            f1 = os.path.join(d1, "x%04d.pdb" % (i))
            f2 = os.path.join(d2, "x%04d.pdb" % (i))
            rms[i], tot[i] = cal_rms(f1, f2)
            sys.stdout.write("Processing %d/%d files\r" % (i, args.frames))
            sys.stdout.flush()

        rms = np.array(rms).T # y-axis is resi number, x-axis is time
        if args.rmsfile.lower() != "none":
            np.save(args.rmsfile, [rms, tot])
    
    else:
        rms, tot = np.load(args.rmsfile)

    # plot rms
    ytick = np.arange(rms.shape[0]) + 1
    xtick = [''] * rms.shape[1]
    xtick_tmp = np.arange(0,rms.shape[1]) * args.dt/1000.0
    for i in np.arange(0,rms.shape[1],rms.shape[1]//20):
        xtick[i] = xtick_tmp[i]
    xtick[-1] = xtick_tmp[-1]
    sns.heatmap(rms, xticklabels=xtick, yticklabels=10, cmap="YlGnBu")
    plt.title("RMS per residue (A)")
    plt.xlabel("time /ns")
    plt.ylabel("residue index")
    plt.show()

    plt.plot(xtick_tmp, tot, '-k')
    plt.title("Total RMS")
    plt.xlabel("time /ns")
    plt.ylabel("RMS /A")
    plt.show()

    # residue rms
    rms_all = rms.flatten()
    res_all,_ = np.indices(rms.shape)
    res_all = res_all.flatten()
    fig = plt.figure(figsize=(15,5))
    sns.boxplot(x="residue-index", y="RMS", \
    data=pd.DataFrame({'residue-index':res_all, 'RMS':rms_all}))
    sns.despine(offset=10, trim=True)
    plt.xticks(rotation=90)
    fig.tight_layout()
    plt.show()