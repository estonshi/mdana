import MDAnalysis as mda
import numpy as np
import sys
import os
import argparse
import matplotlib.pyplot as plt

# maxwell distribution

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--trr", type=str, help="path of trr file")
    parser.add_argument("--tpr", type=str, help="path of tpr file")
    parser.add_argument("--window", type=int, help="frame window to calculate distribution", default=50)
    parser.add_argument("--save", type=str, help="path to save figures/movies, default is none", default="none")
    parser.add_argument("--frame", type=str, help="frame range, e.g 0,1000, default is all frames", default='0,0')
    args = parser.parse_args()

    tprfile = args.tpr
    trrfile = args.trr
    frame_range = list(map(int, args.frame.split(',')))
    if frame_range[1] <= frame_range[0]:
        frame_range[1] = float("inf")
    this = mda.Universe(tprfile, trrfile)
    CA = this.select_atoms("name CA")
    water = this.select_atoms("name OW")
    box = this.dimensions
    vapor_ind = []
    vapor_vel = []

    for wt in this.trajectory:
        i = wt.frame
        t = wt.time
        if i > frame_range[1] or i < frame_range[0]:
            continue
        pos = water.positions
        ca_pos = CA.positions
        vel = np.linalg.norm(water.velocities, 2, axis=1)
        n_water = len(vel)
        # center
        center = np.mean(ca_pos, axis=0)
        # vapor ind
        R = np.linalg.norm(pos-center, axis=1)
        R_sort = np.sort(R)
        R_sort_index = np.argsort(R)
        R_sort_grad = R_sort[1:] - R_sort[:-1]
        R_grad_max = np.max(R_sort_grad[0:int(n_water*0.2)])*2
        cutoff = np.where(R_sort_grad > R_grad_max)[0]
        if len(cutoff) == 0:
        	candidate = np.array([])
        else:
        	candidate = R_sort_index[cutoff[0]+1:]
        
        vapor_ind.extend(list(set(candidate)))
        vapor_vel.extend(vel[list(set(candidate))])
        if (i+1)%args.window == 0:
            print("%d frames are processed" % (i+1))
            vapor_ind = set(vapor_ind)
            cluster_ind = set.difference(set(range(n_water)), vapor_ind)
            vapor_ind = list(vapor_ind)
            cluster_ind = list(cluster_ind)
            # plot
            bins = np.linspace(0,20,20)

            plt.hist(vel, bins=bins, color="r")
            plt.hist(vapor_vel, bins=bins, color="b")
            plt.xlabel("V (nm/ps)")
            plt.ylabel("Molecules")
            plt.ylim([0,2000])
            plt.legend(["Total", "Vapor"])
            #plt.subplot(1,3,3)
            #plt.hist(vel[cluster_ind], bins=bins)
            plt.title("Frame : %d" %(i+1))
            if os.path.isdir(args.save):
                plt.savefig(os.path.join(args.save, "%04d.png" % (i+1)))
                plt.clf()
            else:
                plt.show()
            # refresh
            del vapor_ind[:]
            del vapor_vel[:]