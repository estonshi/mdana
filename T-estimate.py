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
    parser.add_argument("--save", type=str, help="give a file (.npy) to save result", default="none")
    args = parser.parse_args()

    tprfile = args.tpr
    trrfile = args.trr
    this = mda.Universe(tprfile, trrfile)
    CA = this.select_atoms("name CA")
    water = this.select_atoms("name OW")
    box = this.dimensions
    vapor_ind = []
    vapor_vel = []
    mean_V2 = []
    mean_V2_vap = []
    cluster_center = []
    cluster_gr = []

    for wt in this.trajectory:
        i = wt.frame
        t = wt.time
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
	    center_clust = np.mean(pos, axis=0)
            R = np.linalg.norm(pos-center_clust, axis=1)
        else:
            candidate = R_sort_index[cutoff[0]+1:]
	    non_candidate = R_sort_index[:cutoff[0]+1]
	    center_clust = np.mean(pos[non_candidate], axis=0)
            R = np.linalg.norm(pos[non_candidate]-center_clust, axis=1)
        R_sort = np.sort(R)
        gr = np.mean(R_sort[-100:])

        # center
        cluster_center.append(center_clust)
        # gr
        cluster_gr.append(gr)
        # mean(v2)
        vapor_ind = list(set(candidate))
        vapor_vel = vel[list(set(candidate))]
        mv2 = (np.sum(vel**2)-np.sum(vapor_vel**2))/(n_water-len(vapor_ind))
        mean_V2.append(mv2)
        mean_V2_vap.append(np.mean(vapor_vel**2))

    np.save(args.save, {"mean_v2":mean_V2, "mean_v2_vap":mean_V2_vap, "cluster_center":cluster_center, "cluster_gr":cluster_gr})
