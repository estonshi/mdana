import MDAnalysis as mda
import os
import sys
import copy
import numpy as np
import argparse

import matplotlib.pyplot as plt
# Kalvin temperature : T = m<v^2>/3kB

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	parser.add_argument("-r", "--trr", help="path of input .trr file", type=str)
	parser.add_argument("-p", "--tpr", help="path of input .tpr file", type=str)
	#parser.add_argument("-b", "--boxl", help="length of md box before enlarge, in A", type=float)
	parser.add_argument("-f", "--frame", help="frame range, e.g. '0,end' means from 0 to end", type=str, default="0,end")
	parser.add_argument("--tag", help="tag that appended to the name of output files, default is '0'", type=str, default="0")
	parser.add_argument("--noback", help="do not contain backed molecule.", action="store_true", default=False)
	args = parser.parse_args()

	tprfile = args.tpr
	trrfile = args.trr
	#prev_boxl = args.boxl
	frame_range = args.frame.split(',')
	frame_range[0] = int(frame_range[0])
	if frame_range[1].lower() == "end":
		frame_range[1] = "end"
	else:
		frame_range[1] = int(frame_range[1])
	
	#mole_r = (prev_boxl**3*3/4/np.pi)**(1.0/3.0)
	#print("Molecule radius is %.3f nm" % (mole_r/10.0))
	this = mda.Universe(tprfile, trrfile)
	CA = this.select_atoms("name CA")
	water = this.select_atoms("name OW")
	box = this.dimensions     # [a,b,c,al,be,ga], in angstrom

	vapor_vel = []
	vapor_ind = []
	vapor_num_frame = {}
	E_water_frame = {}
	E_vapor_frame = {}
	water_radius = {}
	water_center = {}

	print("Analysing evaporation and kinetic ...")
	for wt in this.trajectory:
		i = wt.frame
		t = wt.time
		if i < frame_range[0]:
			continue
		if frame_range[1] != "end" and i>frame_range[1]:
			break
		# get vapor
		pos = water.positions
		vel = water.velocities
		n_water = len(vel)
		# get center
		prev_vapor_index = copy.copy(vapor_ind)
		center = (np.sum(pos, axis=0) - np.sum(pos[prev_vapor_index], axis=0)) / (len(pos) - len(prev_vapor_index))
		water_center[t] = center
		# calculate new vapor molecues
		E_all = np.sum(np.linalg.norm(vel, axis=1)**2)
		R = np.linalg.norm(pos-center, axis=1)
		R_sort = np.sort(R)
		R_sort_index = np.argsort(R)
		R_sort_grad = R_sort[1:] - R_sort[:-1]
		R_grad_max = np.max(R_sort_grad[0:int(n_water*0.2)])
		candidate = R_sort_index[np.where(R_sort_grad > R_grad_max)[0]+1]
		if not args.noback:
			# contain backed
			vapor_ind = list(set(candidate))
		else:	
			# not contain backed
			vapor_ind = list(set(prev_vapor_index) | set(candidate))
		vapor_vel = vel[vapor_ind,:]
		# energy
		vapor_num_frame[t] = len(vapor_ind)
		tmp = np.sum(np.linalg.norm(vapor_vel, axis=1)**2)
		E_water_frame[t] = (E_all-tmp)/(n_water-len(vapor_ind))
		if len(vapor_ind) == 0:
			E_vapor_frame[t] = 0
		else:
			E_vapor_frame[t] = tmp / len(vapor_ind)
		# gyration
		if len(candidate) > 0:
			R_tail = R_sort[:-len(candidate)]
		else:
			R_tail = R_sort
		r_mean = np.mean(R_tail[-int(0.1*len(R_tail)):])
		r_std = np.std(R_tail[-int(0.1*len(R_tail)):])
		water_radius[t] = [r_mean, r_std]

		# output
		sys.stdout.write("\r- %d frames are processed." %i)
		sys.stdout.flush()
	print("\nDone.")

	np.save('vapor_number_%s.npy' % args.tag, vapor_num_frame)
	np.save('kinetic_water_%s.npy' % args.tag, E_water_frame)
	np.save('kinetic_vapor_%s.npy' % args.tag, E_vapor_frame)
	np.save('water_gyration_%s.npy' % args.tag, water_radius)
	np.save('water_center_%s.npy' % args.tag, water_center)
