import MDAnalysis as mda
import os
import sys
import copy
import numpy as np

# Kalvin temperature : T = m<v^2>/3kB

# box size before enlarger
prev_boxl = 75.046



if __name__ == '__main__':
	
	trrfile = sys.argv[1]
	tprfile = sys.argv[2]
	try:
		mole_r = float(sys.argv[3])  # in angstrom
	except:
		mole_r = (prev_boxl**3*3/4/np.pi)**(1.0/3.0)

	this = mda.Universe(tprfile, trrfile)
	CA = this.select_atoms("name CA")
	center = {}
	print("Analysing centers ...")
	for ct in this.trajectory:
		i = ct.frame
		pos = CA.positions
		center[i] = np.mean(pos, axis=0)
		# output
		sys.stdout.write("\rProcessing %d frames" %i)
		sys.stdout.flush()
	print("\n")

	water = this.select_atoms("name OW")
	box = this.dimensions     # [a,b,c,al,be,ga], in angstrom

	vapor_vel = []
	vapor_ind = []
	vapor_num_frame = {}
	E_water_frame = {}
	E_vapor_frame = {}
	water_radius = {}

	print("Analysing evaporation and kinetic ...")
	for wt in this.trajectory:
		i = wt.frame
		if i>2000:
			break
		# get vapor
		pos = water.positions
		vel = water.velocities
		n_water = len(vel)
		E_all = np.sum(np.linalg.norm(vel, axis=1)**2)
		R = np.linalg.norm(pos-center[i], axis=1)
		candidate = np.where(R > mole_r+20)[0]
		prev_vapor_index = copy.copy(vapor_ind)
		vapor_ind = list(set(prev_vapor_index) | set(candidate))
		vapor_vel = vel[vapor_ind,:]
		# energy
		vapor_num_frame[i] = len(vapor_ind)
		tmp = np.sum(np.linalg.norm(vapor_vel, axis=1)**2)
		E_water_frame[i] = (E_all-tmp)/(n_water-len(vapor_ind))
		if len(vapor_ind) == 0:
			E_vapor_frame[i] = 0
		else:
			E_vapor_frame[i] = tmp / len(vapor_ind)
		# gyration
		water_radius[i] = (np.sum(R)-np.sum(R[vapor_ind])) / (n_water-len(vapor_ind))

		# output
		sys.stdout.write("\rProcessing %d frames" %i)
		sys.stdout.flush()
	print("\nDone.")

	np.save('vapor_number.npy', vapor_num_frame)
	np.save('kinetic_water.npy', E_water_frame)
	np.save('kinetic_vapor.npy', E_vapor_frame)
	np.save('water_gyration.npy', water_radius)
