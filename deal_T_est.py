import numpy as np
import os
import sys
import h5py

def read_topol(topol_file):
	num = 0
	with open(topol_file, 'r') as fp:
		lines = fp.readlines()
	for line in lines:
		if line.startswith("SOL"):
			info = line.strip('\n').strip().split()
			num = int(info[-1])
			break
	return num

if __name__ == "__main__":
	frames = int(sys.argv[1])

	v2_cluster = []
	v2_vapor = []
	num_mol = []
	gr_cluster = []
	total_num = -1

	for i in range(frames):
		ii = i+1
		if i==0:
			total_num = read_topol('../topol/topol_%03d.top' % ii)
		num_prev = read_topol('../topol/topol_%03d.top' % ii)
		num_this = read_topol('../topol/topol_%03d.top' % (ii+1))
		num_mol.append(total_num-num_this)
		d = np.load('./T-estimate/mean_v2_%03d.npy' % ii)[()]
		v2_cluster.append(d['mean_v2'][-1])
		vapor_tmp = d['mean_v2_vap'][-1]
		if not np.isnan(vapor_tmp):
			if len(v2_vapor) == 0:
				v2_vapor.append(vapor_tmp*(num_prev-num_this))
			else:
				v2_vapor.append(v2_vapor[-1]+vapor_tmp*(num_prev-num_this))
		else:
			v2_vapor.append(v2_vapor[-1])
		gr_cluster.append(d['cluster_gr'][-1])

	v2_vapor_ = np.array(v2_vapor) / np.array(num_mol)
	fp = h5py.File('T-estimate-%d.h5' % frames, 'w')
	fp.create_dataset('vapor_k', data=v2_vapor_)
	fp.create_dataset('cluster_k', data=v2_cluster)
	fp.create_dataset('vapor_num', data=num_mol)
	fp.create_dataset('cluster_gr', data=gr_cluster)
	fp.close()
