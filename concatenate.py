import numpy as np
import sys
import argparse
import os

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	parser.add_argument("--folder", type=str, help="the path of folder where data locates, default is './'", default="./")
	parser.add_argument("tags", type=str, nargs=2, help="Tags of input file 1 and file 2")
	args = parser.parse_args()

	tag_1 = args.tags[0]
	tag_2 = args.tags[1]
	file_1 = "md_analytics_%s.npy" % tag_1
	file_2 = "md_analytics_%s.npy" % tag_2
	d1 = np.load(os.path.join(args.folder, file_1))[()]
	d2 = np.load(os.path.join(args.folder, file_2))[()]

	re = {}
	for suff in ['kinetic_vapor', 'kinetic_water', 'vapor_number', 'water_gyration', 'water_center']:
		dx_1 = d1[suff]
		dx_2 = d2[suff]
		v = dict(dx_1, **dx_2)
		re[suff] = v

	newf = 'md_analytics_%s%s.npy' % (tag_1, tag_2)
	newf = os.path.join(args.folder, newf)
	np.save(newf, re)


