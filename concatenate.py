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

	for suff in ['kinetic_vapor', 'kinetic_water', 'vapor_number', 'water_gyration', 'water_center']:
		file_1 = suff + '_%s.npy' % tag_1
		file_2 = suff + '_%s.npy' % tag_2
		file_1 = os.path.join(args.folder, file_1)
		file_2 = os.path.join(args.folder, file_2)
		d1 = np.load(file_1)[()]
		d2 = np.load(file_2)[()]
		v = dict(d1, **d2)
		np.save(suff + '_%s%s.npy' % (tag_1, tag_2), v)


