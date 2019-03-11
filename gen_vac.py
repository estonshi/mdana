import sys
import os
import argparse
import numpy as np


if __name__ == '__main__':
	
	# parser cmd line params
	parser = argparse.ArgumentParser(description="generate vacuum around system box (in .gro)")
	parser.add_argument("-i", "--input", type=str, help="input .gro file")
	parser.add_argument("-l", "--boxlen", type=float, help="new box length, in nm")
	args = parser.parse_args()

	#
	inputfile = args.input
	outputfile = inputfile.split('.')[0] + "-large." + inputfile.split('.')[1]
	bl = float(args.boxlen)

	# read and write file
	sol = 0
	line_num = 0
	atoms = 0
	
	with open(inputfile, "r") as fp:
		lines = fp.readlines()

	mean_pos = np.array([0,0,0])
	for ind, line in enumerate(lines):
		if ind == 1:
			atoms = float(line.strip("\n"))
			continue
		elif ind < 2 or ind == len(lines)-1:
			continue
		xyz = np.array(line[20:44].split()).astype(float)
		mean_pos = mean_pos + xyz
	mean_pos /= atoms

	ori_box = np.array(lines[-1].split()).astype(float)[0:3]
	shift = bl/2.0 - mean_pos

	fp_new = open(outputfile, "w")
	for ind,line in enumerate(lines):
		if ind < 2:
			fp_new.write(line)
			line_num += 1
			continue
		items = line.split()
		if ind == len(lines)-1:
			newline = "%.5f   %.5f   %.5f\n"%(bl, bl, bl)
		else:
			newline = ""+line[0:20]
			xyz = np.array(line[20:44].split()).astype(float)
			xyz = xyz + shift
			newline = newline + "%8.3f%8.3f%8.3f"%(xyz[0], xyz[1], xyz[2])
			newline = newline + line[44:]
		fp_new.write(newline)
		line_num += 1

	fp_new.close()
