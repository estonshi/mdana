import matplotlib.pyplot as plt
import numpy as np
import sys


if __name__ == '__main__':
	try:
		kv = np.load(sys.argv[1])[()]
		kw = np.load(sys.argv[2])[()]
		vn = np.load(sys.argv[3])[()]
		wg = np.load(sys.argv[4])[()]
	except:
		kv = np.load('kinetic_vapor.npy')[()]
		kw = np.load('kinetic_water.npy')[()]
		vn = np.load('vapor_number.npy')[()]
		wg = np.load('water_gyration.npy')[()]

	plt.plot(np.array(kv.keys()[1:])*20.0/1000, kv.values()[1:], 'r-')
	plt.plot(np.array(kw.keys()[1:])*20.0/1000, kw.values()[1:], 'b-')
	plt.legend(['vapor', 'cluster'])
	plt.title('kinetic energy')
	plt.xlabel('time /ns')
	plt.ylabel('kinetic energy (<v^2>)')
	plt.show()

	plt.plot(np.array(vn.keys()[1:])*20.0/1000, vn.values()[1:], 'b-')
	plt.title('Evaporated mocules')
	plt.xlabel('time /ns')
	plt.ylabel('Number')
	plt.show()

	plt.plot(np.array(wg.keys()[1:])*20.0/1000, np.array(wg.values()[1:])/10.0, 'b-')
	plt.title('cluster gyration')
	plt.xlabel('time /ns')
	plt.ylabel('Radius /nm')
	plt.show()