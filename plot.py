import matplotlib.pyplot as plt
import numpy as np
import sys


if __name__ == '__main__':
	try:
		kv = np.load(sys.argv[1])[()]
		kw = np.load(sys.argv[2])[()]
		vn = np.load(sys.argv[3])[()]
		wg = np.load(sys.argv[4])[()]
		ct = np.load(sys.argv[5])[()]
	except:
		tag = sys.argv[1]
		kv = np.load('kinetic_vapor_%s.npy' % tag)[()]
		kw = np.load('kinetic_water_%s.npy' % tag)[()]
		vn = np.load('vapor_number_%s.npy' % tag)[()]
		wg = np.load('water_gyration_%s.npy' % tag)[()]
		ct = np.load('water_center_%s.npy' % tag)[()]

	plt.plot(np.array(kv.keys()[1:])/1000, kv.values()[1:], 'r.')
	plt.plot(np.array(kw.keys()[1:])/1000, kw.values()[1:], 'b.')
	plt.legend(['vapor', 'cluster'])
	plt.title('kinetic energy')
	plt.xlabel('time /ns')
	plt.ylabel('kinetic energy (<v^2>)')
	plt.show()

	plt.plot(np.array(vn.keys()[1:])/1000, vn.values()[1:], 'b.')
	plt.title('Evaporated mocules')
	plt.xlabel('time /ns')
	plt.ylabel('Number')
	plt.show()

	ctv = np.array(ct.values())
	plt.plot(np.array(ct.keys()[1:])/1000, np.array(ctv[1:,0])/10.0, 'b.')
	plt.plot(np.array(ct.keys()[1:])/1000, np.array(ctv[1:,1])/10.0, 'k.')
	plt.plot(np.array(ct.keys()[1:])/1000, np.array(ctv[1:,2])/10.0, 'r.')
	plt.title('Center of water droplet')
	plt.xlabel('time /ns')
	plt.ylabel('Radius /nm')
	plt.legend(['x', 'y', 'z'])
	plt.show()

	plt.errorbar(np.array(wg.keys()[1:])/1000, np.array(wg.values())[1:,0]/10.0,yerr=np.array(wg.values())[1:,1]/10.0,fmt='b.')
	plt.title('cluster gyration')
	plt.xlabel('time /ns')
	plt.ylabel('Radius /nm')
	plt.show()
