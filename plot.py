import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from scipy.optimize import curve_fit

def poly(x,b,c,d,e):
	return b*x+c*x**2+d*x**3+e*x**4

def lnx(x,b,c,d,e):
	return np.log(b*x+1)*c + np.arctan(d*x)*e

def poly_div(x,factor):
	y = np.zeros(x.shape)
	for i, kv in enumerate(factor):
		y = y + kv*x**i*(i+1)
	return y

def lnx_div(x,factor):
	y = factor[0]*factor[1]/(factor[0]*x+1) + factor[2]*factor[3]/(factor[2]**2*x**2+1)
	return y

if __name__ == '__main__':
	try:
		path = sys.argv[1]
		tag = sys.argv[2]
		data = np.load(os.path.join(path, 'md_analytics_%s.npy' % tag))[()]
		kv = data['kinetic_vapor']
		kw = data['kinetic_water']
		vn = data['vapor_number']
		wg = data['water_gyration']
		ct = data['water_center']
	except:
		print("Error occurs.\nUsage: python plot.py <path> <tag>")
		sys.exit(1)

	plt.plot(np.array(kv.keys()[1:])/1000, kv.values()[1:], 'r.')
	plt.plot(np.array(kw.keys()[1:])/1000, kw.values()[1:], 'b.')
	plt.legend(['vapor', 'cluster'])
	plt.title('kinetic Energy')
	plt.xlabel('time (ns)')
	plt.ylabel('kinetic energy (<v^2>) (nm/ps)^2')
	plt.show()

	ax = plt.subplot(111)
	ax.plot(np.array(vn.keys()[1:])/1000, vn.values()[1:], 'r.')
	ax.set_title('Evaporation of Water Mocules')
	ax.set_ylabel('Evaporated Number')
	sort_time = np.arange(0,500)    # in ns
	#sort_time = np.sort(np.array(vn.keys()[1:]))/1000    # in ns
	#popt, pcov = curve_fit(poly, np.array(vn.keys()[1:])/1000, np.array(vn.values()[1:]))
	popt, pcov = curve_fit(lnx, np.array(vn.keys()[1:])/1000, np.array(vn.values()[1:]))
	print("Evaporation description : %s" % popt)
	pp = popt.copy()
	for i in range(len(pp)):
		pp[i] = pp[i] * (i+1)
	pp = pp[::-1].tolist()
	pp.append(0)
	zeros = np.roots(pp)
	print("Root of evaporation rate : %s" % str(zeros))
	for tmp in zeros:
		if np.imag(tmp) == 0 and np.real(tmp)>0:
			sort_time = np.arange(0, np.abs(tmp))
	#fit_evap = poly(sort_time, popt[0], popt[1], popt[2], popt[3])
	fit_evap = lnx(sort_time, popt[0], popt[1], popt[2], popt[3])
	ax.plot(sort_time, fit_evap, '--k', linewidth=3.0)
	#ratio = poly_div(sort_time, popt)
	ratio = lnx_div(sort_time, popt)
	ax2 = ax.twinx()
	ax2.plot(sort_time, ratio, 'm.')
	ax2.set_xlabel('time (ns)')
	ax2.set_ylabel('Evaporation rate (molecules/ns)')
	plt.show()

	ctv = np.array(ct.values())
	plt.plot(np.array(ct.keys()[1:])/1000, np.array(ctv[1:,0])/10.0, 'b.')
	plt.plot(np.array(ct.keys()[1:])/1000, np.array(ctv[1:,1])/10.0, 'k.')
	plt.plot(np.array(ct.keys()[1:])/1000, np.array(ctv[1:,2])/10.0, 'r.')
	plt.title('Center of Cluster')
	plt.xlabel('time (ns)')
	plt.ylabel('Coordinates (nm)')
	plt.legend(['x', 'y', 'z'])
	plt.show()

	plt.errorbar(np.array(wg.keys()[1:])/1000, np.array(wg.values())[1:,0]/10.0,yerr=np.array(wg.values())[1:,1]/10.0,fmt='b.')
	plt.title('Cluster Gyration')
	plt.xlabel('time (ns)')
	plt.ylabel('Radius (nm)')
	plt.show()
