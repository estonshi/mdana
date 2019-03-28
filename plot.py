import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from scipy.optimize import curve_fit

def poly(x,b,c,d,e):
	return b*x+c*x**2+d*x**3+e*x**4

def lnatan(x,b,c,d,e):
	return np.arctan(d*x)*e + np.log(b*x+1)*c

def lnsig(x,b,c,d,e):
	return e*((1-np.exp(-2*d*x))/(1+np.exp(-2*d*x))) + np.log(b*x+1)*c

def poly_div(x,factor):
	y = np.zeros(x.shape)
	for i, kv in enumerate(factor):
		y = y + kv*x**i*(i+1)
	return y

def lnatan_div(x,factor):
	b = factor[0]
	c = factor[1]
	d = factor[2]
	e = factor[3]
	y = d*e/(d**2*x**2+1) + b*c/(b*x+1)
	return y

def lnsig_div(x,factor):
	b = factor[0]
	c = factor[1]
	d = factor[2]
	e = factor[3]
	y = e*d*2/(np.exp(d*x)+np.exp(-1*d*x))**2 + b*c/(b*x+1)
	return y

if __name__ == '__main__':
	try:
		path = sys.argv[1]
		data = np.load(path)[()]
		kv = data['kinetic_vapor']
		kw = data['kinetic_water']
		vn = data['vapor_number']
		wg = data['water_gyration']
		ct = data['water_center']
	except:
		print("Error occurs.\nUsage: python plot.py <path> <tag> [fit type]('poly'/'lnatan'/'lnsig'/'none')")
		sys.exit(1)

	try:
		fittype = sys.argv[2]
	except:
		fittype = "none"

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
	if fittype != "none":
		#sort_time = np.arange(0,500)    # in ns
		sort_time = np.sort(np.array(vn.keys()[1:]))/1000    # in ns
		sort_time = np.arange(0, sort_time.max()*2000)/1000
		if fittype == "poly":
			popt, pcov = curve_fit(poly, np.array(vn.keys()[1:])/1000, np.array(vn.values()[1:]))
			perr = np.sqrt(np.mean((poly(np.array(vn.keys()[1:])/1000, popt[0], popt[1], popt[2], popt[3]) - np.array(vn.values()[1:]))**2))
		elif fittype == "lnatan":
			popt, pcov = curve_fit(lnatan, np.array(vn.keys()[1:])/1000, np.array(vn.values()[1:]))
			perr = np.sqrt(np.mean((lnatan(np.array(vn.keys()[1:])/1000, popt[0], popt[1], popt[2], popt[3]) - np.array(vn.values()[1:]))**2))
		elif fittype == "lnsig":
			popt, pcov = curve_fit(lnsig, np.array(vn.keys()[1:])/1000, np.array(vn.values()[1:]))
			perr = np.sqrt(np.mean((lnsig(np.array(vn.keys()[1:])/1000, popt[0], popt[1], popt[2], popt[3]) - np.array(vn.values()[1:]))**2))
		else:
			raise RuntimeError("I don't know fit type : %s" % fittype)
		print("Evaporation description : %s" % popt)
		print("Fitting rms error       : %f" % perr)
		if fittype != "poly":
			x = np.arange(1e8,5e10,1e8)
			y = 11000-popt[3]-np.log(popt[0]*x+1)*popt[1]
			zero = np.where(y<0)[0]
			if len(zero) > 0:
				print("evaporation stop at     : %f s" % (float(zero[1]*1e8+1e9)/1e9))
		'''
		for tmp in zeros:
			if np.imag(tmp) == 0 and np.real(tmp)>0:
				sort_time = np.arange(0, np.abs(tmp))
		'''
		if fittype == "poly":
			fit_evap = poly(sort_time, popt[0], popt[1], popt[2], popt[3])
		elif fittype == "lnatan":
			fit_evap = lnatan(sort_time, popt[0], popt[1], popt[2], popt[3])
			fit_evap_01 = lnatan(sort_time, popt[0], popt[1], 0, 0)
			fit_evap_02 = lnatan(sort_time, 0, 0, popt[2], popt[3])
		elif fittype == "lnsig":
			fit_evap = lnsig(sort_time, popt[0], popt[1], popt[2], popt[3])
			fit_evap_01 = lnsig(sort_time, popt[0], popt[1], 0, 0)
			fit_evap_02 = lnsig(sort_time, 0, 0, popt[2], popt[3])
		else:
			raise RuntimeError("I don't know fit type : %s" % fittype)
		ax.plot(sort_time, fit_evap, '--k', linewidth=3.0)
		if fittype in ["lnatan", "lnsig"]:
			ax.plot(sort_time, fit_evap_01, '--b', linewidth=3.0)
			ax.plot(sort_time, fit_evap_02, '--g', linewidth=3.0)
			ax.legend(["Raw Data","Total Fit","loge part","epsilon part"])
		if fittype == "poly":
			ratio = poly_div(sort_time, popt)
		elif fittype == "lnatan":
			ratio = lnatan_div(sort_time, popt)
		elif fittype == "lnsig":
			ratio = lnsig_div(sort_time, popt)
		else:
			raise RuntimeError("I don't know fit type : %s" % fittype)
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
