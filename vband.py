# -*- coding: utf-8 -*-
#
# vband
# Copyright (C) 2021  Marco Pizzocaro - Istituto Nazionale di Ricerca Metrologica (INRIM)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# SPDX-License-Identifier: MIT


from numpy import *
from scipy.constants import *

from scipy.special import mathieu_b, mathieu_sem, mathieu_cem, eval_genlaguerre, factorial
import scipy.optimize as opt

import time

# not bugged mathieu functions
# to instal gsl
# apt-get install libgsl-dev
# pip3 instal pygsl
# note that special functions are under testing
import pygsl.testing.sf as sf

from matplotlib.pyplot import *
from matplotlib import cm


import os
import fnmatch
import argparse
import sys
from datetime import datetime

from uncertainties import ufloat, ufloat_fromstr, correlated_values
from uncertainties import umath
from uncertainties import unumpy





__version__="1.1"

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Process IT-Yb1 sideband data based on a Born-Oppenheimer model', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('infile', type=argparse.FileType('r'),  nargs='*', help='File or list of files to be processed', default=[])
	parser.add_argument('--dir',  help='Directory for storing results', default='./Vbands')
	parser.add_argument('--log',  help='Filename where to save the results. Results are appended', default='all.dat')
	parser.add_argument('--tag', nargs='+', help='Tag or list of tags for each file', default='-')
	parser.add_argument('--num', type=str, nargs='*', help='Search files from num range instead of full name, e.g. 1-10')

	parser.add_argument('--Dscale',  type=float, help='Scale for the initial guess of D', default=1.3)
	parser.add_argument('--Tzscale',  type=float, help='Scale for the initial guess of Tz', default=0.2)
	parser.add_argument('--Trscale',  type=float, help='Scale for the initial guess of Tr', default=0.2)
	parser.add_argument('--Tz',  type=float, help='Fix the initial guess of Tz, then Tzscale is ignored', default=None)
	parser.add_argument('--Tr',  type=float, help='Fix the initial guess of Tr, then Trscale is ignored', default=None)
	parser.add_argument('--fac',  type=float, help='Quality of integration parameter', default=10)

	parser.add_argument('--nofit', action='store_true', help='Do not fit for all processed files')
	parser.add_argument('--Tzquad', action='store_true', help='Fit a quadratic Tz vs D')
	parser.add_argument('--Trlin', action='store_true', help='Fit a linear Tr vs D')
	#parser.add_argument('--Nfit', action='store_true', help='Show N vs U for processed files')



	command = sys.argv
	args = parser.parse_args()


	if not os.path.exists(args.dir):
		os.makedirs(args.dir)


	ioff()
	close('all')


m = 2.83846417e-25	# atomic mass
fL = 394798e9		# lattice frequency
lL = c/fL		# lattice wavelength
Er = h**2/(2*m*lL**2)	# recoil energy
vrec = h/(2*m*lL**2)	# recoil frequency
w = 45e-6		# lattice waist # note final result does not depend on w

# store boltzmann for later
kB = k

# problem scale
k = 2*pi/lL
kappa = sqrt(2)/w 

# clock k
wc = 518295836590863.
lc = c/wc
kc = 2*pi/lc



## as reference: Frequencies (Blatt2009 eq. 4 and 5)
#vz = 2*vrec*sqrt(D)
#vr = vz*kappa/k


def U(rho, D, nz):
	"""Return the lattice potential surface as Beloy2020 Eq. D2.
	
	Inputs
	------
	rho : radial coordinates in units of kappa**-1
	D : depth of the lattice in Er
	nz : longitudinal quantum number
	
	Returns
	-------
	U : lattice potential surface in Er
	"""
	Drho = D*exp(-rho**2)
	return (sf.mathieu_b(nz+1,Drho/4) -Drho/2)

def R(E, D, nz):
	"""Inverse uf U (Beloy2020 pag. 4)
	
	Inputs
	------
	E : potential energy in Er
	D : depth of the lattice in Er
	nz : longitudinal quantum number
	
	Returns
	-------
	r : radial coordinates in units of kappa**-1
	"""
	mineps = 0.
	maxeps = 2.7  #max radius for nz = 0 and D=1500
	try:

		res = opt.root_scalar(lambda x, *args: U(x, *args) - E, args=(D, nz), bracket=[mineps,maxeps], method='brentq')
	except ValueError:
		return 0.
	
	return res.root

R = vectorize(R)


def DeltaU(rho, D, nz, dn=1):
	"""Return the difference between lattice potential surfaces.
	
	Inputs
	------
	rho : radial coordinates in units of kappa**-1
	D : depth of the lattice in Er
	nz : longitudinal quantum number of the starting level
	dn : longitudinal quantum number jump (default: 1)
	
	Returns
	-------
	DeltaU : difference U(rho, nz+dn) - U(rho, nz)
	"""

	return U(rho, D, nz+dn) - U(rho, D, nz)


# Harmonic Oscillator
def Omega(rho, D, nz, dn=1):
	"""Normalized Rabi frequency for the harmonic oscillator (Wineland1979 eq. 31)
	
	Inputs
	------
	rho : radial coordinates in units of kappa**-1
	D : depth of the lattice in Er
	nz : longitudinal quantum number of the starting level
	dn : longitudinal quantum number jump (default: 1)
	
	Returns
	-------
	Omega : Normalized Rabi frequency  (between 0 and 1)
	"""

	Drho = D*exp(-rho**2)
	eta = kc/k/sqrt(2)/Drho**0.25 

	
	return exp(-eta**2/2)*sqrt(factorial(nz)/factorial(nz+dn))*eta**(dn)*eval_genlaguerre(nz,dn,eta**2)
	
	
	
	

# Using Mathieu Functions
def OmegaMat(rho, D, nz, dn=1):
	"""Normalized Rabi frequency for Mathieu wavefunctions (Beloy2020 appendix)
	
	Inputs
	------
	rho : radial coordinates in units of kappa**-1
	D : depth of the lattice in Er
	nz : longitudinal quantum number of the starting level
	dn : longitudinal quantum number jump (default: 1)
	
	Returns
	-------
	Omega : Normalized Rabi frequency (between 0 and 1)
	"""

	Drho = D*exp(-rho**2)
	lim = pi/(2*k)
	if dn % 2:
		res2 = integ.quad(lambda z: 2*k/pi * sf.mathieu_se(nz+1,  Drho/4, (k*z + pi/2)) * sin(kc*z) * sf.mathieu_se(nz+1+dn, Drho/4, (k*z + pi/2)), 0,lim) # pygsl -- slower but no bugs!
	else:
		res2 = integ.quad(lambda z: 2*k/pi * sf.mathieu_se(nz+1,  Drho/4, (k*z + pi/2)) * cos(kc*z) * sf.mathieu_se(nz+1+dn, Drho/4, (k*z + pi/2)), 0,lim) # pygsl -- slower but no bugs!
	
	# integral is even
	return 2*abs(res2[0])
OmegaMat = vectorize(OmegaMat)


def Gn(E, D, nz):
	"""Density of states for the lattic trap (Beloy2020 eq. 11)
	
	Inputs
	------
	E : potential energy in Er
	D : depth of the lattice in Er
	nz : longitudinal quantum number of the starting level
		
	Returns
	-------
	G : density of states at energy E (same units of Beloy2020 fig. 3)
	
	"""

	return m/(2*hbar**2)*R(E, D, nz)**2 *Er /k**2 # *Er convert to Er units


def Gr(rc, D, nz):
	"""Density of states for the lattic trap given rc(Beloy2020 eq. 11)
	
	Inputs
	------
	rc : Condon point at energy E
	D : depth of the lattice in Er
	nz : longitudinal quantum number of the starting level
		
	Returns
	-------
	G : density of states at energy E (same units of Beloy2020 fig. 3)
	
	"""

	return m/(2*hbar**2)*rc**2 *Er /k**2 # *Er convert to Er units

def maxnz(D):
	"""Return the maximum nz for a given depth"""
	# ansatz 	twice the harmonic oscillators levels
	maxn = int(-U(0,D,0)/sqrt(D))
	test = arange(maxn)
	return amax(where(U(0,D,test)<0))
	


#lorentzian
def lor(x, x0, w):
	""" Simple lorentzian with center x0 and HWHM w peak 0.5
	"""
	den = 1 + 1/w**2*(x - x0)**2
	return 0.5/den

# both sideband
 #it is faster to calculate sidebands at the same time
def sbands(x, D, Tz, Tr, b, r, wc, dn=1, Emax=0.,  fac=10):
	"""Blue sideband as a sum of lorentzian.
	
	Inputs
	------
	x : frequency in Hz
	D : depth of the lattice in Er
	Tz : longitudinal temperature in K
	Tr : radial temperature in K
	b : blue sidebands scaling
	r : red sidebands scaling 
	wc : carrier half width half maximum
	dn : order of the sideband (default: 1)
	Emax : max energy levels populated (default: 0)
	fac : control the number of lorentzian to be used in the sum
	higher number will give smoother sidebands but take more computational time
	(default: 10)
		
	Returns
	-------
	Both sidebands as excitation.
	
	"""
	Nz = int(maxnz(D)*1.+0.5)
	betar = Er/(kB*Tr)
	betaz = Er/(kB*Tz)

	# simple exp factor for a given Tz
	# will be just used computationally to reduce the number of lorentzian used at high nz
	# r = exp(-betaz)

	tot = zeros(x.shape)
	tnorm = 0

	
	for nz in arange(Nz+1):
		Emin = U(0, D, nz)
		#Emax = 0 # I do not think this should change -- atoms close to this are lost, but I should kept the full normalization
	
		# this just save computational time
		# use less samples for high levels
		# formula with r is normalized exponential distribution (from gemotric series)
		# high dn use higher number because it has sharper lorentzians
		#N = int(Natoms*r**nz/(1-r**Nz)*(1-r)+2.)*dn 
		
		# new way to estabish the number of lorentzians
		N = int(DeltaU(0, D, nz, dn)*fac*(nz+1)**-0.5) # other idea use Omega
	

		
		# Uniform sampling in E
		EE = linspace(Emin, Emax, N)[:, newaxis]
		rc = R(EE, D, nz)
		#dE = (Emax - Emin)/N		
		
		
			
				
		# calc normalization		
		pp = Gr(rc, D, nz) * exp(-(EE-Emin)*betar) * exp(-Emin*betaz) 
		tnorm += trapz(pp, EE, axis=0) #sum(pp, axis=0) *dE #trapz(pp, EE, axis=0) #trapz is a bit slower, but handles better different Ns
		
		# blue		
		x0 = DeltaU(rc, D, nz, dn)*vrec
		ff = Omega(rc, D, nz, dn)*wc
		
		# sum lorentzian for blue sideband - note cutoff on energy
		blue = pp*lor(x, x0, ff) *(U(rc, D, nz+dn)<Emax)

		res = b*trapz(blue, EE, axis=0) #sum(yy, axis=0)*dE #trapz(yy, EE, axis=0) # speed sum > trapz > simps
		
		tot += res

		
		
		# red
		if nz >= dn:		
			# rc = R(EE, D, nz)  # same as blue
			x0 = DeltaU(rc, D, nz, -dn)*vrec
			ff = Omega(rc, D, nz-dn, dn)*wc
			
			# sum lorentzian on red sideband					
			red = pp*lor(x, x0, ff)
			
			res = r*trapz(red, EE, axis=0) #trapz(yy, EE, axis=0) # sum(yy, axis=0)*median(diff(EE, axis=0)) #
			
			tot += res
		
		
	return tot/tnorm


def distrib(E, D, Tz, Tr, Emax=0.):
	""" Distribution of atoms in the trap.	
	"""
	Nz = maxnz(D)
	betar = Er/(kB*Tr)
	betaz = Er/(kB*Tz)
	nz = arange(Nz+1)[:, newaxis]

	Emin = U(0, D, nz)
	pp = Gn(E, D, nz) * exp(-(E-Emin)*betar) * exp(-Emin*betaz) *(E < Emax)
	
	
	AEmin = amin(Emin)
	
	EE = arange(AEmin, Emax+1, 1.)
	ppNorm = Gn(EE, D, nz) * exp(-(EE-Emin)*betar) * exp(-Emin*betaz) 
	znorm = trapz(ppNorm, EE, axis=1)
	norm = sum(znorm)
	
		
	return pp/norm, znorm/norm







def fit_lor(x, A, x0, w):
	return A * lor(x, x0, w)


def fit_sbands(x, A, D, Tz, Tr, Sb=0.): # note misisng w from parameters
	return  sbands(x, D, Tz, Tr, A, A*(1+Sb), w0, dn=1, Emax=0,  fac=args.fac)


if __name__ == '__main__':
	allfigname = os.path.join(args.dir,os.path.splitext(args.log)[0] + '.png')
	alltxtname = os.path.join(args.dir,args.log)
	fitbasename = os.path.splitext(args.log)[0]


	alltxt = open(alltxtname, 'a')
	labels = ['Ac','fc/Hz','wc/Hz', 'A', 'D/Er','Tz/K','Tr/K', 'Nat/mV']
	tit = "#file" + " "*32 + "\ttag\tl/mV\tspin\t" + '\t'.join(labels) + '\n'
	logfmt="{}\t{}\t{}\t{}"+  "\t{:.2uS}"*len(labels) + "\n"
	alltxt.write(tit)


	if args.num:
		def myint(x):
			return int(x) if x else None
	
		lim_list = [tuple(map(myint, i.split('-'))) for i in args.num]
		corrected_lim_list = [(x[0],x[0]+1) if len(x)==1  else (x[0],x[1]+1) for x in lim_list]
		num_list = concatenate([arange(*x) for x in corrected_lim_list])
		tag_list = ["{:03d}".format(x) for x in num_list]
		
		for tag in tag_list:
			# serach for the tag  file
			# and return the first occurrence
			# https://stackoverflow.com/questions/1724693/find-a-file-in-python
			def find(pattern, path):
				result = []
				for root, dirs, files in os.walk(path):
					for name in files:
						if fnmatch.fnmatch(name, pattern):
							result.append(os.path.join(root, name))
				return result
			
			
			df = find(tag + "_*.dat", ".")
			if not df:	
				print("WARNING: I have not found a file starting with {}".format(tag))
			else:
				args.infile += [open(f,'r') for f in df]

	Nf = len(args.infile)



	for i, fil in enumerate(args.infile):
		basename = os.path.basename(fil.name)[:-4]
		numname = os.path.basename(fil.name)[:3]
		figname = os.path.join(args.dir,'fig_' + basename+'.png')
		txtname = os.path.join(args.dir,basename+'.txt')
		
		
		# parse metadata
			
		# rewind the file
		fil.seek(0)
		# read first 15 lines (mostly comments) 
		# https://stackoverflow.com/questions/1767513/read-first-n-lines-of-a-file-in-python
		comments = [next(fil) for x in range(15)]
		
		
		# regex l=... for tag and s=... for spin
		# https://stackoverflow.com/questions/46690385/how-to-read-numerical-data-from-the-comment-of-a-txt-file-into-numpy
		# regex inspired from https://stackoverflow.com/questions/63544447/python-regex-parse-string-and-extract-key-value-pairs
		
		ll = []
		ss = []
		
		for line in comments:
			ll += re.findall(r'l\s?=\s?(\".*?\"|\S*)', line) 
			ss += re.findall(r's\s?=\s?(\".*?\"|\S*)', line) 
			
		l = ','.join(ll)
		spin = ','.join(ss)
		
		tag = args.tag[min(i, len(args.tag)-1)]
		
		
		
		
		
		
		print('file: '+fil.name)
		print('tag: '+tag)
		print('l: '+l)
		print('spin: '+spin)
		
		# fitting procedure
		# =================
		
		data = genfromtxt(fil, usecols=(2,3,4,5,6,7),  skip_footer=1, skip_header=1)
		
		
		freq = data[:, 0]
		freq = around(freq, 4)
		funique = unique(freq) #find scan frequency values

		# average data
		ave = []
		for f in funique:
			ave += [mean(data[where(freq==f)], axis=0)]

		
		ave = array(ave)	
		#freq = funique
		
		
		# unpack!
		freq, g,b,e,num,exc = ave.T
		G = g-b
		E = e-b
		
		
		
		# lets get down some scales from labedit exc
		medex = mean(exc)
		test = where(exc > medex)[0]
		
		fmi = freq[test[0]]
		fma = freq[test[-1]]
		fce = (fmi+fma)/2
		fxc = (fma-fmi)/2
		
		
		
		# now recalculate excitation eta on the red side of the spectrum
		nscale = 5e3 # Hz
		mask = funique < fce + nscale
		
		def findeta(eta, g, e):
			return std(g+e/eta)
			
		res = opt.minimize(findeta, 0.97, args=(G[mask],E[mask]))
		eta = res.x[0]	# we found abug for wich the first blue pulse gives less signal than the third, so eta can be > 1 !?!
		
			
		N = G + E/eta
		baseN = median(N[mask])
		
		
		
		
		
		
		# now use the new data
		# this substantially fit just th excitation
		exc = E/baseN


		# fit carrier
		cscale = 10e3 # Hz
		mask = abs(funique - fce) < cscale
		
		cfreq = freq[mask]
		cexc = exc[mask]
		
		a = amax(cexc)*2
		fc = cfreq[argmax(cexc)]
		
		copt, ccov = opt.curve_fit(fit_lor, cfreq, cexc, p0=[a, fc+100, cscale/10])

		# save carrier width and translate to center
		w0  = copt[2]
		freq -= copt[1]
		copt[1] = 0
		
		# remove carrier for future fits
		sexc = exc - fit_lor(freq, *copt)
		
		# guess sbands parameters
		A = sum(cexc)*median(diff(freq))/w0
		D = fxc**2/(4*vrec**2)*args.Dscale
		
		if args.Tz:
			Tz =  args.Tz
		else:
			Tz =  D*Er/kB*args.Tzscale
				
		if args.Tr:		
			Tr = args.Tr
		else:
			Tr = D*Er/kB*args.Trscale
		

		# fit
		sopt, scov = opt.curve_fit(fit_sbands, freq, sexc, p0=[A,D,Tz,Tr], bounds=([0., 0., 1e-7, 1e-7],[300, 2*D, 1e-4, 1e-4]))
		

		
		
		
		# Plot fit
		ff = linspace(amin(freq),amax(freq),400)
		fit = fit_sbands(ff, *sopt) + fit_lor(ff, *copt)
		
		
		# single plot
		figure('comp', figsize=(12.8,4.8))
		plot(freq,exc, 'o', label=numname, color=cm.viridis(i/Nf))
		plot(ff,fit, color=cm.viridis(i/Nf))
		grid(True)
		ylabel('Excitation probability')
		xlabel('Frequency detuning /Hz')
		legend(loc=0)
		
		
		
		# for each file
		figure(figsize=(12.8,4.8))
		title(basename)
		plot(freq,exc, 'o')
		plot(ff, fit)
		grid(True)
		ylabel('Excitation probability')
		xlabel('Frequency detuning /Hz')
		savefig(figname)
		#pause(0.001) # matplotlib need a pause to catch up -- see https://stackoverflow.com/questions/28269157/plotting-in-a-non-blocking-way-with-matplotlib
		# for many figure do not show
		if Nf>9:
			close()
			
		
		
		# extract parameters
		c_par = correlated_values(copt, ccov)
		s_par = correlated_values(sopt, scov)

		Natoms = ufloat(mean(data[:,-2]) ,std(data[:,-2]))

		print('A = {:.2uS}'.format(s_par[0]))
		print('D = {:.2uS} Er'.format(s_par[1]))
		print('vz = {:.2uS} kHz'.format(s_par[1]**0.5*vrec*2*1e-3))
		print('Tz = {:.2uS} K'.format(s_par[2]))
		print('Tr = {:.2uS} K'.format(s_par[3]))
		print('wc = {:.2uS} Hz'.format(c_par[2]))
		print('Nat = {:.2uS} mV'.format(Natoms))
		print('eta = {:.3}'.format(eta))
		print('\n')




			
			
		# save data in tabular form on the common file
		alltxt.write(logfmt.format(basename, tag, l, spin, *c_par, *s_par, Natoms))

		
	if Nf>1:
		figure('comp')
		savefig(allfigname)

	if Nf>0:
		alltxt.write("# Above data generated on: " + str(datetime.now()) + "\n")		
		alltxt.write("# With the command line: " + " ".join(command) + '\n')		
		#alltxt.write("# Script version: " +  __version__ + '\n')		

	alltxt.close()

	

	if not args.nofit:
		# deletechars = set() avoid genfromtxt stripping away my "/" in the column names
		data = genfromtxt(alltxtname, names=True, dtype=None, converters={i: ufloat_fromstr for i in arange(4,4+len(labels))}, encoding=None, deletechars=set())
		
		D = unumpy.nominal_values(data['D/Er'])
		uD = unumpy.std_devs(data['D/Er'])
		
			
		Tz = unumpy.nominal_values(data['Tz/K'])
		uTz = unumpy.std_devs(data['Tz/K'])


		Tr = unumpy.nominal_values(data['Tr/K'])
		uTr = unumpy.std_devs(data['Tr/K'])
			
		
		Na = unumpy.nominal_values(data['Nat/mV'])
		uNa = unumpy.std_devs(data['Nat/mV'])
		
		
		l = data['l/mV'].astype(float)
		
		
		
		# quadratic fit
		def fun(x, a, b):
			return a*(x*Er)/kB + b*x**2*Er/kB
		
		eps = 1e-6
		
		if args.Tzquad:
			zopt, zcov = opt.curve_fit(fun, D, Tz, sigma=uTz, p0=[0.2, 0.001])
		else:	
			zopt, zcov = opt.curve_fit(fun, D, Tz, sigma=uTz, p0=[0.2, 0.], bounds=([0,-eps],[10., eps]))

		if not args.Trlin:
			ropt, rcov = opt.curve_fit(fun, D, Tr, sigma=uTr, p0=[0.2, 0.001])
		else:
			ropt, rcov = opt.curve_fit(fun, D, Tr, sigma=uTr, p0=[0.2, 0.], bounds=([0,-eps],[10., eps]))	

		z_par = correlated_values(zopt, zcov)
		r_par = correlated_values(ropt, rcov)
		
		print('Tfit\n')
		print('Az = {:.2uS}'.format(z_par[0]))
		if args.Tzquad:
			print('Bz = {:.2uS}'.format(z_par[1]))
		print('Ar = {:.2uS}'.format(r_par[0]))
		if not args.Trlin:
			print('Br = {:.2uS}'.format(r_par[1]))
		
		
		
		
		
		
		
			
		
		figure()
		scale=1e6

		for i, x in enumerate(set(data['spin'])):
			mask = where(data['spin']==x)
			errorbar(D[mask], Tz[mask]*scale, yerr=uTz[mask]*scale, xerr=uD[mask], fmt='o', color=cm.tab20(1-i), label='Tz, s={}'.format(x))
			errorbar(D[mask], Tr[mask]*scale, yerr=uTr[mask]*scale, xerr=uD[mask], fmt='d', color=cm.tab20(3-i), label='Tr, s={}'.format(x))

		uu = linspace(0, max(D)*1.1, 100)

		plot(uu, fun(uu, *zopt)*scale, color=cm.tab20(4))
		plot(uu, fun(uu, *ropt)*scale,  color=cm.tab20(6))
		
		legend(loc=0)
		xlabel('D/ Er')
		ylabel('T /µK')
		grid(True)

		left, right = xlim()
		if right < 550:
			xlim(left, 550)


		tfitname = os.path.join(args.dir,fitbasename + '-tfit.png')
		savefig(tfitname)


		figure()
		scale=1
		
		for i, x in enumerate(set(data['spin'])):
			mask = where(data['spin']==x)
			errorbar(D[mask], Na[mask]*scale, yerr=uNa[mask]*scale, xerr=uD[mask], color=cm.tab20(1-i), fmt='o', label='Natoms s={}'.format(x))
		
		legend(loc=0)
		xlabel('D/ Er')
		ylabel('Natoms /mV')
		grid(True)
		xlim(0,None)
		
		left, right = xlim()
		if right < 550:
			xlim(left, 550)

		
		nfitname = os.path.join(args.dir,fitbasename + '-N.png')
		savefig(nfitname)
		
		
		# also fit tags
		
		def fun(x, A):
			return A*x
			
		lopt, lcov = opt.curve_fit(fun, l, D, sigma=uD)
		
		
		figure()
		

		errorbar(l, D, yerr=uD, fmt='o')

		ll = linspace(0, max(l)*1.1, 100)

		plot(ll, fun(ll, *lopt))
		
		xlabel('l /mV')
		ylabel('D/ Er')
		
		grid(True)
		
		Dfitname = os.path.join(args.dir,fitbasename + '-Dfit.png')
		savefig(Dfitname)
		
		
		l_par = correlated_values(lopt, lcov)
		
		print('Dfit\n')
		print('k = {:.2uS}'.format(l_par[0]))
					
		Ttxtname = os.path.join(args.dir,fitbasename + '-fit.dat')
		Ttxt = open(Ttxtname, 'w')
		Tlabels = ['Az', 'Bz', 'Ar', 'Br', 'k']
		Ttit = '#' + '\t'.join(Tlabels) + '\n'
		Tlogfmt="\t".join(["{:.2uS}"]*len(Tlabels)) + "\n"
		
		Ttxt.write(Ttit)
		Ttxt.write('#\n# Temperature fits\n')
		Ttxt.write('# kB*T/Er = A*(D/Er) + B*(D/Er)**2\n')
		Ttxt.write('#\n# Depth fit\n')
		Ttxt.write('# D/Er = k*l/mV\n#\n')
		Ttxt.write(Tlogfmt.format(*z_par,*r_par,*l_par))
		Ttxt.close()
	

		
		
		
	show(block=False) # note i prefer to be ioff for this script
		

