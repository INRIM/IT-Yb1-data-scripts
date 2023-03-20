# -*- coding: utf-8 -*-
#
# sproc
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



import time
import datetime
import os
import os.path
import shutil
import argparse

import allantools

import re

from numpy import *
from uncertainties import ufloat, ufloat_fromstr, UFloat
from uncertainties.umath import *
from uncertainties import unumpy

from scipy.interpolate import interp1d
from scipy.ndimage import uniform_filter1d, maximum_filter1d
import scipy.stats
from matplotlib.pyplot import *
ion()
close('all')
import matplotlib.offsetbox as offsetbox


import tintervals as ti

import sys

from scipy.constants import k, h
kB = k
Er = 2024.*h

fYb = 518295836590863.61
N0 = 40 #mV ref for density shift ~ 500 atoms




parser = argparse.ArgumentParser(description='Process IT-Yb1 data for interleaved measurements', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('infile', type=argparse.FileType('r'), nargs='+', help='File or list of files to be processed')
parser.add_argument('--dir',  help='Directory for storing results', default='./Proc')
parser.add_argument('-g',  type=str, nargs='+', help='Good data ranges as list of hyphen separated ints. Unbound is empty, e.g. 40- ', default=['40-'])
parser.add_argument('--density','-d', action='store_true', help='Interpret data as density data')
parser.add_argument('--triple', '-t', action='store_true', help='Interpret data as the interleave of 3 cycles')
parser.add_argument('--single', '-s', action='store_true', help='Interpret data as only one cycle')

parser.add_argument('--lattice_l',  type=float, nargs='+', help='lattice depth signal in mV')
parser.add_argument('--lattice_f',  type=float, nargs='+', help='lattice frequency as deviation from 394 798 000 MHz')
parser.add_argument('--Toven',  type=float, nargs='+', help='oven temperature in degree Celsius')
parser.add_argument('--trabi',  type=float, nargs='+', help='rabi time in ms')
parser.add_argument('--Voffsetxp',  type=float, nargs='+', help='Voltage offset along the xp direction.')
parser.add_argument('--noTroom', action='store_true', help='Load room temperature data')
parser.add_argument('--Tdir',  help='Folder of temperature data', default="../../Temperature Data/")

parser.add_argument('--sbfile',  nargs='*', help='sidebands file (used to automatically populate lattice info)', default=['./Vbands/all.dat'])
parser.add_argument('--dfile',  nargs='*', help='density file or files. Automatically strip down to the first 19 chars: num + date + hour.', default=None)



command = sys.argv
args = parser.parse_args()

if not os.path.exists(args.dir):
	os.makedirs(args.dir)

files = [os.path.basename(f.name) for f in args.infile]


def myint(x):
	return int(x) if x else None

good_data = [tuple(map(myint, i.split('-'))) for i in args.g]


if args.triple:
	# number of cycles
	cycles = [1, 2, 3, 4, 5, 6]
	Ncycles = 6

	# name High and Low locks
	H1 = [0,1]
	H2 = [2,3]
	L = [4,5]
	locks = [H1, H2, L]
	names = ["H1", "H2", "L"]

elif args.single:	
	cycles = [1, 2]
	Ncycles = 2
	# name High and Low locks
	L = [0,1]
	locks = [L]
	names = ["L"]
else:
	# number of cycles
	cycles = [1, 2, 3,4]
	Ncycles = 4

	# name High and Low locks
	H = [0,1]
	L = [2,3]
	locks = [H, L]
	names = ["H", "L"]


# prepare to save some infos
conds = [{} for _ in range(len(locks))] # note that [{}]*n create list of n times the same dict!
res = [{} for _ in range(len(locks))] 

# give a proper name
# now of only the first file!
shortname = files[0][:19]
name_with_description = files[0][:-4]
subdir = os.path.join(args.dir,name_with_description)


# create a subfolder
if not os.path.exists(subdir):
	os.makedirs(subdir)

basename = os.path.join(subdir,shortname)
logname = basename + ".txt"
log2name = basename + ".dat"
filo = open(logname, "w")
filo2 = open(log2name, "w")

filo.write("run " + " ".join(command))
filo.write("\n\n")	

filo.write(os.path.relpath(subdir,'../../..'))
filo.write("\n\n")     

filo.write("Data generated with the above command, results in the above folder on: " + str(datetime.datetime.now()) + "\n") 




# load all files in a single array.
arrays = []
aoms = []


for f in args.infile:
	data = genfromtxt(f,skip_footer=1)
	# ignore the first Ncycles lines (usually bugged)
	data = data[Ncycles:]
	# consider only a proper number of lines (multiple of Ncycles)
	x = (len(data)%Ncycles)
	if x != 0:
		data = data[:-x]
	
	arrays += [data]
	
	
	# parse the AOM frequency
	# rewind the file
	f.seek(0)
	# read frist 15 lines (comments) 
	# https://stackoverflow.com/questions/1767513/read-first-n-lines-of-a-file-in-python
	comments = [next(f) for x in range(30)]
	
	
	# regex the aom frequency
	# https://stackoverflow.com/questions/46690385/how-to-read-numerical-data-from-the-comment-of-a-txt-file-into-numpy
	for line in comments:
		aoms += [float(x) for x in re.findall('# Synth [\w\s]+ (\d+[.]\d+)', line)]
		
		# regex some metadata, in the form of list of numbers
		# note only the first file is read 
		if not args.lattice_l:
			temp = re.findall("# lattice_l\s*=\s*([\d.,]+[\d.])", line)
			if temp:
				args.lattice_l = [float(x) for x in temp[0].split(',')]
	
		if not args.lattice_f:
			temp = re.findall("# lattice_f\s*=\s*([\d.,]+[\d.])", line)
			if temp:
				args.lattice_f = [float(x) for x in temp[0].split(',')]
	
		if not args.trabi:
			temp = re.findall("# trabi\s*=\s*([\d.,]+[\d.])", line)
			if temp:
				args.trabi = [float(x) for x in temp[0].split(',')]
	
		if not args.Toven:
			temp = re.findall("# Toven\s*=\s*([\d.,]+[\d.])", line)
			if temp:
				args.Toven = [float(x) for x in temp[0].split(',')]

		if not args.Voffsetxp:
			temp = re.findall("# Voffset_xp\s*=\s*([\d.,]+[\d.])", line)
			if temp:
				args.Voffsetxp = [float(x) for x in temp[0].split(',')]
		
	
	
	# I am done with the files, lets close them
	f.close()

# get the AOM freq
faom = mean(aoms)
if (aoms != faom).any():
	print("WARNING! AOM FREQ NOT CONSISTENT! Sproc is not yet programmed to handle that")


		
# get the the points of file separation
splits = [len(d)/Ncycles for d in arrays]
splits = cumsum(splits)[:-1]

# concatenate all the data in a single array
data = concatenate(arrays, axis=0)		

# get rid of big numbers
epoch0 = data[0,0]



# separate the infos for each cycle in a list
ldata = []
for c in cycles:
	# 1 is the column with the cycle name
	thiscycle = data[data[:,1]==c]	
	ldata += [thiscycle]



# now just take the good data
cdata = []

for d in ldata:
	cdd = []
	for interval in good_data:
		cdd+= [d[slice(*interval)]]
		
	cdd = concatenate(cdd, axis=0)	
	cdata += [cdd]







#========================================
# Plot the data for review of good ranges
#========================================
# freq
figure()
title(shortname)
xlabel('Data point')
ylabel('Frequency /Hz')
# shade good regions
for xmin, xmax in good_data:
	if xmin == None:
		xmin = 0.
	if xmax == None:
		xmax = len(ldata[0])	
	axvspan(xmin,xmax,alpha=0.2,color="yellow")
# mark file separation
for x in splits:
	axvline(x,color='red')
for (d,t) in zip(ldata,cycles):
	plot(d[:,8],label=t)
legend(loc=0)	
pause(0.001)


# natoms
figure()
title(shortname)
xlabel('Data point')
ylabel('Natoms /mV')
# shade good regions
for xmin, xmax in good_data:
	if xmin == None:
		xmin = 0.
	if xmax == None:
		xmax = len(ldata[0])	
	axvspan(xmin,xmax,alpha=0.2,color="yellow")
# mark file separation
for x in splits:
	axvline(x,color='red')
for (d,t) in zip(ldata,cycles):
	plot((d[:,5]+d[:,6])/2,label=str(t))
	#plot(d[:,5],label=str(t) + "L")
	#plot(d[:,6],label=str(t) + "R")
legend(loc=0)
pause(0.001)
	
#excitations
figure()
title(shortname)
xlabel('Data point')
ylabel('Excitations (Maximum of both line sides)')
# shade good regions
for xmin, xmax in good_data:
	if xmin == None:
		xmin = 0.
	if xmax == None:
		xmax = len(ldata[0])	
	axvspan(xmin,xmax,alpha=0.2,color="yellow")
# mark file separation
for x in splits:
	axvline(x,color='red')
for (d,t) in zip(ldata,cycles):
	plot(maximum(d[:,3],d[:,4]), '.', label=str(t))
	#plot(uniform_filter1d(d[:,3],1000), '-', label=str(t) + "L")
	#plot(uniform_filter1d(d[:,4],1000), '-', label=str(t) + "R")
legend(loc=0)
pause(0.001)



# correct freq vs time	
figure()
title(shortname)
xlabel('MJD')
ylabel('Frequency /Hz')
for (d,t) in zip(cdata,cycles):
	plot(ti.mjd_from_epoch(d[:,0]), d[:,8],label=t)
legend(loc=0)
pause(0.001)


#========================================
# Timescale
#========================================
# lets looks at the timetags
tts = [d[:,0] for d in cdata]
# this is the "average" timescale
t = mean(tts,axis=0)
t = t[1:-1]

# hack time -> datapoints
l_tts = [d[:,0] for d in ldata]
l_t = mean(tts,axis=0)
t_datapoints = arange(len(l_t))[isin(l_t,t)]



#========================================
# Analyze each cycle indipendently
#========================================


# calc t0 for each cycle, before removing bad data points
t0 = mean([median(diff(d[:,0])) for d in ldata]) 
rate = 1./t0
#tau = t0 * arange(1,int(len(t)/4))
tau = t0*2**arange(log(int(len(t)/4),2))
taut = t0 * arange(1,int(len(t)/2))

tottime = len(t)*t0


fmt = """
Time infos:
t0 = {:.2f} s
Run duration = {:.2f} s
Meas. time = {:.2f} s
Meas. points = {}"""
out = fmt.format(t0, max(t)-min(t), len(t)*t0, len(t))
filo.write(out+'\n')
print(out)


# calc some useful means

means = [mean(d, axis=0) for d in cdata]
stds = [std(d, axis=0)/sqrt(len(d)-1) for d in cdata]


out = """
Averages per cycle
#Cycle\tExc     \tNum /mV \tExcL    \tExcR    \tNumL /mV \tNumR /mV \tErr /Hz \tGood points"""
filo.write(out+'\n')
print(out)
fmt = "{}" + "\t{:S}"*6 + "\t{:.2uS}\t{:.2f}"


for c,d in zip(cycles,cdata):
	means = mean(d[:,3:8], axis=0)
	stds = std(d[:,3:8], axis=0)/sqrt(len(d)-1)
	mea = unumpy.uarray(means, stds)
	
	out = [c, mean(mea[0:2]), mean(mea[2:4])] + list(mea) + [mean(d[:,9])]
	out = fmt.format(*out)
	filo.write(out+'\n')
	print(out)

 


#========================================
# Now separate H and L loops
#========================================

# now lets do some math
# gets some functions to interpolate the frequencies
ff = [interp1d(d[:,0],d[:,8]) for d in cdata]


# process separately the H and L loops (if interleaving, just one otherwise)
pdata = []
pLOdev = []
pBdev = []
pNdev = []



out = """
Averages per lock
#Lock\tExc     \tNum /mV \tBsplit /Hz\tErr /Hz"""
filo.write(out + '\n')
print (out)
fmt = "{}" + "\t{:S}"*3 + "\t{:.2uS}"



savenum = []

# import tintervals as ti

# figure()
# title(shortname + " - Line Diff ")
# xlabel('MJD')
# ylabel('Rel Frequency Diff.')
# plot(ti.mjd_from_epoch(t),uniform_filter1d((ff[2](t)-ff[0](t))/fYb,20))
# plot(ti.mjd_from_epoch(t),uniform_filter1d((ff[3](t)-ff[1](t))/fYb,20))



for i, lock in enumerate(locks):
	what = names[i]

	# proc base name
	procname = basename + "_proc" + what

	# calculate the frequencies interpoalting
	freqR = ff[lock[0]](t)
	freqL = ff[lock[1]](t)

	freq = (freqR + freqL)/2
	Bsplit = (freqR - freqL)/2

	data = column_stack((t,freqR,freqL,freq,Bsplit))
	pdata += [data]

	# calc allan deviation vs cavity
	tau2, adev, aderr, adn = allantools.oadev(freq/fYb, data_type='freq', rate=rate, taus='all')
	pLOdev += [(tau2, adev, aderr)]	
	

	# calc allan deviation of the Bsplit using allantools!
	tau2, adev, aderr, adn = allantools.oadev(Bsplit, data_type='freq', rate=rate, taus=tau)
	pBdev += [(tau2, adev, aderr)]	
	mBsplit = ufloat(mean(Bsplit), adev[-1])
		

	# calc averages values
	Excs = array([cdata[lock[0]][:,3:5], cdata[lock[1]][:,3:5]])	
	Nums = array([cdata[lock[0]][:,5:7], cdata[lock[1]][:,5:7]])
	Errs = array([cdata[lock[0]][:,7], cdata[lock[1]][:,7]])
	
	
	# new, calc allan variance of number of atoms
	tau2, adev, aderr, adn = allantools.oadev(mean(mean(Nums, axis=0), axis=1), data_type='freq', rate=rate, taus=tau)
	pNdev += [(tau2, adev, aderr)]
	
	
	mExc = ufloat(mean(Excs), std(Excs)/sqrt(Excs.size-1))
	mNum = ufloat(mean(Nums), adev[-1])
	mErr = ufloat(mean(Errs), std(Errs)/sqrt(Errs.size-1))
	

	savenum += [mNum]

	out = fmt.format(what, mExc, mNum, mBsplit, mErr)	
	filo.write(out + '\n')
	print (out)
	

	# save data in my conds list of dict
	conds[i]['Exc'] = mExc
	conds[i]['Num/mV'] = mNum
	conds[i]['Bsplit/Hz'] = mBsplit
	conds[i]['Err/Hz'] = mErr

	# read sidebands data
	if args.lattice_l and args.sbfile:
		# old code to read directly sideband data
		sb = []
		for sbf in args.sbfile:

			# encoding =None avoid a warning
			# deletechars = set() avoid genfromtxt stripping away my "/" in the column names
			try:
				temp = genfromtxt(sbf, names=True, dtype=None, converters={i: ufloat_fromstr for i in arange(4,11)}, encoding=None, deletechars=set())
				sb += [temp]
			except:
				print('WARNING: No sideband data found.')
			
			if len(sb)>1:
				sbdata =concatenate(sb, axis=0)
			else:
				sbdata = array(sb)

		# chose the only tag or the corresponding tag
		lattice_l = args.lattice_l[min(i, len(args.lattice_l)-1)]
		conds[i]['l/mV'] = lattice_l



		

		

		if size(sbdata) > 2:
			# get unique l
			unique_l = unique(sbdata['l/mV'])

			to_extract = ['D/Er', 'Tz/K', 'Tr/K']
			for par in to_extract:
				sbl = sbdata['l/mV']
				var = unumpy.nominal_values(sbdata[par])
				unc = unumpy.std_devs(sbdata[par])


				# before interpolation I need to average different measurement
				# I do not think it is correct to reduce the uncertainty by taking averages, so the uncertainty is also the mean uncertainty!
				# Note that uncertainti in D, Tz and Tr is not really used
				mvar, munc = [], []
				for l in unique_l:
					x = var[where(sbl == l)]
					m = mean(x)
					s = std(x)/sqrt(len(x))*2 # *2 account for low number of points without putting exact t-student values
					

					mvar += [m]
					munc += [sqrt(mean(unc[where(sbl == l)])**2 + s**2)] # final unc is mean of unc summed in quadrature with stddev

				mvar, munc = array(mvar), array(munc)


				# interpolate!
				# if the value is in the table, this is just the average 
				# otherwise give an interpolation of both value and uncertainty
				ivar = interp1d(unique_l, mvar, kind='linear', bounds_error=False, fill_value='extrapolate')
				iunc = interp1d(unique_l, munc, kind='linear', bounds_error=False, fill_value='extrapolate')

				conds[i][par] = ufloat(ivar(lattice_l), iunc(lattice_l))


		# old code to read directly sideband data
#		sb = []
#		for sbf in args.sbfile:
#
#		# encoding =None avoid a warning
#		# deletechars = set() avoid genfromtxt stripping away my "/" in the column names
#		try:
#			temp = genfromtxt(sbf, names=True, dtype=None, converters={i: ufloat_fromstr for i in arange(4,11)}, encoding=None, deletechars=set())
#			sb += [temp]
#		except:
#			print('WARNING: No sideband data found.')
#		
#		if len(sb)>1:
#			sbdata =concatenate(sb, axis=0)
#		else:
#			sbdata = array(sb)
		
		# code suing fit
		# no robust enought to problems on sidebands data

		# # chose the only tag or the corresponding tag
		# lattice_l = args.lattice_l[min(i, len(args.lattice_l)-1)]
		# conds[i]['l/mV'] = lattice_l

		# # new code: read fit result
		# try:
		# 	sbfit = genfromtxt(args.sbfile, names=True, dtype=None, converters={i: ufloat_fromstr for i in arange(0,5)}, encoding=None, deletechars=set())
		# 	filefound = True
		# except:
		# 	print('WARNING: No sideband fit found.')
		# 	filefound = False
		


		
		
		# if filefound and size(sbfit)>0:
		# 	sbfit = atleast_1d(sbfit)		
		# 	D = sbfit['k']*lattice_l
		# 	Tz = sbfit['Az']*D+sbfit['Bz']*D**2	
		# 	Tr = sbfit['Ar']*D+sbfit['Br']*D**2
			
		# 	# convert in K
		# 	Tz = Tz*Er/kB
		# 	Tr = Tr*Er/kB			

		# 	# mean convert to ufloat (from array of ufloat) and is robust  in case more fit results are given
		# 	conds[i]['D/Er'] = D.mean()
		# 	conds[i]['Tz/K'] = Tz.mean()	
		# 	conds[i]['Tr/K'] = Tr.mean()
	
		
			

		
	# other conds
	if args.lattice_f:
		conds[i]['f/MHz'] = args.lattice_f[min(i, len(args.lattice_f)-1)]

	if args.trabi:
		conds[i]['trabi/ms'] = args.trabi[min(i, len(args.trabi)-1)]

	if args.Toven:
		conds[i]['Toven/*C'] = args.Toven[min(i, len(args.Toven)-1)]

	if args.Voffsetxp:
		conds[i]['Voffsetxp/V'] = args.Voffsetxp[min(i, len(args.Toven)-1)]



	conds[i]['faom/Hz'] = faom
	


	
# save data to file
for data,what in zip(pdata, names):
	fmt="%.2f" + "\t%.6f"*4
	header = "Time\tfreqR\tfreqL\tfreq\tBsplit"
	savetxt(basename + "_proc" + what + ".dat", data, fmt=fmt, header = header)


# plot
# mean freq
# figure()
# title(shortname + " - Freq")
# xlabel('Time /s')
# ylabel('Frequency /Hz')
# for data,what in zip(pdata, names):
# 	plot(data[:,0]-epoch0, data[:,3],label=what)
# legend(loc=0)
# pause(0.001)
figure()
title(shortname + " - Clock vs Cavity")
xlabel('Tau /s')
ylabel('Allan dev.')
for data,what in zip(pLOdev, names):
	loglog(data[0], data[1], label=what)

wpm=3e-13
wfm=4e-14
ffm=3e-16
rwm=2e-19 
loglog(tau, (wpm**2*tau**-2 + wfm**2*tau**-1 + ffm**2 + rwm**2*tau)**0.5, label='maser')
grid(which="both")
legend(loc=0)
pause(0.001)




# Bsplit
figure()
title(shortname + " - Bsplit")
xlabel('Time /s')
ylabel('Frequency /Hz')
for data,what in zip(pdata, names):
	plot(data[:,0]-epoch0, data[:,4],label=what)
legend(loc=0)
pause(0.001)

# Bsplit allan
figure()
title(shortname + " - Bsplit")
xlabel('Tau /s')
ylabel('Allan dev. /Hz')
for data,what in zip(pBdev, names):
	loglog(data[0], data[1],label=what)
grid(which="both")
legend(loc=0)
pause(0.001)

# Num allan
figure()
title(shortname + " - Natoms")
xlabel('Tau /s')
ylabel('Allan dev. /mV')
for data,what in zip(pNdev, names):
	loglog(data[0], data[1],label=what)
grid(which="both")
legend(loc=0)
pause(0.001)





#========================================
#get the difference between H and L loops
#========================================


if len(locks) > 1:
	# gets a new functions to interpolate the frequencies
	pff = [interp1d(d[:,0],d[:,3]) for d in pdata]





	for j in arange(len(locks))[:-1]:
		j = -1-j
		for i in arange(len(locks))[:j]:

			# calculate the frequencies interpoalting
			diffe = pff[i](t) - pff[j](t)
			diff_rel = diffe/fYb


			data = column_stack((t,diffe, diff_rel))
			diffname = "{}-{}".format(names[i],names[j])
		
		
			fmt="%.2f\t%.6f\t%e"
			header = "Time\tdiff\trel_diff"
		
			savetxt(basename + "_proc_diff_" + diffname + ".dat" , data, fmt=fmt, header = header)


			figure()
			title(shortname + " - Diff. Shift " + diffname)
			xlabel('Approx data point')
			ylabel('Rel Frequency Diff.')
			plot((t-epoch0)/t0, diff_rel, label=diffname)
			plot((t-epoch0)/t0, uniform_filter1d(diff_rel,1000), label='Rolling mean')
			legend(loc=0)
			pause(0.001)
		
			# calc allan deviation using allantools!	
			tau2, adev, aderr, adn = allantools.oadev(diff_rel, data_type='freq', rate=rate, taus='octave')	
			taut2, tadev, tderr, tdn = allantools.totdev(diff_rel, data_type='freq', rate=rate, taus=taut)
		
			# allantools octave put a point too much
			tau2 = tau2[:-1]
			adev = adev[:-1]
			aderr = aderr[:-1]
			adn = adn[:-1]	


			# calc uncertainty
			# https://tf.boulder.nist.gov/general/pdf/666.pdf
			N = adn
			m = tau2/t0	
			edf = (3*(N-1)/(2*m)-2*(N-2)/N)*4*m**2/(4*m**2+5)
			


			#https://tf.boulder.nist.gov/general/pdf/666.pdf
			admax = adev * (edf/scipy.stats.chi2.ppf(0.1, edf))**.5
			admin = adev * (edf/scipy.stats.chi2.ppf(0.9, edf))**.5


			taufit = append(taut2, tottime)	


			data_out = column_stack((tau2, adev, aderr))
			fmt="%.2f\t%.3e\t%.3e"
			header = "Tau /s\toadev\toadev_unc"
			savetxt(basename + "_proc_diff_" + diffname + "_adev.dat", data_out, fmt=fmt, header = header)
		
			# simple fit
			white = adev*tau2**0.5
			#if len(white) > 10.:			
			#	white = white[where((tau2>10.) & (tau2<tottime/8.))]
			if len(white) > 4.:			
				white = white[where((tau2>10.) & (tau2<tottime/8.))]
			a = mean(white)	

			whiteunc = a*tottime**-0.5

			udiff_rel = ufloat(mean(diff_rel), whiteunc)
			udiffe = udiff_rel*fYb

			figure()
			title(shortname + " - Interl. Stability " + diffname)
			xlabel('Tau /s')
			ylabel('Allan Dev.')
			loglog(taufit, a*taufit**-0.5, "r-", label = f'White = {a:.02}')
			loglog(taut2, tadev, "-", label = "Total")	
			#loglog(tau2, adev, "o", label = "Overl.")
			errorbar(tau2, adev,  yerr = [adev-admin, admax-adev], fmt='o')
			#plot(tau2, adev, fmt='.', ecolor='g')
			grid(which="both")
			legend(loc=0)
			
			text = f'Relative Diff = {udiff_rel:.2uS}'
			ob = offsetbox.AnchoredText(text, loc=3)
			gca().add_artist(ob)

			pause(0.001)



		

			# print the final message
			



			fmt = """
	#Interleaved {} - {}
	Difference = {:.2uS} Hz
	Relative Diff = {:.2uS}

	White noise = {:.02} at 1 s fitted at {:.1f} s.
	"""
		
			out = fmt.format(names[i],names[j],udiffe,udiff_rel, a, tottime)
			print(out)
			filo.write(out + '\n')
		
			res[i]['Diff/Hz'] = udiffe
			res[i]['relDiff'] = udiff_rel
			res[i]['white'] = a
			res[i]['totTime/s'] = tottime
			
			res[-1]['Diff/Hz'] = 0.
			res[-1]['relDiff'] = 0.
			res[-1]['white'] = a
			res[-1]['totTime/s'] = tottime
			
			
			if args.density or abs(savenum[i].n-savenum[j].n) > 15:
				# is density if made in similar conditions
				if(conds[i]['l/mV'] == conds[j]['l/mV'] and conds[i]['f/MHz'] == conds[j]['f/MHz'] and conds[i]['trabi/ms'] == conds[j]['trabi/ms']):
				
					
					fmt = """Density shift = {:.2uS}     (at N0 = {} mV)
					"""
					dcoeff =  udiff_rel/(savenum[i]-savenum[j])*N0
					out = fmt.format(dcoeff, N0)
					print(out)
					filo.write(out + '\n')
					
							
								
					conds[i]['dcoeff'] = dcoeff

					conds[j]['dcoeff'] = dcoeff # TODO: handle multiple cycles
				



# read density info
for i, lock in enumerate(locks):
	# process density tag
	
	if args.density:
		pass
	elif args.dfile:
		dtags = [x[:19] for x in args.dfile]
		dtag = ','.join(dtags)
		
		
		# recover density shift files
		dfiles = []
		for tag in dtags:
			# serach for the tag + ".dat" file
			# and return the first occurrence
			# https://stackoverflow.com/questions/1724693/find-a-file-in-python
			def find(name, path):
				for root, dirs, files in os.walk(path):
					if name in files:
						return os.path.join(root, name)
			
			
			df = find(tag + ".dat", "./Proc")
			if not df:	
				# search for the file i the Data folder (slow!)
				df = find(tag + ".dat", "..")
			
			if df:
				dfiles += [df]
		
		# then load the data
		dd = []
		for df in dfiles:
			# encoding =None avoid a warning
			# deletechars = set() avoid genfromtxt stripping away my "/" in the column names
			try:
				cols = ['lock','l/mV','f/MHz','trabi/ms','dcoeff']
				temp = genfromtxt(df, names=True, dtype=None, usecols=cols, converters={'dcoeff': ufloat_fromstr}, encoding=None, deletechars=set())
				dd += [temp]
			except:
				print('WARNING: No density data found.')
		
		if len(dd)>1:
			ddata =concatenate(dd, axis=0)
		else:
			ddata = array(dd)
		
		
	
		
		# then extract only the relevant data
		a = ddata['lock'] != 'L' # only non "low" cycles (otherwise will double count)
		b = ddata['l/mV'] == conds[i]['l/mV'] # same lattice depth
		c = ddata['f/MHz'] == conds[i]['f/MHz'] # same lattice freq
		d = ddata['trabi/ms'] == conds[i]['trabi/ms'] # same rabi time
		mask = (a & b & c & d)
		dcoeff = ddata[mask]['dcoeff']
					
		if len(dcoeff)>0:
			dcoeff = array(dcoeff)
			weights = array([x.std_dev**-2 for x in dcoeff])
			conds[i]['dcoeff'] = sum(weights*dcoeff)/sum(weights)
		


# look for room temoperature data
if not args.noTroom:
	startdate = datetime.date.fromtimestamp(min(t))
	stopdate = datetime.date.fromtimestamp(max(t))
	step = datetime.timedelta(days=1)
	tfiles = []
	x = startdate
	while x <= stopdate:
		tfiles += [args.Tdir + x.isoformat() + ".dat"]
		x += step
	#print tfiles

	# get my tinetrvals blocks a sstart stop
	#Tvals = ti.array2intervals(t, tgap=30., tblock=30.)
	# for now, rather than importing tintervals here is array2intervals
	def array2intervals(t, tgap=1., tblock=0.):
		""" Calculate from an array of timetags t a 2-d array in the form (start,stop),
	including gaps > tgap and removing intervals < tblock
	"""
		t2 = roll(t,1)
		t3 = roll(t,-1)
		tstarts = t[abs(t-t2) > tgap]
		tstops = t[abs(t-t3) > tgap]

		mask = (tstops - tstarts) >= tblock
		tstarts = tstarts[mask]
		tstops = tstops[mask]

		out = column_stack((tstarts, tstops))
		return out	

	Tvals = array2intervals(t, tgap=90., tblock=90.)
	timesT = []
	minT = []
	maxT = []


	# calculate temperture spread
	for f in tfiles:
		if os.path.isfile(f):
			tdata = genfromtxt(f, skip_footer=1, invalid_raise=False)
			# only loook for temperature in the (start, stop) ranges
			for mint, maxt in Tvals:
				times = tdata[:,0]
		
		
				T = array(tdata[(times > mint) & (times < maxt),1:])
				times2 = array(times[(times > mint) & (times < maxt)])
				
				if times2.size > 0:
					timesT += [times2]
					minT += [T.min(axis=1)]
					maxT += [T.max(axis=1)]
				
				
			


	if len(timesT) >0:
		timesT2 = concatenate(timesT, axis=0)
		minT = concatenate(minT, axis=0)
		maxT = concatenate(maxT, axis=0)

		T = 0.5*(minT+maxT)
		uncT = (maxT - minT)/sqrt(12.)

		
		# take the average temperature
		# note that the mean of the uncertainties is also correct, assuming uniform averaging in time and correlated uncertainties
		mT = mean(T)
		uT = mean(uncT)
		mTroom = ufloat(mT, uT)

		for i, lock in enumerate(locks):
			conds[i]['Troom/*C'] = mTroom
		
		check = amax(t)-amax(timesT2)
		if check > 180.:
			print('WARNING: Temperature data from {} not up to date up to {:.1f} s!'.format(os.path.basename(f), check))	

	else:
		print('WARNING: No corresponding temperature data found.')





# save data to the new .dat file
condkeys = ['Exc', 'Num/mV', 'Bsplit/Hz', 'Err/Hz', 'D/Er', 'Tz/K', 'Tr/K', 'l/mV', 'f/MHz', 'Troom/*C', 'Toven/*C', 'trabi/ms', 'Voffsetxp/V', 'faom/Hz', 'dcoeff']
reskeys = ['Diff/Hz','relDiff','white','totTime/s']


tit = ['key                ', 'lock'] + condkeys + reskeys + ['Files']
#fmt="{}\t{}\t{}"+ "\t{:.2uS}"*(len(condkeys)-3) + "\t{}"*2 + "\t{:.2uS}"*2 + "\t{:.2}"*2 + "\t{}\n"
filo2.write("#" + "\t".join(tit) + "\n")


for i, lock in enumerate(locks):
	# note i wrote zero for the difference of the low lock
	# the get method on the conds dictionary default to None
	out = [shortname, names[i]] + [conds[i].get(x, nan) for x in condkeys] + [res[i].get(x, nan) for x in reskeys] + [",".join(files)]
	
	# if the values in the conds dictionary are sepcified or not change the format
	
	fmt = ["{:.2uS}" if isinstance(x, UFloat) else "{}" for x in out]
	# special format for some thing
	fmt[-3] = "{:.3}"
	fmt[-2] = "{:.1f}"
	
	fmt = "\t".join(fmt) + "\n"
	
	
	textout = fmt.format(*out)
	
	filo2.write(textout)
	
	
	
	
# print all figures		
for i in get_fignums():
	figure(i)
	savefig(basename + "_figure" + str(i) + ".png")
		



filo.close()
filo2.close()


show()
