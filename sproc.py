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
import scipy.stats
from matplotlib.pyplot import *
ion()
close('all')


import sys


fYb = 518295836590865.





parser = argparse.ArgumentParser(description='Process IT-Yb1 data for interleaved measurements')
parser.add_argument('infile', type=argparse.FileType('r'), nargs='+', help='File or list of files to be processed')
parser.add_argument('--dir',  help='Directory for storing results', default='./Proc')
parser.add_argument('-g',  type=str, nargs='+', help='Good data ranges as list of hyphen separated ints. Unbound is empty, e.g. 40- ', default=['40-'])
parser.add_argument('--density','-d', action='store_true', help='Interpret data as density data')
parser.add_argument('--triple', '-t', action='store_true', help='Interpret data as the interleave of 3 cycles')
parser.add_argument('--single', '-s', action='store_true', help='Interpret data as only one cycle')

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


# give a proper name
# now of only the first file!
name = files[0][:19]
name_with_description = files[0][:-4]
subdir = os.path.join(args.dir,name_with_description)


# create a subfolder
if not os.path.exists(subdir):
	os.makedirs(subdir)

basename = os.path.join(subdir,name)
logname = basename + ".txt"
log2name = basename + ".dat"
filo = open(logname, "w")
filo2 = open(log2name, "w")

filo.write(" ".join(command))
filo.write("\n")	

# load all files in a single array.
arrays = []
aoms = []

parse_list = ["lattice_high", "lattice_low", "lattice_freq", "trabi", "Toven"]

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
	comments = [next(f) for x in range(20)]
	
	
	# regex the aom frequency
	# https://stackoverflow.com/questions/46690385/how-to-read-numerical-data-from-the-comment-of-a-txt-file-into-numpy
	for line in comments:
		aoms += [float(x) for x in re.findall('# Synth [\w\s]+ ([0-9]+[.][0-9]+)', line)]
	
	
	
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
title("data vs points")
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



# natoms
figure()
title("data vs points")
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
	
#excitations
figure()
title("data vs points")
xlabel('Data point')
ylabel('Excitations')
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
	plot(d[:,3], label=str(t) + "L")
	plot(d[:,4], label=str(t) + "R")
legend(loc=0)



# correct freq vs time	
figure()
title("data vs time")
xlabel('Time /s')
ylabel('Frequency /Hz')
for (d,t) in zip(cdata,cycles):
	plot(d[:,0]-epoch0, d[:,8],label=t)
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
Meas. points = {}

"""
out = fmt.format(t0, max(t)-min(t), len(t)*t0, len(t))
filo.write(out)
print(out)


# calc some useful means

means = [mean(d, axis=0) for d in cdata]
stds = [std(d, axis=0)/sqrt(len(d)-1) for d in cdata]


out = """Averages per cycle
# Cycle	Exc	Num /mV	ExcL	ExcR	NumL /mV	NumR /mV	Err /Hz	Good points"""
filo.write(out)
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
pBdev = []
pNdev = []



out = """
Averages per lock
# Lock	Exc	Num /mV	Bsplit /Hz	Err /Hz"""
filo.write(out)
print (out)
fmt = "{}" + "\t{:S}"*3 + "\t{:.2uS}\n"



savenum = []





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
	filo.write(out)
	print (out)
	

	# save data in my conds list of dict
	conds[i]['Exc'] = mExc
	conds[i]['Num/mV'] = mNum
	conds[i]['Bsplit/Hz'] = mBsplit
	conds[i]['Err/Hz'] = mErr



	


	
# save data to file
for data,what in zip(pdata, names):
	fmt="%.2f" + "\t%.6f"*4
	header = "Time\tfreqR\tfreqL\tfreq\tBsplit"
	savetxt(basename + "_proc" + what + ".dat", data, fmt=fmt, header = header)


# plot
# mean freq
figure()
title("Freq")
xlabel('Time /s')
ylabel('Frequency /Hz')
for data,what in zip(pdata, names):
	plot(data[:,0]-epoch0, data[:,3],label=what)
legend(loc=0)

# Bsplit
figure()
title("Bsplit")
xlabel('Time /s')
ylabel('Frequency /Hz')
for data,what in zip(pdata, names):
	plot(data[:,0]-epoch0, data[:,4],label=what)
legend(loc=0)

# Bsplit allan
figure()
title("Bsplit")
xlabel('Tau /s')
ylabel('Allan dev. /Hz')
for data,what in zip(pBdev, names):
	loglog(data[0], data[1],label=what)
grid(which="both")
legend(loc=0)

# Num allan
figure()
title("Natoms")
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


	for i in arange(len(locks))[:-1]:

		# calculate the frequencies interpoalting
		diffe = pff[i](t) - pff[-1](t)
		diff_rel = diffe/fYb


		data = column_stack((t,diffe, diff_rel))
		name = "{}-{}".format(names[i],names[-1])
	
	
		fmt="%.2f\t%.6f\t%e"
		header = "Time\tdiff\trel_diff"
	
		savetxt(basename + "_proc_diff_" + name + ".dat" , data, fmt=fmt, header = header)


		figure()
		title("Diff. Shift " + name)
		xlabel('Approx data point')
		ylabel('Rel Frequency Diff.')
		plot((t-epoch0)/t0, diff_rel)
		legend(loc=0)
	
	
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
		savetxt(basename + "_proc_diff_" + name + "_adev.dat", data_out, fmt=fmt, header = header)
	
		# simple fit
		white = adev*tau2**0.5
		#if len(white) > 10.:			
		#	white = white[where((tau2>10.) & (tau2<tottime/8.))]
		if len(white) > 4.:			
			white = white[where((tau2>10.) & (tau2<tottime/8.))]
		a = mean(white)	

		figure()
		title("Interl. Stability " + name)
		xlabel('Tau /s')
		ylabel('Allan Dev.')
		loglog(taufit, a*taufit**-0.5, "r-", label = "White")
		loglog(taut2, tadev, "-", label = "Total")	
		#loglog(tau2, adev, "o", label = "Overl.")
		errorbar(tau2, adev,  yerr = [adev-admin, admax-adev], fmt='o')
		#plot(tau2, adev, fmt='.', ecolor='g')
		grid(which="both")
		legend(loc=0)
	

		whiteunc = a*tottime**-0.5

		udiff_rel = ufloat(mean(diff_rel), whiteunc)
		udiffe = udiff_rel*fYb

	

		# print the final message
		



		fmt = """
#Interleaved {} - {}
Difference = {:.2uS} Hz
Relative Diff = {:.2uS}

White noise = {:.02} at 1 s fitted at {:.1f} s.


"""
	
		out = fmt.format(names[i],names[-1],udiffe,udiff_rel, a, tottime)
		print (out)
		filo.write(out)
	
		
				
			
	
	

#udiffe/(savenum[0]-savenum[1])*1e3


#for i, lock in enumerate(locks):
#	# process density tag
#	
#	if args.dself:
#		density = udiffe/(conds[0]['Num/mV'] - conds[1]['Num/mV'])
#		conds[i]['dcoeff/(Hz/mV)'] = density
#		tags[i]['dtag'] = name
#	elif args.dtag:
#		# chose the only tag or the corresponding tag
#		given = args.dtag[min(i, len(args.dtag)-1)]
#		mores = given.split(',')
#		dtags = [x[:19] for x in mores]
#		dtag = ','.join(dtags)
#		
#		
#		# can handle multiple density files
#		density = []
#		
#		
#		# recover density shift
#		for tag in dtags:
#			# serach for the tag + ".dat" file
#			# and return the first occurrence
#			# https://stackoverflow.com/questions/1724693/find-a-file-in-python
#			def find(name, path):
#				for root, dirs, files in os.walk(path):
#					if name in files:
#						return os.path.join(root, name)
#			
#			
#			df = find(tag + ".dat", ".")
#			if not df:	
#				# search for the file i the Data folder (slow!)
#				df = find(tag + ".dat", "..")
#			
#			if df:
#				#cols to read
#				# encoding =None avoid a warning
#				# deletechars = set() avoid genfromtxt stripping away my "/" in the column names
#				cols = ("Num/mV", "Diff/Hz")
#				dd = genfromtxt(df, names=True, dtype=None, encoding=None, deletechars = set(), usecols=cols, converters={i: ufloat_fromstr for i in cols})
#				
#				diff = dd['Diff/Hz']
#				nums = dd['Num/mV']
#			
#				density += [diff[0]/(nums[0] - nums[1])]
#				
#				 
#				
#				
#		if density:
#			density = array(density)
#			weights = array([x.std_dev**-2 for x in density])
#			conds[i]['dcoeff/(Hz/mV)'] = sum(weights*density)/sum(weights)
#		
#		tags[i]['dtag'] = dtag

#	


## save data to the new .dat file
#condkeys = ['Exc', 'Num/mV', 'Bsplit/Hz', 'Err/Hz', 'U/Er', 'Chi', 'nz', 'fbridge/MHz', 'Troom/*C', 'Toven/*C', 'trabi/ms', 'faom/Hz', 'dcoeff/(Hz/mV)']
#tagkeys = ['sbtag', 'dtag']
#tit = ['key                ', 'lock'] + condkeys + ['Diff/Hz','relDiff','white','toTime/s'] + tagkeys + ['Files']
##fmt="{}\t{}\t{}"+ "\t{:.2uS}"*(len(condkeys)-3) + "\t{}"*2 + "\t{:.2uS}"*2 + "\t{:.2}"*2 + "\t{}\n"
#filo2.write("#" + "\t".join(tit) + "\n")


#for i, lock in enumerate(locks):
#	# note i wrote zero for the difference of the low lock
#	# need to be adapeted if the interleaved locks are not 2
#	# the get method on the conds dictionary default to None
#	out = [name, names[i]] + [conds[i].get(x, None) for x in condkeys] + [udiffe*(i==0), udiff_rel*(i==0), a/sqrt(2), tottime] + [tags[i].get(x, None) for x in tagkeys] + [",".join(files)]
#	
#	# if the values in the conds dictionary are sepcified or not change the format
#	
#	fmt = ["{:.2uS}" if isinstance(x, UFloat) else "{}" for x in out]
#	# special format for some thing
#	fmt[-5] = "{:.3}"
#	fmt[-4] = "{:.1f}"
#	
#	fmt = "\t".join(fmt) + "\n"
#	
#	
#	textout = fmt.format(*out)
#	
#	filo2.write(textout)
	
	
	
	
# print all figures		
for i in get_fignums():
	figure(i)
	savefig(basename + "_figure" + str(i) + ".png")
		



filo.close()
filo2.close()
show()
