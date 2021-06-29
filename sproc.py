
# Sproc3
# process Yb data for interleaved measurements

# changelog
# 20210115 - swapped from pdf figures to png

# based on sproc.py, lates version found on Data/20181017
# changes from sproc
# *now its is based on argparse 
# *lives in the daily folder (and not in the ana folder)
# *while analyzing multiple files only keep the first as name
# -can point to other files to automatically populate density and sidebands data

__version__="1.0"

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

from configparser import SafeConfigParser

import sys

ion()
close('all')
fYb = 518295836590865.





parser = argparse.ArgumentParser(description='Process the Yb data for interleaved measurements')
parser.add_argument('infile', type=argparse.FileType('r'), nargs='+', help='File or list of files to be processed')
parser.add_argument('--dir',  help='Directory for storing results', default='./ana3')
parser.add_argument('-g',  type=str, nargs='+', help='good data ranges as list of hyphen separated ints. Unbound is empty, e.g. 40- ', default=['40-'])
parser.add_argument('-Toven',  type=ufloat_fromstr, help='explicit Toven with uncertanty in degree Celsius')
parser.add_argument('-trabi',  type=float, help='explicit rabi time in ms', default = 80.)
parser.add_argument('-sbfile',  nargs='*', help='associated processed sidebands file or files', default=['./sidebands3/last.dat'])
parser.add_argument('-sbtag', nargs='*', help='tag or tags of the sidebands data. If more than one, assumed it is one per lock.', default=None)
parser.add_argument('-dtag', nargs='*', help='tag or tags of the density data. If more than one, assumed is one per lock. Automatically strip down to the first 19 chars: num + date + hour. ', default=None)
parser.add_argument('-dself', action='store_true', help='interpret analyzed data as density data. overwrite dtag')


command = sys.argv
args = parser.parse_args()


if not os.path.exists(args.dir):
	os.makedirs(args.dir)




files = [os.path.basename(f.name) for f in args.infile]


def myint(x):
	return int(x) if x else None



good_data = [tuple(map(myint, i.split('-'))) for i in args.g]





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
tags = [{} for _ in range(len(locks))]


# VS PTB here


#files = [ "010_20160510_125804_vsptb.dat","012_20160510_174326_vsptb.dat"]
#good_data = [(40,2091),(2319,16967),(17811,22325),(22438,23325),(23468,None)]
## number of cycles
#cycles = [1, 2]
#Ncycles = 2
## name High and Low locks
#H = [0,1]
#locks = [H]
#names = ["U"]

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
filo.write("\nSproc3 version " + str(__version__) + " processing\n")
filo.write("\n".join(files) + "\n")



print ("Processing")
print ("\n".join(files) + "\n")


# load all files in a single array.
arrays = []
aoms = []
for f in args.infile:
	# path = os.path.join(".",f)
	# copy the file in the current dir
	# shutil.copyfile(path,f)
	
	
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
	comments = [next(f) for x in range(15)]
	
	
	# regex the frequency
	# https://stackoverflow.com/questions/46690385/how-to-read-numerical-data-from-the-comment-of-a-txt-file-into-numpy
	
	for line in comments:
		match1 = re.search('# Synth  with center freq [0-9]+[.][0-9]+', line)
		if match1:
			# save the freq already as float float
			aoms += [float(match1.group(0)[26:])]
	
	
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
#data[:,0] = data[:,0] - data[0,0]
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
#savetxt("witchhunting.dat", column_stack((ldata[0][:,0], ldata[0][:,8])))


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
	plot(d[:,5],label=str(t) + "L")
	plot(d[:,6],label=str(t) + "R")
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

#========================================
# Timescale
#========================================
# lets looks at the timetags
tts = [d[:,0] for d in cdata]
# this is the "average" timescale
t = mean(tts,axis=0)
t = t[1:-1]



# I removed the tmeperature calculations here because pcloud do not sync the temp data very quickly
# and it is kind of annoying to handle all possible exceptions of missing data
# the script to calculate the uncertainty should in any case be able to recover temp data on his own, withouth sproc3 intervention

## from tintervals.py
#def array2intervals(t, tgap=1., tblock=0.):
#	""" Calculate from an array of timetags t a 2-d array in the form (start,stop),
#including gaps > tgap and removing intervals < tblock
#"""
#	t2 = roll(t,1)
#	t3 = roll(t,-1)
#	tstarts = t[abs(t-t2) > tgap]
#	tstops = t[abs(t-t3) > tgap]

#	mask = (tstops - tstarts) >= tblock
#	tstarts = tstarts[mask]
#	tstops = tstops[mask]

#	out = column_stack((tstarts, tstops))
#	return out

## intervals as (start, stop)	
#tintervals = array2intervals(t, tgap=30., tblock=30.)	


## get the average room temperature from Tmonitor folder

## get the date and the corresponding temp files
#startdate = datetime.date.fromtimestamp(min(t))
#stopdate = datetime.date.fromtimestamp(max(t))
#step = datetime.timedelta(days=1)
#tfiles = []
#x = startdate
#while x <= stopdate:
#	tfiles += ["../../../Yblab/Tmonitor/" + x.isoformat() + ".txt"]
#	x += step
#print tfiles



#timesT = []
#minT = []
#maxT = []


#got_some_data = False

## calculate temperture spread
#for f in tfiles:
#	if os.path.isfile(f):
#		# only loook for temoperature sin the (start, stop) ranges
#		for mint, maxt in tintervals:
#			tdata = genfromtxt(f, skip_footer=1)
#			times = tdata[:,0]
#		
#		
#			T = array(tdata[(times > mint) & (times < maxt),1:])
#			times2 = array(times[(times > mint) & (times < maxt)])
#		
#			timesT += [times2]
#			minT += [T.min(axis=1)]
#			maxT += [T.max(axis=1)]
#			
#			if times2:
#				gto_some_data=True


## open the figure here to save fig numbering even if missing T data
#figure()





#if got_some_data:
#	timesT2 = concatenate(timesT, axis=0)
#	minT = concatenate(minT, axis=0)
#	maxT = concatenate(maxT, axis=0)
#	



#	T = 0.5*(minT+maxT)
#	uncT = (maxT - minT)/sqrt(12.)

#	
#	# take the average temperature
#	# note that the mean of the uncertainties is also correct, assuming uniform averaging in time and correlated uncertainties
#	mT = mean(T)
#	uT = mean(uncT)
#	mTroom = ufloat(mT, uT)

#	for i, lock in enumerate(locks):
#		conds[i]['Troom/*C'] = mTroom


#	# plot	
#	#figure()
#	title("Temp vs time")
#	xlabel('Time /s')
#	ylabel('Temperature /degree C')
#	plot(timesT2 - epoch0, maxT, label = 'Max T')
#	plot(timesT2 - epoch0, minT, label = 'Min T')

#	axhspan(mT - uT,mT+uT,alpha=0.2,color="red")
#	axhline(mT,color='red')

#	legend(loc=0)




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
print (out)


# calc some useful means

means = [mean(d, axis=0) for d in cdata]
stds = [std(d, axis=0)/sqrt(len(d)-1) for d in cdata]


out = """Averages per cycle
# Cycle	Exc	Num /mV	ExcL	ExcR	NumL /mV	NumR /mV	Err /Hz	Good points
"""
filo.write(out)
fmt = "{}" + "\t{:S}"*6 + "\t{:.2uS}\t{:.2f}\n"


for c,d in zip(cycles,cdata):
	means = mean(d[:,3:8], axis=0)
	stds = std(d[:,3:8], axis=0)/sqrt(len(d)-1)
	mea = unumpy.uarray(means, stds)
	
	out = [c, mean(mea[0:2]), mean(mea[2:4])] + list(mea) + [mean(d[:,9])]
	out = fmt.format(*out)
	filo.write(out)
	print (out)

 


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
# Lock	Exc	Num /mV	Bsplit /Hz	Err /Hz
"""
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

	# calc allan deviation of the Bsplit suing allantools!
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


	# read sidebands data
	if args.sbtag:
		sb = []
		for sbf in args.sbfile:
			# encoding =None avoid a warning
			# deletechars = set() avoid genfromtxt stripping away my "/" in the column names
			temp = genfromtxt(sbf, names=True, dtype=None, converters={i: ufloat_fromstr for i in arange(2,10)}, encoding=None, deletechars=set())
			sb += [temp]
			
		if len(sb)>1:
			sb1 =concatenate(sb, axis=0)
		else:
			sb1 = array(sb)

	
		
		# chose the only tag or the corresponding tag
		tag = args.sbtag[min(i, len(args.sbtag)-1)]

		# use only data with the right tag
		sb = sb1[where(sb1['tag'].astype(str) == tag)]
		
		
		tags[i]['sbtag'] = tag
		
		
		if len(sb)>0:
			conds[i]['U/Er'] = sb['U/Er'].mean()
			conds[i]['Chi'] = sb['Chi'].mean()	
			conds[i]['nz'] = sb['nz'].mean()
			conds[i]['fbridge/MHz'] = sb['fbridge/MHz'].mean()	
	
		else:
			print("WARNING: sideband tag not found")
		
		
		

		
	# other conds
	conds[i]['Toven/*C'] = args.Toven
	#conds[i]['Troom/*C'] = mTroom
	conds[i]['trabi/ms'] = args.trabi
	conds[i]['faom/Hz'] = faom
	
	
#	# save conditions config
#	config = SafeConfigParser()
#	# preserve case
#	config.optionxform=str
#	config.add_section(cond)
#	config.add_section(cond2)
#	
#	# set dafaults for user values
#	config.set(cond2, 'Ldepth', '')
#	config.set(cond2, 'Lfreq', '')
#	config.set(cond2, 'Toven', '')


#	# try to read existing user values
#	try:
#		config.read(procname + ".cfg")
#	except:
#		pass
#	
#	
#	# set new automatic user values
#	config.set(cond, 'exc', '{:S}'.format(mExc))
#	config.set(cond, 'num', '{:S}'.format(mNum))
#	config.set(cond, 'Bsplit', '{:S}'.format(mBsplit))
#	config.set(cond, 'Err', '{:.2uS}'.format(mErr))
#	#config.set(cond, 'Troom', '{:S}'.format(mTroom))
#	
#	with open(procname + ".cfg", 'wb') as cfg:
#		config.write(cfg)



	


	
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





#========================================
#get the difference between H and L loops
#========================================


if len(locks) == 2:
	# gets a new functions to interpolate the frequencies
	pff = [interp1d(d[:,0],d[:,3]) for d in pdata]

	# calculate the frequencies interpoalting
	diffe = pff[0](t) - pff[1](t)
	diff_rel = diffe/fYb


	data = column_stack((t,diffe, diff_rel))
	
	fmt="%.2f\t%.6f\t%e"
	header = "Time\tdiff\trel_diff"
	
	savetxt(basename + "_proc_diff.dat", data, fmt=fmt, header = header)


	figure()
	title("Diff. Shift")
	xlabel('Approx data point')
	ylabel('Frequency /Hz')
	plot((t-epoch0)/t0, diffe)
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
		
#	edg = vectorize(allantools.edf_greenhall)
#	edf = edg(0,2,m,N+1,True)

	#https://tf.boulder.nist.gov/general/pdf/666.pdf
	admax = adev * (edf/scipy.stats.chi2.ppf(0.1, edf))**.5
	admin = adev * (edf/scipy.stats.chi2.ppf(0.9, edf))**.5


	taufit = append(taut2, tottime)	


	data_out = column_stack((tau2, adev, aderr))
	fmt="%.2f\t%.3e\t%.3e"
	header = "Tau /s\toadev\toadev_unc"
	savetxt(basename + "_proc_diff_adev.dat", data_out, fmt=fmt, header = header)
	
	# simple fit
	white = adev*tau2**0.5
	#if len(white) > 10.:			
	#	white = white[where((tau2>10.) & (tau2<tottime/8.))]
	if len(white) > 4.:			
		white = white[where((tau2>10.) & (tau2<tottime/8.))]
	a = mean(white)	

	figure()
	title("Interl. Stability")
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

	

	# print some useful line to be copy pasted	
	condkeys = ['Exc', 'Num/mV', 'Bsplit/Hz', 'U/Er', 'Chi', 'nz']
	
	tit = [names[0] + x for x in condkeys] + [names[1] + x for x in condkeys] + ['Diff/Hz','relDiff']
	tit =  '#' + '\t'.join(tit) + '\n'
	print(tit)
	filo.write("\n\nA useful line to be copy pasted\n")
	filo.write(tit)
	
		
	out = [conds[0].get(x, None) for x in condkeys] + [conds[1].get(x, None) for x in condkeys] + [udiffe, udiff_rel]
	
	# if the values in the conds dictionary are sepcified or not change the format
	fmt = ["{:.2uS}" if isinstance(x, UFloat) else "{}" for x in out]
	fmt = "\t".join(fmt) + "\n"
	textout = fmt.format(*out)
	
	
	print(textout)
	filo.write(textout)






	# print the final message
	fmt = """
#Interleaved {} - {}
Difference = {:.2uS} Hz
Relative Diff = {:.2uS}

Density shift = {:.2uS} mHz/mV

White noise = {:.02} at 1 s fitted at {:.1f} s.


"""
	
	out = fmt.format(names[0],names[1],udiffe,udiff_rel, udiffe/(savenum[0]-savenum[1])*1e3, a, tottime)
	print (out)
	filo.write(out)
	

for i, lock in enumerate(locks):
	# process density tag
	
	if args.dself:
		density = udiffe/(conds[0]['Num/mV'] - conds[1]['Num/mV'])
		conds[i]['dcoeff/(Hz/mV)'] = density
		tags[i]['dtag'] = name
	elif args.dtag:
		# chose the only tag or the corresponding tag
		given = args.dtag[min(i, len(args.dtag)-1)]
		mores = given.split(',')
		dtags = [x[:19] for x in mores]
		dtag = ','.join(dtags)
		
		
		# can handle multiple density files
		density = []
		
		
		# recover density shift
		for tag in dtags:
			# serach for the tag + ".dat" file
			# and return the first occurrence
			# https://stackoverflow.com/questions/1724693/find-a-file-in-python
			def find(name, path):
				for root, dirs, files in os.walk(path):
					if name in files:
						return os.path.join(root, name)
			
			
			df = find(tag + ".dat", ".")
			if not df:	
				# search for the file i the Data folder (slow!)
				df = find(tag + ".dat", "..")
			
			if df:
				#cols to read
				# encoding =None avoid a warning
				# deletechars = set() avoid genfromtxt stripping away my "/" in the column names
				cols = ("Num/mV", "Diff/Hz")
				dd = genfromtxt(df, names=True, dtype=None, encoding=None, deletechars = set(), usecols=cols, converters={i: ufloat_fromstr for i in cols})
				
				diff = dd['Diff/Hz']
				nums = dd['Num/mV']
			
				density += [diff[0]/(nums[0] - nums[1])]
				
				 
				
				
		if density:
			density = array(density)
			weights = array([x.std_dev**-2 for x in density])
			conds[i]['dcoeff/(Hz/mV)'] = sum(weights*density)/sum(weights)
		
		tags[i]['dtag'] = dtag

	


# save data to the new .dat file
condkeys = ['Exc', 'Num/mV', 'Bsplit/Hz', 'Err/Hz', 'U/Er', 'Chi', 'nz', 'fbridge/MHz', 'Troom/*C', 'Toven/*C', 'trabi/ms', 'faom/Hz', 'dcoeff/(Hz/mV)']
tagkeys = ['sbtag', 'dtag']
tit = ['key                ', 'lock'] + condkeys + ['Diff/Hz','relDiff','white','toTime/s'] + tagkeys + ['Files']
#fmt="{}\t{}\t{}"+ "\t{:.2uS}"*(len(condkeys)-3) + "\t{}"*2 + "\t{:.2uS}"*2 + "\t{:.2}"*2 + "\t{}\n"
filo2.write("#" + "\t".join(tit) + "\n")


for i, lock in enumerate(locks):
	# note i wrote zero for the difference of the low lock
	# need to be adapeted if the interleaved locks are not 2
	# the get method on the conds dictionary default to None
	out = [name, names[i]] + [conds[i].get(x, None) for x in condkeys] + [udiffe*(i==0), udiff_rel*(i==0), a/sqrt(2), tottime] + [tags[i].get(x, None) for x in tagkeys] + [",".join(files)]
	
	# if the values in the conds dictionary are sepcified or not change the format
	
	fmt = ["{:.2uS}" if isinstance(x, UFloat) else "{}" for x in out]
	# special format for some thing
	fmt[-5] = "{:.3}"
	fmt[-4] = "{:.1f}"
	
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
