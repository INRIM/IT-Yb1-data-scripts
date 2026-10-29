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


import argparse
import fnmatch
import os
import re
import sys
from datetime import datetime

import large_lattice_model as llm
import scipy.optimize as opt
from matplotlib import cm
from matplotlib.pyplot import *
from numpy import *
from uncertainties import correlated_values, ufloat, ufloat_fromstr, unumpy

__version__ = "1.1"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process IT-Yb1 sideband data based on a Born-Oppenheimer model",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("infile", nargs="*", help="File or list of files to be processed", default=[])
    parser.add_argument("--dir", help="Directory for storing results", default="./Vbands")
    parser.add_argument("--log", help="Filename where to save the results. Results are appended", default="all.dat")
    parser.add_argument("--tag", nargs="+", help="Tag or list of tags for each file", default="-")
    parser.add_argument(
        "--num", type=str, nargs="*", help="Search files from num range instead of full name, e.g. 1-10"
    )

    parser.add_argument("--Dscale", type=float, help="Scale for the initial guess of D", default=0.9)
    parser.add_argument("--Tzscale", type=float, help="Scale for the initial guess of Tz", default=0.4)
    parser.add_argument("--Trscale", type=float, help="Scale for the initial guess of Tr", default=0.4)
    parser.add_argument("--Tz", type=float, help="Fix the initial guess of Tz, then Tzscale is ignored", default=None)
    parser.add_argument("--Tr", type=float, help="Fix the initial guess of Tr, then Trscale is ignored", default=None)
    parser.add_argument("--fac", type=float, help="Quality of integration parameter", default=10)

    parser.add_argument("--nofit", action="store_true", help="Do not fit for all processed files")
    parser.add_argument("--Tzlin", action="store_true", help="Fit a linear Tz vs D")
    parser.add_argument("--Trlin", action="store_true", help="Fit a linear Tr vs D")
    # parser.add_argument('--Nfit', action='store_true', help='Show N vs U for processed files')

    command = sys.argv
    args = parser.parse_args()

    if not os.path.exists(args.dir):
        os.makedirs(args.dir)
    figdir = os.path.join(args.dir, "Figures")
    if not os.path.exists(figdir):
        os.makedirs(figdir)

    ioff()
    close("all")


Er = llm.settings.Er
kB = llm.kB
vrec = Er / llm.h


if __name__ == "__main__":
    allfigname = os.path.join(args.dir, os.path.splitext(args.log)[0] + ".png")
    alltxtname = os.path.join(args.dir, args.log)
    fitbasename = os.path.splitext(args.log)[0]

    with open(alltxtname, "a") as alltxt:
        labels = ["Ac", "fc/Hz", "wc/Hz", "A", "D/Er", "Tz/K", "Tr/K", "Nat/mV"]
        tit = "#file" + " " * 32 + "\ttag\tl/mV\tspin\t" + "\t".join(labels) + "\n"
        logfmt = "{}\t{}\t{}\t{}" + "\t{:.2uS}" * len(labels) + "\n"
        alltxt.write(tit)

    if args.num:

        def myint(x):
            return int(x) if x else None

        lim_list = [tuple(map(myint, i.split("-"))) for i in args.num]
        corrected_lim_list = [(x[0], x[0] + 1) if len(x) == 1 else (x[0], x[1] + 1) for x in lim_list]
        num_list = concatenate([arange(*x) for x in corrected_lim_list])
        tag_list = [f"{x:03d}" for x in num_list]

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
                print(f"WARNING: I have not found a file starting with {tag}")
            else:
                args.infile = df

    Nf = len(args.infile)

    for i, filname in enumerate(args.infile):
        with open(filname) as fil:
            basename = os.path.basename(fil.name)[:-4]
            numname = os.path.basename(fil.name)[:3]
            figname = os.path.join(figdir, basename + ".png")
            txtname = os.path.join(args.dir, basename + ".txt")

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
                ll += re.findall(r"\bl\s?=\s?(\".*?\"|\S*)", line)
                ss += re.findall(r"\bs\s?=\s?(\".*?\"|\S*)", line)

            if ll:
                l = ",".join(ll)
            else:
                l = nan
            if ss:
                spin = ",".join(ss)
            else:
                spin = nan

            tag = args.tag[min(i, len(args.tag) - 1)]

            print("file: " + fil.name)
            print("tag: " + tag)
            print("l: ", l)
            print("spin: ", spin)

            # fitting procedure
            # =================

            data = genfromtxt(fil, usecols=(2, 3, 4, 5, 6, 7), skip_footer=1, skip_header=1)

            freq = data[:, 0]
            freq = around(freq, 4)
            funique = unique(freq)  # find scan frequency values

            # average data
            ave = []
            for f in funique:
                ave += [mean(data[where(freq == f)], axis=0)]

            ave = array(ave)
            # freq = funique

            # unpack!
            freq, g, b, e, num, exc = ave.T
            G = g - b
            E = e - b

            # lets get down some scales from labedit exc
            # medex = mean(exc)
            # test = where(exc > medex)[0]

            # fmi = freq[test[0]]
            # fma = freq[test[-1]]
            # fce = (fmi+fma)/2
            # fxc = (fma-fmi)/2

            # nowadays sidebands are always pretty centered
            # so, no need to test where the sidebands are
            fmi = amin(freq)
            fma = amax(freq)
            fce = (fmi + fma) / 2
            fxc = (fma - fmi) / 2

            # now recalculate excitation eta on the red side of the spectrum
            nscale = 5e3  # Hz
            mask = funique < fce + nscale

            def findeta(eta, g, e):
                return std(g + e / eta)

            res = opt.minimize(findeta, 0.97, args=(G[mask], E[mask]))
            eta = res.x[
                0
            ]  # we found abug for wich the first blue pulse gives less signal than the third, so eta can be > 1 !?!

            N = G + E / eta
            baseN = median(N[mask])

            # now use the new data
            # this substantially fit just th excitation
            exc = E / baseN

            # fit carrier
            cscale = 10e3  # Hz
            mask = abs(funique - fce) < cscale

            cfreq = freq[mask]
            cexc = exc[mask]

            a = amax(cexc) * 2
            fc = cfreq[argmax(cexc)]

            copt, ccov = opt.curve_fit(llm.fit_lorentzian, cfreq, cexc, p0=[a, fc + 100, cscale / 10])

            # save carrier width and translate to center
            w0 = copt[2]
            freq -= copt[1]
            copt[1] = 0

            # remove carrier for future fits
            sexc = exc - llm.fit_lorentzian(freq, *copt)

            # guess sbands parameters
            A = sum(cexc) * median(diff(freq)) / w0
            D = fxc**2 / (4 * vrec**2) * args.Dscale
            # D=58.34
            if args.Tz:
                Tz = args.Tz
            else:
                Tz = D * Er / kB * args.Tzscale

            if args.Tr:
                Tr = args.Tr
            else:
                Tr = D * Er / kB * args.Trscale

            # fit
            fit_sbands = llm.get_fit_sidebands(w0, dn=1, E_max=0, fac=args.fac)
            sopt, scov = opt.curve_fit(
                fit_sbands, freq, sexc, p0=[A, D, Tz, Tr], bounds=([0.0, 0.0, 1e-7, 1e-7], [300, 2 * D, 1e-4, 1e-4])
            )

            # Plot fit
            ff = linspace(amin(freq), amax(freq), 400)
            fit = fit_sbands(ff, *sopt) + llm.fit_lorentzian(ff, *copt)

            # single plot
            figure("comp", figsize=(12.8, 4.8))
            plot(freq, exc, "o", label=numname, color=cm.viridis(i / Nf))
            plot(ff, fit, color=cm.viridis(i / Nf))
            grid(True)
            ylabel("Excitation probability")
            xlabel("Frequency detuning /Hz")
            legend(loc=0)

            # for each file
            figure(figsize=(12.8, 4.8))
            title(basename)
            plot(freq, exc, "o")
            plot(ff, fit)
            grid(True)
            ylabel("Excitation probability")
            xlabel("Frequency detuning /Hz")
            savefig(figname)
            # pause(0.001) # matplotlib need a pause to catch up -- see https://stackoverflow.com/questions/28269157/plotting-in-a-non-blocking-way-with-matplotlib
            # for many figure do not show
            if Nf > 9:
                close()

            # extract parameters
            c_par = correlated_values(copt, ccov)
            s_par = correlated_values(sopt, scov)

            Natoms = ufloat(mean(data[:, -2]), std(data[:, -2]))

            print(f"A = {s_par[0]:.2uS}")
            print(f"D = {s_par[1]:.2uS} Er")
            print(f"vz = {s_par[1] ** 0.5 * vrec * 2 * 1e-3:.2uS} kHz")
            print(f"Tz = {s_par[2]:.2uS} K")
            print(f"Tr = {s_par[3]:.2uS} K")
            print(f"wc = {c_par[2]:.2uS} Hz")
            print(f"Nat = {Natoms:.2uS} mV")
            print(f"eta = {eta:.3}")
            print("\n")

            # save data in tabular form on the common file
            with open(alltxtname, "a") as alltxt:
                alltxt.write(logfmt.format(basename, tag, l, spin, *c_par, *s_par, Natoms))

        if Nf > 1:
            figure("comp")
            savefig(allfigname)

        if Nf > 0:
            with open(alltxtname, "a") as alltxt:
                alltxt.write("# Above data generated on: " + str(datetime.now()) + "\n")
                alltxt.write("# With the command line: " + " ".join(command) + "\n")
                # alltxt.write("# Script version: " +  __version__ + '\n')

        if not args.nofit and Nf > 1:
            # deletechars = set() avoid genfromtxt stripping away my "/" in the column names
            data = genfromtxt(
                alltxtname,
                names=True,
                dtype=None,
                converters={i: ufloat_fromstr for i in arange(4, 4 + len(labels))},
                encoding=None,
                deletechars=set(),
            )

            D = unumpy.nominal_values(data["D/Er"])
            uD = unumpy.std_devs(data["D/Er"])

            Tz = unumpy.nominal_values(data["Tz/K"])
            uTz = unumpy.std_devs(data["Tz/K"])

            Tr = unumpy.nominal_values(data["Tr/K"])
            uTr = unumpy.std_devs(data["Tr/K"])

            Na = unumpy.nominal_values(data["Nat/mV"])
            uNa = unumpy.std_devs(data["Nat/mV"])

            l = data["l/mV"].astype(float)

            # quadratic fit
            def fun(x, a, b):
                return a * (x * Er) / kB + b * x**2 * Er / kB

            eps = 1e-6

            if not args.Tzlin:
                zopt, zcov = opt.curve_fit(fun, D, Tz, sigma=uTz, p0=[0.2, 0.001])
            else:
                zopt, zcov = opt.curve_fit(fun, D, Tz, sigma=uTz, p0=[0.2, 0.0], bounds=([0, -eps], [10.0, eps]))

            if not args.Trlin:
                ropt, rcov = opt.curve_fit(fun, D, Tr, sigma=uTr, p0=[0.2, 0.001])
            else:
                ropt, rcov = opt.curve_fit(fun, D, Tr, sigma=uTr, p0=[0.2, 0.0], bounds=([0, -eps], [10.0, eps]))

            z_par = correlated_values(zopt, zcov)
            r_par = correlated_values(ropt, rcov)

            print("Tfit\n")
            print(f"Az = {z_par[0]:.2uS}")
            if not args.Tzlin:
                print(f"Bz = {z_par[1]:.2uS}")
            print(f"Ar = {r_par[0]:.2uS}")
            if not args.Trlin:
                print(f"Br = {r_par[1]:.2uS}")

            figure()
            scale = 1e6

            for i, x in enumerate(set(data["spin"])):
                mask = where(data["spin"] == x)
                errorbar(
                    D[mask],
                    Tz[mask] * scale,
                    yerr=uTz[mask] * scale,
                    xerr=uD[mask],
                    fmt="o",
                    color=cm.tab20(1 - i),
                    label=f"Tz, s={x}",
                )
                errorbar(
                    D[mask],
                    Tr[mask] * scale,
                    yerr=uTr[mask] * scale,
                    xerr=uD[mask],
                    fmt="d",
                    color=cm.tab20(3 - i),
                    label=f"Tr, s={x}",
                )

            uu = linspace(0, max(D) * 1.1, 100)

            plot(uu, fun(uu, *zopt) * scale, color=cm.tab20(4))
            plot(uu, fun(uu, *ropt) * scale, color=cm.tab20(6))

            legend(loc=0)
            xlabel("D/ Er")
            ylabel("T /µK")
            grid(True)

            left, right = xlim()
            if right < 550:
                xlim(left, 550)

            tfitname = os.path.join(args.dir, fitbasename + "-tfit.png")
            savefig(tfitname)

            figure()
            scale = 1

            for i, x in enumerate(set(data["spin"])):
                mask = where(data["spin"] == x)
                errorbar(
                    D[mask],
                    Na[mask] * scale,
                    yerr=uNa[mask] * scale,
                    xerr=uD[mask],
                    color=cm.tab20(1 - i),
                    fmt="o",
                    label=f"Natoms s={x}",
                )

            legend(loc=0)
            xlabel("D/ Er")
            ylabel("Natoms /mV")
            grid(True)
            xlim(0, None)

            left, right = xlim()
            if right < 550:
                xlim(left, 550)

            nfitname = os.path.join(args.dir, fitbasename + "-N.png")
            savefig(nfitname)

            # also fit tags

            def fun(x, A):
                return A * x

            figure()

            errorbar(l, D, yerr=uD, fmt="o")

            mask = ~isnan(l)
            if sum(mask) > 2:
                lopt, lcov = opt.curve_fit(fun, l[mask], D[mask], sigma=uD[mask])
                ll = linspace(0, max(l[mask]) * 1.1, 100)
                plot(ll, fun(ll, *lopt))
            else:
                lopt, lcov = zeros(1), atleast_2d(ones(1))

            xlabel("l /mV")
            ylabel("D/ Er")

            grid(True)

            Dfitname = os.path.join(args.dir, fitbasename + "-Dfit.png")
            savefig(Dfitname)

            l_par = correlated_values(lopt, lcov)

            print("Dfit\n")
            print(f"k = {l_par[0]:.2uS}")

            Ttxtname = os.path.join(args.dir, fitbasename + "-fit.dat")

            with open(Ttxtname, "w") as Ttxt:
                Tlabels = ["Az", "Bz", "Ar", "Br", "k"]
                Ttit = "#" + "\t".join(Tlabels) + "\n"
                Tlogfmt = "\t".join(["{:.2uS}"] * len(Tlabels)) + "\n"

                Ttxt.write(Ttit)
                Ttxt.write("#\n# Temperature fits\n")
                Ttxt.write("# kB*T/Er = A*(D/Er) + B*(D/Er)**2\n")
                Ttxt.write("#\n# Depth fit\n")
                Ttxt.write("# D/Er = k*l/mV\n#\n")
                Ttxt.write(Tlogfmt.format(*z_par, *r_par, *l_par))

    show(block=False)  # note i prefer to be ioff for this script
