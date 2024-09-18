# IT-Yb1 data scripts

Low-level processing scripts for the data of the optical clock IT-Yb1.

The two scripts are:

- `sproc.py` to process clock data: it provides analisys of interleaved measurements, instability and interpolation to common timetags
- `vband.py` to process sideband data: it provides sidebands fits and temperature vs depth analisys

## New development
New development of the sidebands fit has been moved to https://github.com/INRIM/large-lattice-model

## Basic usage

In a Ipython console

`run sproc clock-data-file-to-be-processed` or

`run vband sideband-data-file-to-be-processed`

For more options see

`run sproc -h` or

`run vband -h`

## Requirements

The scripts require the following python packages:
- `numpy`, `scipy`, `matplotlib`
- `uncertainties`
- `allantools`

`vband.py` requires also:
- `pygsl`

To instal gsl (under Linux Mint):

`apt-get install libgsl-dev`

`pip3 instal pygsl`


## License

[MIT](https://opensource.org/licenses/MIT)

## Authors

(c) 2021 Marco Pizzocaro - Istituto Nazionale di Ricerca Metrologica (INRIM)
