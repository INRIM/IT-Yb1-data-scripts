# IT-Yb1 data scripts

Low-level processing scripts for the data of the optical clock IT-Yb1.

The two scripts are:

- `sproc.py` to process clock data: it provides analisys of interleaved measurements, instability and interpolation to common timetags
- `vband.py` to process sideband data: it provides sidebands fits and temperature vs depth analisys


## Basic usage

In a Ipython console

`run sproc clock-data-file-to-be-processed` or

`run vband sideband-data-file-to-be-processed`

For more options see

`run sproc -h` or

`run vband -h`

## License

[MIT](https://opensource.org/licenses/MIT)

## Authors

(c) 2021 Marco Pizzocaro - Istituto Nazionale di Ricerca Metrologica (INRIM)
