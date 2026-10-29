# IT-Yb1 data scripts

Low-level processing scripts for the data of the optical clock IT-Yb1.

The two scripts are:

- `sproc.py` to process clock data: it provides analisys of interleaved measurements, instability and interpolation to common timetags
- `vband.py` to process sideband data: it provides sidebands fits and temperature vs depth analisys

## New development
New development of the sidebands fit has been moved to https://github.com/INRIM/large-lattice-model

## Basic usage

Make the module available using uv directly from github, for example:

    uvx --with "git+https://github.com/INRIM/IT-Yb1-data-scripts.git" ipython

Run the scripts in iptyhon with:

    run -m sproc -- --help
    run -m vband -- --help

Note the `-m` and the need to escape options with the first `--` 

## License

[MIT](https://opensource.org/licenses/MIT)

## Authors

(c) 2021-2026 Marco Pizzocaro - Istituto Nazionale di Ricerca Metrologica (INRIM)
