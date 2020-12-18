# GA W-Phase

Software used for running wphase moment tensor + centroid inversions in the
National Earthquake Alert Centre.

Based on the paper [*Source inversion of W phase: speeding up seismic tsunami
warning* by Kanamori and Rivera, 2008](https://doi.org/10.1111/j.1365-246X.2008.03887.x)

## Requirements

- Python 2.7 (with various libraries)
- Greens functions. TODO: Do we have some reference explaining the required format?
  We use a single-file database called `gfs_1.hdf5`.
- Optional: SeisComP3 with python libraries (to use the fully integrated script)


## Direct Installation

First, install the matplotlib basemap toolkit following [their
instructions](https://github.com/matplotlib/basemap). The rest of the python
dependencies can be installed automatically from PyPI.

From the same directory as this README lives in:

```sh
sudo pip install .
```

If you want to use it 'in place', use develop mode:

```sh
sudo pip install -e .
```

## Running wphase

Required environment variables:

- `WPHASE_HOME`: absolute path of a working directory forwphase, e.g. `/home/you/wphase`
- `WPHASE_GREENS_FUNCTIONS`: absolute path to greens functions (either a hdf5
  file or a directory in SAC structure)

You can then run W-Phase in a couple of ways:

- Directly from python, with basic results returned in a nested dict and
  visualizations written to a provided output directory. See e.g.
  [`test_runner_fdsn`](tests/test_runner_fdsn.py) for a minimal example.
- Using the command-line script [`wphase`](scripts/wphase), which takes its
  input as (a lot of!) command-line arguments and can send the basic results
  directly to a SeisComP3 messaging system, and the visualizations to an Amazon
  S3 bucket.  (Note that the docker container below does not have SeisComP
  included, so it won't work there.)

  TODO: Provide an example of using this script to push data into sc3.


## Using Docker

If you're running Linux with Docker and curl installed, you should be able to
build a docker container for development simply by running (as root)
`./run-or-build-container.sh build`.

The resulting container is designed for development use; it does not have wphase
installed, just all the dependencies. I mount the current folder at */wphase*.
when launching the container (which is done for you if you use
`run-or-build-container.sh`).

When you launch an interactive bash session on the container (e.g. using the
`run-or-build-container.sh` script), if wphase is mounted at `/wphase` then it
will be automatically built and installed in development mode.

To use this script as-is, you'll need your greens functions to be stored at
`~/wphase/greens`; but you should be able to modify the script to suit your
needs.


To start some process inside the container:

```sh
sudo -E ./run-or-build-container.sh run <command>

# for example, to start a bash shell:
sudo -E ./run-or-build-container.sh run bash

# to run the test suite:
sudo -E ./run-or-build-container.sh run pytest
```


## TODO

I extracted this code from our private repository and ran a quick clean-up
pass, but many more improvements are possible:

- Migrate to Python 3 and add SeisComP 4+ support to [scripts/wphase](scripts/wphase).
- Improve the Python API to return the results (focal mech,
  derived origin, etc) in ObsPy or SeisComP data structures instead of the
  current nested dictionaries.
- Migrate from basemap to cartopy.
- Look into replacing the libtau fortran code with obspy's travel time model.
    - When this code was first developed, apparently the latter was too slow; but
      things might have improved or there might be optimizations we can make.
- If both of the two above points are achieved, I think we could relicense
  under the more permissive Apache license (GA's preference).
- Create a fully-baked docker image for production use.
- Strip out any remaining dead code, refactor and improve documentation of the
  remainder.
- Create a fully-integrated SeisComP plugin/daemon that watches the messaging
  bus for events matching certain criteria and runs automatically (like
  `scautomt` does).

## License

Copyright (C) 2020 Geoscience Australia

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
