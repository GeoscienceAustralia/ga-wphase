# GA W-Phase

Software used for running wphase moment tensor + centroid inversions in the
National Earthquake Alert Centre at Geoscience Australia.

Originally implemented by Roberto Benavente.

Based on the paper [*Source inversion of W phase: speeding up seismic tsunami
warning* by Kanamori and Rivera, 2008](https://doi.org/10.1111/j.1365-246X.2008.03887.x)


## Requirements

- Python 2.7 (with various libraries from PyPI).
  - Python 3 is now also tentatively supported.
- Green's functions (see the section below)
- Optional:
  - SeisComP3 with python libraries (to use the CLI script to push results into
    a SeisComP3 system/convert results to SCML)
  - The basemap matplotlib toolkit (to plot maps showing station distribution
    and illustrating the grid search)


## Green's Functions

The Green's Function database used for W-Phase at GA is the same as those used
by Kanamori and Rivera in their original implementation - see
[here](http://wphase.unistra.fr/wiki/doku.php/wphase#green_s_function_databases)
for their information regarding these. You can contact them to attempt to
obtain a copy of their original database; or attempt to generate a new one
yourself using the PREM model. The database is stored as a directory structure
containing SAC files, following this template:

    H{DEPTH}/{MT_COMPONENT}/GF.{DISTANCE}.SY.LH{C}.SAC

Where:

- `DEPTH` is the depth in km, which must end in `.5`
- `MT_COMPONENT` is the moment tensor component, one of `PP`, `RP`, `RR`, `RT`,
  `TP`, `TT` where `P` is &phi;, `T` is &theta;, `R` is r.
- `DISTANCE` is the distance in 10ths of a degree, padded to 4 characters.
  (In the current implementation, your database *must* include every 10th of a
  degree from 0 to 90 degrees.)
- `C` is the waveform component

For example, the vertical (`Z`) Green's function for M<sub>&phi;&phi;</sub>=1
(`PP`) at a depth of 500.5km and distance of 90 degrees is stored at
`H500.5/PP/GF.0900.SY.LHZ.SAC`.

As alternative to this very large file structure, an equivalent structure can be
stored inside a HDF5 file.


## Direct Installation

First, install the matplotlib basemap toolkit following [their
instructions](https://github.com/matplotlib/basemap).  You can skip this step
if you don't want to use the map plotting functionality.

The rest of the python dependencies can be installed automatically from PyPI.

From the same directory as this README lives in:

```sh
pip install .
```

If you want to use it 'in place', use develop mode:

```sh
pip install -e .
```

If obspy has issues installing, you might need to install numpy first:

```sh
pip install wheel
pip install numpy
pip install .
```

## Running wphase

Required environment variables:

- `WPHASE_HOME`: absolute path of a working directory for wphase, e.g.
  `/home/you/wphase`
- `WPHASE_GREENS_FUNCTIONS`: absolute path to greens functions (either a
  directory or a hdf5 file, as described above)

You can then run W-Phase in a couple of ways:

- Directly from python, with basic results returned in a nested dict and
  visualizations written to a provided output directory. See e.g.
  [`test_runwphase`](tests/test_runwphase.py) for a minimal example.

- Using the command-line script [`wphase`](scripts/wphase), which takes its
  input as (a lot of!) command-line arguments and can send the basic results
  directly to a SeisComP3 messaging system, and the visualizations to an Amazon
  S3 bucket.

  For example, if you have an FDSN server at http://localhost:8081 and a
  SeisComP3 messaging system (spread) at localhost:4803, the following script
  would attempt to solve for the centroid moment tensor of an event at the
  given location and time, and if successful, send the results to seiscomp
  under the given event ID.

```sh
wphase \
    --evid ga2020yskxhe \
    --lat 5.2 \
    --lon 125.47 \
    --depth 42 \
    --time '2020-12-15T23:22:01Z' \
    --outputs /tmp/wphase-outputs \
    --server http://localhost:8081 \
    --host localhost:4803
```

## Using Docker

If you're running Linux with Docker and curl installed, you should be able to
build docker containers for development and production simply by running (as root)
`./run-or-build-container.sh build`. By default, this will build two very
similar containers: `wphase` and `wphase-dev`.

### `wphase-dev`

This container is designed for development use; it does not have
wphase installed, just all the dependencies. I mount the current folder at
*/wphase*.  when launching the container (which is done for you if you use
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

# to run the same example as above:
sudo -E ./run-or-build-container.sh run wphase \
    --evid ga2020yskxhe \
    --lat 5.2 \
    --lon 125.47 \
    --depth 42 \
    --time '2020-12-15T23:22:01Z' \
    --outputs /tmp/wphase-outputs \
    --server http://localhost:8081 \
    --host localhost:4803
```

Since the `./run-or-build-container.sh` script uses [host
networking](https://docs.docker.com/network/host/), this last example would
still be able to contact the FDSN and Spread servers running on your host.
Alternatively, you could of course use some other Docker networking scheme and
swap out localhost for an appropriate hostname or IP.

### `wphase-prod`

This container is designed for production use: it includes a fully-built
version of wphase and directly invokes the wphase script by default. It does
*not* include the greens functions database (since these are typically huge),
so you'll still need to mount these; and you'll also need to mount a directory
to `/outputs`, where the output files will be stored.

You can use the `run-wphase` directive to invoke the production container,
automatically handling the directory mounting:

```
# to run the same example as above:
sudo -E ./run-or-build-container.sh run-wphase \
    --evid ga2020yskxhe \
    --lat 5.2 \
    --lon 125.47 \
    --depth 42 \
    --time '2020-12-15T23:22:01Z' \
    --server http://localhost:8081 \
    --host localhost:4803
```

In this example, the output files would then be available at
`./outputs/ga2020ykxhe/`.

All the usual seiscomp configuration flags are available:

- `--debug` to see log output
- `--user` to set the messaging username (default is `gawphase`)
- `-g` to set the primary messaging group on which results are sent (default is
  `FOCMECH`)

You can get a full list of options by invoking `run-wphase` with `--help`.


## TODO

I extracted this code from our private repository and ran a quick clean-up
pass, but many more improvements are possible:

- Add SeisComP 4+ support to [scripts/wphase](scripts/wphase).
- Improve the Python API to return the results (focal mech,
  derived origin, etc) in ObsPy or SeisComP data structures instead of the
  current nested dictionaries?
- Migrate from basemap to cartopy.
- Look into replacing the libtau fortran code with obspy's travel time model.
    - When this code was first developed, apparently the latter was too slow; but
      things might have improved or there might be optimizations we can make.
- If both of the two above points are achieved, I think we could relicense
  under the more permissive Apache license (GA's preference).
- Strip out any remaining dead code, then refactor what's left and improve the
  documentation.
- Create a fully-integrated SeisComP plugin/daemon that watches the messaging
  bus for events matching certain criteria and runs automatically (like
  `scautomt` does).
- Reimplement the production container on top of a slimmer base OS to decrease
  the filesize.


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
