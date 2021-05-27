# GA W-Phase

Software used for running wphase moment tensor + centroid inversions in the
National Earthquake Alert Centre at Geoscience Australia.

Originally implemented by Roberto Benavente.

Based on the paper [*Source inversion of W phase: speeding up seismic tsunami
warning* by Kanamori and Rivera, 2008](https://doi.org/10.1111/j.1365-246X.2008.03887.x)


## Requirements

- Python 2.7+ or 3.6+ (with various PyPI dependencies)
- Green's functions (see the section below)
- Optional:
  - SeisComP with python libraries (to use the CLI script to push results into
    a SeisComP system/convert results to SCML)


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

From the same directory this README lives in:

```sh
pip install .
```

If you want to use it 'in place', use develop mode:

```sh
pip install -e .
```

If obspy has issues installing, you might need to install some dependencies first:

```sh
pip install wheel Cython numpy
pip install .
```

## Running wphase

Required environment variables:

- `WPHASE_OUTPUT_DIR`: absolute path to a directory in which W-Phase will write
  its outputs.
- `WPHASE_GREENS_FUNCTIONS`: absolute path to greens functions (either a
  directory or a hdf5 file, as described above)

You can then run W-Phase in a couple of ways:

- Directly from python, with basic results returned in a nested dict and
  visualizations written to a provided output directory. See e.g.
  [`test_runwphase`](tests/test_runwphase.py) for a minimal example.

- Using the command-line script [`wphase`](scripts/wphase), which takes its
  input as (a lot of!) command-line arguments and can send the basic results
  directly to a SeisComP messaging system, and the visualizations to an Amazon
  S3 bucket.

  For example, if you have an FDSN server at http://localhost:8081 and a
  SeisComP4+ messaging system at localhost:18180, the following script
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
    --host localhost:18180 \
    --debug
```

All the usual seiscomp configuration flags are available:

- `--console=1` to see errors on the command line (otherwise they are only
  written to the logfile!)
- `--debug` to see all log messages on the command line
- `--user` to set the messaging username (default is `gawphase`)
- `-g` to set the primary messaging group on which results are sent (default is
  `FOCMECH`)

You can get a full list of options by invoking `run-wphase` with `--help`.



## Using Docker

If you're running Linux with Docker and curl installed, you should be able to
build a docker container simply by running (as root)
`./run-or-build-container.sh build`.

By default, the image is configured to run as a non-root user with uid 1000;
hopefully this maps on to your host user for easy development. To change this
uid (or gid), provide the build-args `DOCKER_USER_UID` and/or
`DOCKER_USER_GID`.

### Developing on Docker

While the container image includes the complete W-Phase application, it can
easily be used for development by mounting in your working copy of the W-Phase
source.

When you launch an interactive bash session on the container (e.g. using the
`run-or-build-container.sh` script), if wphase is mounted at `/wphase` then it
will be automatically built and installed in development mode.

This is done for you if you use `./run-or-build-container.sh run`.

To use this script as-is, you'll need your greens functions to be stored at
`~/wphase/greens`; but you should be able to modify the script to suit your
needs.

To start some process inside the container:

```sh
./run-or-build-container.sh run <command>

# for example, to start a bash shell:
./run-or-build-container.sh run bash

# to run the test suite:
./run-or-build-container.sh run pytest

# to run the same example as above:
./run-or-build-container.sh run wphase \
    --evid ga2020yskxhe \
    --lat 5.2 \
    --lon 125.47 \
    --depth 42 \
    --time '2020-12-15T23:22:01Z' \
    --outputs /tmp/wphase-outputs \
    --server http://localhost:8081 \
    --host localhost:18180 \
    --debug
```

Since the `./run-or-build-container.sh` script uses [host
networking](https://docs.docker.com/network/host/), this last example would
still be able to contact the FDSN and Spread servers running on your host.
Alternatively, you could of course use some other Docker networking scheme and
swap out localhost for an appropriate hostname or IP.

### Using Docker in production

The container image is also ready for production use: it includes a fully-built
version of wphase and directly invokes the wphase script as its default
entrypoint default. It does *not* include the greens functions database (since
these are typically huge), so you'll still need to mount these; and you'll also
need to mount a directory to `/outputs`, where the output files will be stored.

You can use the `run-wphase` directive to invoke the production container,
automatically handling the directory mounting:

```
# to run the same example as above:
./run-or-build-container.sh run-wphase \
    --evid ga2020yskxhe \
    --lat 5.2 \
    --lon 125.47 \
    --depth 42 \
    --time '2020-12-15T23:22:01Z' \
    --server http://localhost:8081 \
    --host localhost:18180
```

In this example, the output files would then be available at
`./outputs/ga2020ykxhe/`.

## TODO

I extracted this code from our private repository and ran a quick clean-up
pass, but many more improvements are possible:

- Improve the Python API to return the results (focal mech,
  derived origin, etc) in ObsPy or SeisComP data structures instead of the
  current nested dictionaries?
- Look into replacing the libtau fortran code with obspy's travel time model.
    - When this code was first developed, apparently the latter was too slow; but
      things might have improved or there might be optimizations we can make.
    - I think this would allows us to relicense
      under the more permissive Apache license (GA's preference).
- General refactoring and code cleanup. I've done some, but there is a lot left
  to do (particularly in the `psi/` directory).
- Create a fully-integrated SeisComP plugin/daemon that watches the messaging
  bus for events matching certain criteria and runs automatically (like
  `scautomt` does).


## License

Copyright (C) 2021 Geoscience Australia

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
