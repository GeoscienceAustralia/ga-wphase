#!/bin/bash

set -e

pip3 install 'numpy>=1.16.5'

if [ "$1" == "dev" ]; then
    # When developing by mounting the source directory into a container, we need
    # to run build to compile the fortran code every time the container is
    # launched:
    echo 'cd "$WPHASE_HOME" && pip3 install -e .[aws,plotting]' >> ~/.bashrc
else
    echo 'cd "$WPHASE_HOME" && pip3 install -e .[aws,plotting]' > /reinstall-wphase
    echo '/reinstall-wphase && "$@"' > /reinstall-wphase-and-run
    chmod +x /reinstall-wphase /reinstall-wphase-and-run

    # For the production container, we can do the build immediately.
    /reinstall-wphase

    # And we want to preload cartopy features:
    python3 - <<PYTHON
from cartopy.feature import NaturalEarthFeature as NEF
for feature in (('physical', 'coastline'),
                ('physical', 'land'),
                ('physical', 'ocean'),
                ('cultural', 'admin_0_boundary_lines_land')):
    NEF(*feature, '110m').geometries() # triggers the download
PYTHON

fi
