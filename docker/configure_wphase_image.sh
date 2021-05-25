#!/bin/bash

set -e

pip3 install --no-cache-dir --user 'numpy>=1.16.5'

echo 'cd "/home/wphase/app" && pip3 install --no-cache-dir --user -e .[aws,plotting] "$@"' > /home/wphase/reinstall-wphase
echo '/home/wphase/reinstall-wphase --quiet && "$@"' > /home/wphase/reinstall-wphase-and-run
chmod +x /home/wphase/reinstall-wphase /home/wphase/reinstall-wphase-and-run

# Do the build once to get dependencies installed.
/home/wphase/reinstall-wphase

# And we want to preload cartopy features:
python3 - <<PYTHON
from cartopy.feature import NaturalEarthFeature as NEF
for feature in (('physical', 'coastline'),
                ('physical', 'land'),
                ('physical', 'ocean'),
                ('cultural', 'admin_0_boundary_lines_land')):
    NEF(*feature, '110m').geometries() # triggers the download
PYTHON
