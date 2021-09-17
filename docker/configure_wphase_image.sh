#!/bin/bash

set -e

pip3 install --no-cache-dir --user 'numpy>=1.16.5'

echo 'cd "/home/wphase/app" && pip3 install --no-cache-dir --user -e .[aws,plotting] "$@"' > /home/wphase/reinstall-wphase
echo '/home/wphase/reinstall-wphase --quiet && "$@"' > /home/wphase/reinstall-wphase-and-run
chmod +x /home/wphase/reinstall-wphase /home/wphase/reinstall-wphase-and-run

# Do the build once to get dependencies installed.
/home/wphase/reinstall-wphase

# Configure cartopy to fetch downloads from S3:
usersitedir=$(python -c 'from __future__ import print_function; import site; print(site.getusersitepackages())')
mkdir -p "$usersitedir"/cartopy_userconfig
cat > "$usersitedir"/cartopy_userconfig/__init__.py <<'PYTHON'

_SOURCE_TEMPLATE = 'https://naturalearth.s3.amazonaws.com/{resolution}_{category}/ne_{resolution}_{name}.zip'

def update_config(config):
    """Configures cartopy to download NaturalEarth shapefiles from S3 instead
    of naciscdn."""
    from cartopy.io.shapereader import NEShpDownloader
    target_path_template = NEShpDownloader.default_downloader().target_path_template
    downloader = NEShpDownloader(url_template=_SOURCE_TEMPLATE,
                                 target_path_template=target_path_template)
    config['downloaders'][('shapefiles', 'natural_earth')] = downloader

PYTHON

# Finally, pre-download the polygons we need:
python3 - <<PYTHON
from cartopy.feature import NaturalEarthFeature as NEF
for feature in (('physical', 'coastline'),
                ('physical', 'land'),
                ('physical', 'ocean'),
                ('cultural', 'admin_0_boundary_lines_land')):
    NEF(*feature, '110m').geometries() # triggers the download
PYTHON
