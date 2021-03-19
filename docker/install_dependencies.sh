#!/bin/bash

set -e

#-------------------------------------------------------------------------------
# update system and install standard stuff
#-------------------------------------------------------------------------------
# for gcc, g++, ...
yum group install -y "Development Tools"
yum install -y curl python-devel.x86_64 python-setuptools tkinter
yum clean all

# get python major version
v="$(python --version 2>&1)"
v="${v#Python }"
v="${v:0:1}"
# install pip
[ "$v" == 2 ] \
    && url=https://bootstrap.pypa.io/pip/2.7/get-pip.py \
    || url=https://bootstrap.pypa.io/get-pip.py
curl -o get-pip.py $url
python get-pip.py
rm get-pip.py

# install EXACT VERSIONS of python dependencies.
# (since the py2 EOL, more and more PyPI libraries have broken compatibility, and
# the basemap version dependencies do NOT reflect this.)
pip install --no-cache-dir --force-reinstall --ignore-installed six==1.10.0 numpy==1.16.6
pip install --no-cache-dir --force-reinstall --ignore-installed -r /requirements.txt
rm /requirements.txt

#-------------------------------------------------------------------------------
# install basemap
#-------------------------------------------------------------------------------
LOCAL_BASEMAP_ARCHIVE_NAME="$LOCAL_INSTALL_FOLDER"/"$BASEMAP_ARCHIVE_NAME"
basemaptmpdir=/tmp/basemapinstall
mkdir "$basemaptmpdir"
tar -xzf "$LOCAL_BASEMAP_ARCHIVE_NAME" -C "$basemaptmpdir"

# install geos
cd "$basemaptmpdir"
cd basemap* # maybe the archive has a top-level folder
cd geos-*
./configure
make && make install

# install basemap itself
cd ..
python setup.py install

# cleanup
rm -rf "$LOCAL_BASEMAP_ARCHIVE_NAME"
cd && rm -rf "$basemaptmpdir"
