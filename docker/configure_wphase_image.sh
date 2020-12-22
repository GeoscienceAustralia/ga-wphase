#!/bin/bash

set -e

if [ "$1" == "dev" ]; then
    # When developing by mounting the source directory into a container, we need
    # to run build to compile the fortran code every time the container is
    # launched:
    echo 'cd "$WPHASE_HOME" && pip install -e .' >> ~/.bashrc
else
    # For the production container, we can do the build immediately.
    cd "$WPHASE_HOME" && pip install -e .
fi
