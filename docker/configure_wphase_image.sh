#!/bin/bash

#-------------------------------------------------------------------------------
# create required folders for wphase
#-------------------------------------------------------------------------------
install -g "$(id -gn)" -o "$(id -un)" -d \
    "$WPHASE_GREENS_FUNCTIONS_DIR" \
    "$WPHASE_WEB_OUTPUTS_ROOT"

# When developing by mounting the source directory into a container, we need
# to run build to compile the fortran code:
echo 'cd "$WPHASE_HOME" && pip install -e .' > ~/.bashrc
