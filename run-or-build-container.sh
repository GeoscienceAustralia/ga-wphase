#!/bin/bash

#-------------------------------------------------------------------------------
# script for playing with docker.
#
# $1: "build" (to build a wphase container with ready for wphase) or "run" to run
#     a container for development.
# Rest of args: The command to run when launching the container. Not providing
#               this will launch a terminal session.
#-------------------------------------------------------------------------------

mode="$1"
shift
declare -a cmds=("$@")

IMAGE_NAME=wphase
if [ "$WPHASE_OUTPUT_DIR" = "" ]; then 
    WPHASE_OUTPUT_DIR="$HOME/wphase-outputs"
fi
HOST_WPHASE_OUTPUT_DIR="$WPHASE_OUTPUT_DIR"
HOST_WPHASE_GREENS_FUNCTIONS="$HOME/wphase/greens/gfs_1.hdf5"

mkdir -p $HOST_WPHASE_OUTPUT_DIR
if [ "$mode" == "build" ]; then
    docker build \
        -t "$IMAGE_NAME" \
        -f ./docker/Dockerfile \
        .

elif [ "$mode" == "run" ]; then
    # Run container for development
    if [ "${#cmds[@]}" != 0 ]; then
        cmd="${cmds[@]}"
        cmds=(bash -lc "$cmd") # this is an array, not a subshell!
    fi
    docker run -it --entrypoint /bin/bash --rm \
        --network=host \
        --mount type=bind,source=$HOST_WPHASE_GREENS_FUNCTIONS,target=/home/wphase/host_gf,readonly=true \
        --mount type=bind,source=`pwd`,target="/home/wphase/app",readonly=false \
        --mount type=bind,source="$HOST_WPHASE_OUTPUT_DIR",target=/outputs,readonly=false \
        --env WPHASE_GREENS_FUNCTIONS=/home/wphase/host_gf \
        --env WPHASE_OUTPUT_DIR=/outputs \
        "$IMAGE_NAME" \
        /home/wphase/reinstall-wphase-and-run "${cmds[@]}"
elif [ "$mode" == "run-wphase" ]; then
    # A single, production-style run
    docker run -it --rm \
        --network=host \
        --mount type=bind,source=$HOST_WPHASE_GREENS_FUNCTIONS,target=/home/wphase/host_gf,readonly=true \
        --mount type=bind,source="$HOST_WPHASE_OUTPUT_DIR",target=/outputs,readonly=false \
        --env WPHASE_GREENS_FUNCTIONS=/home/wphase/host_gf \
        --env WPHASE_OUTPUT_DIR=/outputs \
        "$IMAGE_NAME" \
        "${cmds[@]}"
else
    echo 'First argument must be either "build", "run" or "run-wphase".'
fi
