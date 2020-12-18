#!/bin/bash

#-------------------------------------------------------------------------------
# script for playing with docker.
#
# $1: "build" (to build a wphase container with ready for wphase) or "run" to run
#     a container for development.
# $2: The command to run when launching the container. Not providing this will
#     launch a terminal session.
#-------------------------------------------------------------------------------

mode="$1"
cmd="$2"
if [ "$cmd" == "" ]; then
    cmd="bash"
else
    cmd="bash -lc '$cmd'"
fi

CONTAINER_NAME=wphase
BASEMAP_ARCHIVE_NAME=basemap.tar.gz
WPHASE_CONFIG_BUCKET_URL=
WPHASE_ROOT=/wphase
WPHASE_HOME="$WPHASE_ROOT"/eatws-wphase
WPHASE_WORKING_ROOT=/var/tmp/wphase
WPHASE_GREENS_FUNCTIONS_DIR="$WPHASE_ROOT"/greens
WPHASE_GREENS_FUNCTIONS_FILE=gfs_1.hdf5

if [ "$mode" == "build" ]; then
    if [ ! -f ./docker/"$BASEMAP_ARCHIVE_NAME" ]; then
        # get the archive of basemap, either from the specified S3 bucket or from github
        if [ "$WPHASE_CONFIG_BUCKET_URL" == "" ]; then
            curl -L https://github.com/matplotlib/basemap/archive/v1.1.0.tar.gz \
                -o docker/"$BASEMAP_ARCHIVE_NAME" \
                || exit 1
        else
            aws s3 cp \
                "$WPHASE_CONFIG_BUCKET_URL/$BASEMAP_ARCHIVE_NAME" \
                ./docker/"$BASEMAP_ARCHIVE_NAME"
        fi
    fi

    docker build \
        -t "$CONTAINER_NAME" \
        -f ./docker/Dockerfile \
        --build-arg BASEMAP_ARCHIVE_NAME="$BASEMAP_ARCHIVE_NAME" \
        --build-arg WPHASE_GREENS_FUNCTIONS_DIR="$WPHASE_GREENS_FUNCTIONS_DIR" \
        --build-arg WPHASE_GREENS_FUNCTIONS_FILE="$WPHASE_GREENS_FUNCTIONS_FILE" \
        --build-arg ARG_WPHASE_WEB_OUTPUTS_ROOT="$WPHASE_WORKING_ROOT"/web_outputs \
        --build-arg ARG_WPHASE_SAVED_DATASETS_ROOT="$WPHASE_WORKING_ROOT"/saved_datasets \
        --build-arg ARG_WPHASE_TEST_DATASETS_ROOT="$WPHASE_WORKING_ROOT"/test_data \
        --build-arg ARG_WPHASE_GREENS_FUNCTIONS="$WPHASE_GREENS_FUNCTIONS_DIR"/"$WPHASE_GREENS_FUNCTIONS_FILE" \
        ./docker

elif [ "$mode" == "run" ]; then
    docker run -it --rm \
        --mount type=bind,source=$HOME/wphase/greens,target="$WPHASE_GREENS_FUNCTIONS_DIR",readonly=true \
        --mount type=bind,source=`pwd`,target="$WPHASE_HOME" \
        -e "WPHASE_HOME=$WPHASE_HOME" \
        -e "WPHASE_HOST_NAME=0.0.0.0" \
        -p 5000:9999 \
        -w "$WPHASE_HOME"/api \
        "$CONTAINER_NAME" \
        $cmd
else
    echo 'First argument must be either "build" or "run".'
fi
