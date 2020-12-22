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
cmd="$@"

DEV_CONTAINER_NAME=wphase-dev
PROD_CONTAINER_NAME=wphase
BASEMAP_ARCHIVE_NAME=basemap.tar.gz
WPHASE_CONFIG_BUCKET_URL=
WPHASE_ROOT=/wphase
WPHASE_HOME="$WPHASE_ROOT"/app
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

    build() {
        target="$1"
        name="$2"
        docker build \
            -t "$name" \
            -f ./docker/Dockerfile \
            --build-arg BASEMAP_ARCHIVE_NAME="$BASEMAP_ARCHIVE_NAME" \
            --build-arg WPHASE_GREENS_FUNCTIONS_DIR="$WPHASE_GREENS_FUNCTIONS_DIR" \
            --build-arg WPHASE_GREENS_FUNCTIONS_FILE="$WPHASE_GREENS_FUNCTIONS_FILE" \
            --build-arg ARG_WPHASE_WEB_OUTPUTS_ROOT="$WPHASE_WORKING_ROOT"/web_outputs \
            --build-arg ARG_WPHASE_SAVED_DATASETS_ROOT="$WPHASE_WORKING_ROOT"/saved_datasets \
            --build-arg ARG_WPHASE_TEST_DATASETS_ROOT="$WPHASE_WORKING_ROOT"/test_data \
            --build-arg ARG_WPHASE_GREENS_FUNCTIONS="$WPHASE_GREENS_FUNCTIONS_DIR"/"$WPHASE_GREENS_FUNCTIONS_FILE" \
            --build-arg WPHASE_HOME="$WPHASE_HOME" \
            --target "$target" \
            .
    }
    build dev $DEV_CONTAINER_NAME
    build prod $PROD_CONTAINER_NAME

elif [ "$mode" == "run" ]; then
    if [ "$cmd" != "" ]; then
        cmd="bash -lc '$cmd'"
    fi
    docker run -it --rm \
        --mount type=bind,source=$HOME/wphase/greens,target="$WPHASE_GREENS_FUNCTIONS_DIR",readonly=true \
        --mount type=bind,source=`pwd`,target="$WPHASE_HOME" \
        --network=host \
        -e "WPHASE_HOME=$WPHASE_HOME" \
        "$DEV_CONTAINER_NAME" \
        $cmd
elif [ "$mode" == "run-wphase" ]; then
    mkdir outputs
    docker run -it --rm \
        --mount type=bind,source=$HOME/wphase/greens,target="$WPHASE_GREENS_FUNCTIONS_DIR",readonly=true \
        --mount type=bind,source=`pwd`/outputs,target="/outputs",readonly=false \
        --network=host \
        -e "WPHASE_HOME=$WPHASE_HOME" \
        "$PROD_CONTAINER_NAME" \
        --outputs outputs \
        $cmd
else
    echo 'First argument must be either "build", "run" or "run-wphase".'
fi
