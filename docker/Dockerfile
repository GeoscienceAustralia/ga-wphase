################################################
# Base image with all dependencies and seiscomp
FROM ubuntu:20.04 AS base

# SeisComP + W-Phase dependencies from apt:
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        build-essential \
        curl \
        gfortran \
        libboost-filesystem1.71.0 \
        libboost-iostreams1.71.0 \
        libboost-program-options1.71.0 \
        libboost-regex1.71.0 \
        libboost-system1.71.0 \
        libboost-thread1.71.0 \
        libgeos-dev \
        libmariadb3 \
        libmysqlclient21 \
        libncurses5 \
        libpq5 \
        libproj-dev \
        libssl1.1 \
        libxml2 \
        proj-bin \
        proj-data \
        python3-dev \
        python3-pip \
        python3-tk \
        python3-twisted \
        python3-dateutil \
    && rm -rf /var/lib/apt/lists/*


# Install seiscomp3:
ARG SEISCOMP_TARBALL_URL=https://www.seiscomp.de/downloader/seiscomp-4.5.0-ubuntu20.04-x86_64.tar.gz
RUN cd /opt/ && curl -s $SEISCOMP_TARBALL_URL | tar xfz -

ARG LOCAL_INSTALL_FOLDER=/tmp
ARG BASEMAP_ARCHIVE_NAME
ARG WPHASE_GREENS_FUNCTIONS_DIR
ARG WPHASE_GREENS_FUNCTIONS_FILE
ARG ARG_WPHASE_WEB_OUTPUTS_ROOT
ARG ARG_WPHASE_SAVED_DATASETS_ROOT
ARG ARG_WPHASE_TEST_DATASETS_ROOT
ARG ARG_WPHASE_GREENS_FUNCTIONS
ENV WPHASE_WEB_OUTPUTS_ROOT="$ARG_WPHASE_WEB_OUTPUTS_ROOT"
ENV WPHASE_SAVED_DATASETS_ROOT="$ARG_WPHASE_SAVED_DATASETS_ROOT"
ENV WPHASE_TEST_DATASETS_ROOT="$ARG_WPHASE_TEST_DATASETS_ROOT"
ENV WPHASE_GREENS_FUNCTIONS="$ARG_WPHASE_GREENS_FUNCTIONS"

# Allow seiscomp to run as root:
RUN sed -i 's/asRoot = False/asRoot = True/' /opt/seiscomp/bin/seiscomp-control.py

# Configure environment for seiscomp:
ENV SEISCOMP_ROOT="/opt/seiscomp" \
    PATH="/opt/seiscomp/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
    LD_LIBRARY_PATH="/opt/seiscomp/lib" \
    PYTHONPATH="/opt/seiscomp/lib/python" \
    MANPATH="/opt/seiscomp/share/man"

# Configure seiscomp completion in bash:
RUN echo "source /opt/seiscomp/share/shell-completion/seiscomp.bash" >> /root/.bashrc

# Let seiscomp utilities run without requiring database dependencies:
RUN echo 'core.plugins = ""' > /opt/seiscomp/etc/global.cfg

# create required folders for wphase
RUN install -g "$(id -gn)" -o "$(id -un)" -d \
    "$WPHASE_GREENS_FUNCTIONS_DIR" \
    "$WPHASE_WEB_OUTPUTS_ROOT"

# This script is used by both dev and prod to set up for wphase
COPY docker/configure_wphase_image.sh "$LOCAL_INSTALL_FOLDER"/configure_wphase_image.sh

#############################################################################
# Dev image with auto-installation, to be used with run-or-build-container.sh
FROM base AS dev
ARG LOCAL_INSTALL_FOLDER=/tmp
RUN "$LOCAL_INSTALL_FOLDER"/configure_wphase_image.sh dev

##############################
# Fully-baked production image
FROM base AS prod
ARG LOCAL_INSTALL_FOLDER=/tmp
ARG WPHASE_HOME=/wphase/app
ENV WPHASE_HOME="$WPHASE_HOME"
COPY wphase "$WPHASE_HOME"/wphase
COPY scripts "$WPHASE_HOME"/scripts
COPY tests "$WPHASE_HOME"/tests
COPY setup.py "$WPHASE_HOME"/setup.py
RUN "$LOCAL_INSTALL_FOLDER"/configure_wphase_image.sh prod
ENTRYPOINT ["wphase"]
