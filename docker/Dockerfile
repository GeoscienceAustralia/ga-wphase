################################################
# Base image with all dependencies and seiscomp3
FROM centos:centos7 AS base

# Seiscomp3 dependencies:
RUN yum install -y \
    boost \
    flex \
    libxml2 \
    realpath \
    && yum clean all

# Install seiscomp3:
ARG SEISCOMP_TARBALL_URL=https://www.seiscomp.de/seiscomp3/downloader/seiscomp3-jakarta-2020.330.02-centos7-x86_64.tar.gz
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

COPY docker/$BASEMAP_ARCHIVE_NAME "$LOCAL_INSTALL_FOLDER"/"$BASEMAP_ARCHIVE_NAME"
COPY docker/requirements.txt /requirements.txt

# The installation scripts for dependencies are separated here because they
# take a while and don't change very often, so it's easiest to keep them
# separate for quick development.
COPY docker/install_dependencies.sh "$LOCAL_INSTALL_FOLDER"/install_dependencies.sh
RUN "$LOCAL_INSTALL_FOLDER"/install_dependencies.sh

# Allow seiscomp to run as root:
RUN sed -i 's/asRoot = False/asRoot = True/' /opt/seiscomp3/bin/seiscomp-control.py

# Configure environment for seiscomp:
ENV SEISCOMP_ROOT="/opt/seiscomp3" \
    PATH="/opt/seiscomp3/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
    LD_LIBRARY_PATH="/opt/seiscomp3/lib" \
    PYTHONPATH="/opt/seiscomp3/lib/python" \
    MANPATH="/opt/seiscomp3/share/man"

# Configure seiscomp completion in bash:
RUN echo "source /opt/seiscomp3/share/shell-completion/seiscomp.bash" >> /root/.bashrc

# Let seiscomp utilities run without requiring database dependencies:
RUN echo 'core.plugins = ""' > /opt/seiscomp3/etc/global.cfg

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
#ARG WPHASE_HOME=/wphase/app
#ENV WPHASE_HOME="$WPHASE_HOME"
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
