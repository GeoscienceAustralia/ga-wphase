FROM centos:centos7

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

COPY $BASEMAP_ARCHIVE_NAME "$LOCAL_INSTALL_FOLDER"/"$BASEMAP_ARCHIVE_NAME"
COPY requirements.txt /requirements.txt

# The installation scripts are separated into two layers because
# install_depedencies is very slow and doesn't change very often.
COPY install_dependencies.sh "$LOCAL_INSTALL_FOLDER"/install_dependencies.sh
RUN "$LOCAL_INSTALL_FOLDER"/install_dependencies.sh

COPY configure_wphase_image.sh "$LOCAL_INSTALL_FOLDER"/configure_wphase_image.sh
RUN "$LOCAL_INSTALL_FOLDER"/configure_wphase_image.sh
