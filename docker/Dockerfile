################################################
# Base image with all dependencies and seiscomp
FROM ubuntu:20.04

# SeisComP + W-Phase dependencies from apt:
RUN apt-get update &&  \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
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
        python-is-python3 \
        python3-dev \
        python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Set up non-root user
ARG DOCKER_USER_GID=1000
ARG DOCKER_USER_UID=1000

RUN groupadd -r -g $DOCKER_USER_GID wphase && \
    useradd -r -g wphase -u $DOCKER_USER_UID -s /bin/bash -d /home/wphase -c "docker user" wphase && \
    mkdir -p /home/wphase/app /opt/seiscomp && \
    chown -R wphase:wphase /home/wphase /opt/seiscomp

USER wphase
WORKDIR /home/wphase

# Install pinned dependencies from PyPI
COPY requirements.txt ./requirements.txt
RUN pip3 install --user --no-cache-dir --no-deps -r requirements.txt

# Install seiscomp 4.5:
ARG SEISCOMP_TARBALL_URL=https://www.seiscomp.de/downloader/seiscomp-4.5.0-ubuntu20.04-x86_64.tar.gz
RUN cd /opt/ && curl -s $SEISCOMP_TARBALL_URL | tar xfz -

# Configure environment for seiscomp:
ENV SEISCOMP_ROOT="/opt/seiscomp" \
    PATH="/home/wphase/.local/bin:/opt/seiscomp/bin:${PATH}" \
    LD_LIBRARY_PATH="/opt/seiscomp/lib" \
    PYTHONPATH="/opt/seiscomp/lib/python"
    # PYTHONPATH is tolerable here because we only have one python version installed in the container

WORKDIR /home/wphase/app
COPY --chown=wphase:wphase wphase ./wphase
COPY --chown=wphase:wphase scripts ./scripts
COPY --chown=wphase:wphase tests ./tests
COPY --chown=wphase:wphase setup.py ./setup.py

COPY --chown=wphase:wphase docker/configure_wphase_image.sh /tmp/configure_wphase_image.sh
RUN /tmp/configure_wphase_image.sh

ENTRYPOINT ["wphase", "--console=1"]
