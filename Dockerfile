FROM ubuntu:24.04 AS seiscomp
ADD https://github.com/GeoscienceAustralia/ga-wphase/releases/download/v0.3.2/seiscomp-6.5.1-ubuntu24.04-x86_64.tar.gz /opt/
RUN cd /opt && tar xfz seiscomp*.tar.gz

FROM ubuntu:24.04 AS wphase
COPY --from=seiscomp /opt/seiscomp/ /opt/seiscomp/
RUN apt-get update &&  \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
        build-essential \
        ca-certificates \
        gfortran \
        libboost-filesystem1.83.0 \
        libboost-iostreams1.83.0 \
        libboost-program-options1.83.0 \
        libboost-regex1.83.0 \
        libboost-system1.83.0 \
        libboost-thread1.83.0 \
        libmariadb3 \
        libmysqlclient21 \
        libncurses6 \
        libpq5 \
        libproj-dev \
        libssl3 \
        libxml2 \
        python-is-python3 \
        python3-dev \
    && apt-get clean autoclean -y && rm -rf /var/lib/apt/lists/*

COPY --from=ghcr.io/astral-sh/uv:0.5.6 /uv /uvx /bin/

# Set up non-root user
ARG DOCKER_USER_GID=1000
ARG DOCKER_USER_UID=1000

RUN userdel -r ubuntu && \
    groupadd -r -g $DOCKER_USER_GID wphase && \
    useradd -r -g wphase -u $DOCKER_USER_UID -s /bin/bash -d /home/wphase -c "docker user" wphase && \
    mkdir -p /home/wphase/app /opt/seiscomp && \
    chown -R wphase:wphase /home/wphase /opt/seiscomp

USER wphase
WORKDIR /home/wphase/app

ENV UV_LINK_MODE=copy

# Install pinned dependencies from PyPI
COPY --chown=wphase:wphase pyproject.toml uv.lock ./
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --no-install-project --locked

# Configure environment for seiscomp:
ENV SEISCOMP_ROOT="/opt/seiscomp" \
    PATH="/home/wphase/.local/bin:/opt/seiscomp/bin:${PATH}" \
    LD_LIBRARY_PATH="/opt/seiscomp/lib" \
    PYTHONPATH="/opt/seiscomp/lib/python"

COPY --chown=wphase:wphase wphase ./wphase
COPY --chown=wphase:wphase tests ./tests
COPY --chown=wphase:wphase README.md CMakeLists.txt ./

# Install wphase package
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked

# Default mountpoints:
ENV WPHASE_OUTPUT_DIR=/outputs
ENV WPHASE_GREENS_FUNCTIONS=/greens
ENTRYPOINT ["uv", "run", "wphase", "--console=1"]
