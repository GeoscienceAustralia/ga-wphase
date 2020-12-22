# Docker container

The Dockerfile in this directory uses a
[multistage](https://docs.docker.com/develop/develop-images/multistage-build/)
build:

- The stage `base` builds a docker container with the *dependencies* of
  w-phase installed, but not w-phase itself.
- The stage `dev` is a very thin layer on top of `base` that configures
  it to automatically build wphase and launch an interactive shell on each run.
  It is intended to be used for development purposes by mounting the w-phase
  source tree inside it.  See
  [run-or-build-container.sh](../run-or-build-container.sh) for how to do this.
- The stage `prod` installs wphase into the container and sets the entrypoint
  to run the wphase script. This is intended to be used for one-off invocations
  in production.
