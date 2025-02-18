# syntax=docker/dockerfile:1.4
FROM debian:sid-slim

ENV RUNLEVEL 1
ENV DEBIAN_FRONTEND noninteractive
ENV LC_ALL C.UTF-8

RUN <<EOF
#!/usr/bin/env bash
set -exo pipefail

# update OS
apt-get -y update
apt-get -y upgrade

# install cctbx
declare -ar pkgs=(
    # this is huge, qt, hdf...
    'python3-cctbx'
    'g++'
)
apt-get -y --no-install-recommends install "${pkgs[@]}"

# run fable to check works
fable.cout --example

# docker cleanup
apt-get -y autoremove --purge
rm -rf /var/lib/apt/lists/*
EOF

# /usr/bin
ENTRYPOINT ["fable.cout"]
