# Fedora 30 image with gdb-heap and cctbx for debugging heap memory

FROM fedora:30

# Install basic utilities, devel packages, debuginfo packages, and gdb-heap
RUN dnf install -y make which dnf-plugins-core && \
  dnf install -y python3-devel glibc-devel glib2-devel openssl-devel && \
  dnf debuginfo-install -y python3 glibc ncurses-libs readline openssl-libs zlib glib2 && \
  dnf install -y gdb-heap

# Pull in some gdb-heap changes on personal fork
RUN dnf install -y git && \
  cd /usr/share/gdb-heap/ && \
  rm heap/*.py && \
  git init && \
  git remote add origin https://github.com/dwpaley/gdb-heap && \
  git fetch origin && \
  git checkout origin/master


# CCTBX system dependencies and directory structure
RUN dnf install -y mesa-libGL-devel
RUN mkdir /img
RUN mkdir /mc3

# Conda
RUN cd /mc3 && \
  curl -LO https://repo.anaconda.com/miniconda/Miniconda3-py37_4.10.3-Linux-x86_64.sh && \
  /bin/sh Miniconda3-*.sh -b -u -p /mc3

# Bootstrap.py and conda yml file
RUN cd /img && \
  curl -LO https://raw.githubusercontent.com/cctbx/cctbx_project/master/xfel/conda_envs/psana_environment.yml && \
  curl -LO  https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py

# Set up conda_base environment using mamba
RUN cd /img && \
  source /mc3/etc/profile.d/conda.sh && \
  conda activate base && \
  conda install mamba -c conda-forge -y && \
  mamba env create -f psana_environment.yml -p $PWD/conda_base && \
  conda install -y -c conda-forge --no-deps -p $PWD/conda_base gemmi

# Get sources and build
RUN cd /img && \
  source /mc3/etc/profile.d/conda.sh && \
  conda activate $PWD/conda_base && \
  python bootstrap.py \
    --builder=xfel \
    --use_conda=$PWD/conda_base \
    --nproc=12 \
    --python=37 \
    --build-dir=build_debug \
    --config-flags="--build=debug" \
    hot update build

# Modify conda installation to use system python with debug symbols
RUN ln -sf /usr/bin/python3 /img/conda_base/bin/python3.7
