# Debugging and memory analysis for CCTBX applications

Daniel Paley, 2021  
dwpaley@lbl.gov

The CCTBX `bootstrap.py` installer provides support for building with debug symbols. This permits interactive debugging in `gdb` or similar.
Heap memory analysis using the gdb extension [gdb-heap](https://github.com/rogerhu/gdb-heap) is also possible.

## Building

Build as usual via `bootstrap.py` but use the additional flag `--config-flags="--build=debug`. You can also set the build directory using
`--build-dir=<path>`. To create a second debug build alongside an existing one, this could look like:
```
$ python bootstrap.py \
   --builder=xfel \
   --use_conda=$PWD/conda_base \
   --nproc=16 --build-dir=build_debug \
   --config-flags="--build=debug" \
   build
```

For a new build on a Centos 7.4 server I do the following, which should be adapted for your environment:
```
$ curl -LO https://raw.githubusercontent.com/cctbx/cctbx_project/master/xfel/conda_envs/psana_environment.yml
$ curl -LO  https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py
$ conda activate base
$ conda install mamba -c conda-forge -y
$ mamba env create -f psana_environment.yml -p $PWD/conda_base
$ conda activate $PWD/conda_base
$ python bootstrap.py \
   --builder=xfel \
   --use_conda=$PWD/conda_base \
   --nproc=60 \
   --python=37 \
   --build-dir=build_debug \
   --config-flags="--build=debug" \
   hot update build
```

## Debugging

I have not had a lot of success running `$ gdb <command>`, which I assume is connected to the environment setup that
the cctbx dispatchers make before starting Python. Instead, we use Python breakpoints and then connect gdb to the running
process.

Start two shells on the same machine. Modify your Python program with this line where you want to begin:
```
import os; print(os.getpid()); import pdb; pdb.set_trace()
```
In the other terminal, start gdb. When you have the PID, do `attach <PID>`. Now you have two debuggers running simultaneously.
In general you will want to control execution in gdb using `<ctrl-c>` and `continue`. The Python session won't accept any input
until you `continue` in gdb.

## Memory analysis

The gdb extension `gdb-heap` is a powerful tool for examining stack memory. 
