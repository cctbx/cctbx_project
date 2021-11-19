# Debugging and memory analysis for CCTBX applications

Daniel Paley, 2021  
dwpaley@lbl.gov

The CCTBX `bootstrap.py` installer provides support for building with debug symbols. This permits
interactive debugging in `gdb` or similar.  Heap memory analysis using the gdb extension
[gdb-heap](https://github.com/rogerhu/gdb-heap) is also possible.

## Building

Build as usual via `bootstrap.py` but use the additional flag `--config-flags="--build=debug`. You
can also set the build directory using `--build-dir=<path>`. To create a second debug build
alongside an existing one, this could look like:
```
$ python bootstrap.py \
   --builder=xfel \
   --use_conda=$PWD/conda_base \
   --nproc=16 --build-dir=build_debug \
   --config-flags="--build=debug" \
   build
```

For a new build on a Centos 7.4 server I do the following, which should be adapted for your
environment:
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

I have not had a lot of success running `$ gdb <command>`, which I assume is connected to the
environment setup that the cctbx dispatchers make before starting Python. Instead, we use Python
breakpoints and then connect gdb to the running process.

Start two shells on the same machine. Modify your Python program with this line where you want to
begin:
```
import os; print(os.getpid()); import pdb; pdb.set_trace()
```
In the other terminal, start gdb. When you have the PID, do `attach <PID>`. Now you have two
debuggers running simultaneously.  In general you will want to control execution in gdb using
`<ctrl-c>` and `continue`. The Python session won't accept any input until you `continue` in gdb.

## Memory analysis

The gdb extension `gdb-heap` is a powerful tool for examining stack memory. An example is given
here in which we diagnose a problem with excessive memory usage in the `dxtbx::ImageSet` class.

### Configuration

The configuration is a bit of an ordeal; running in a Docker container seems to work well. A
prebuilt image on Fedora 30 with a CCTBX build and the necessary debug symbols for Python, glibc,
etc. is available from Dockerhub as `dwpaley/cctbx-gdbheap`. The corresponding Dockerfile is also
provided in this directory. A few modifications were needed to make gdb-heap compatible with
Python3 and are available from https://github.com/dwpaley/gdb-heap. Also note that a python3.7 with
debug symbols provided by Fedora was symlinked in the cctbx installation at
`/img/conda_base/bin/python3.7`.

### Example

It was observed that construction of dxtbx.model.ExperimentList objects seemed to consume excessive
memory for some serial diffraction applications. The problem was observed for single-image
experiment lists with the underlying data coming from a multi-image container file. The file
`split_00.expt` in this directory is a good example; by doing the following in `libtbx.python` we
can see that about 1.5 MB is used by simply constructing this single-experiment ExperimentList.

```
>>> from dxtbx.model import ExperimentList
>>> el1 = ExperimentList.from_file('split_00.expt') # After this step, `top` shows the process at 245,884 kB memory
>>> el2 = ExperimentList.from_file('split_00.expt') # Now we are using 247,204 kB
```

We will load the experiment using our Docker image and explore the heap.  The example expt file is
in `cctbx_project/dox/rst/data`, which I will assume is the working directory here. The data is
exposed via a bind mount: `-v $PWD:/data:ro` mounts the working directory as a read-only filesystem
at `/data` in the container. After starting the Python interpreter, we do a throwaway load of the
experiment file to ensure any needed imports have been done.

```
$ docker clone dwpaley/cctbx-gdbheap
$ docker run -it --cap-add=SYS_PTRACE -v $PWD:/data:ro dwpaley/cctbx-gdbheap

# source /img/build_debug/conda_setpaths.sh
# libtbx.python

>>> from dxtbx.model import ExperimentList as el
>>> el1 = el.from_file('/data/split_00.expt', check_format=False)
>>> import os; os.getpid()
400
```
Open a second shell in the running container. For brevity, from this point forward, some terminal
input will be labeled with a Python `>>>` prompt and some will have a `(gdb)` prompt. These inputs
will be entered respectively in the first and second shells we started in the container.

Once you have started the new shell, attach gdb to the running Python process and save the current
state of the heap:
```
$ docker ps
$ docker exec -it <container id> /bin/bash

# gdb -p 400
(gdb) heap label before
Blocks retrieved 10000
Blocks retrieved 20000
Blocks retrieved 30000
24,280,536 allocated, in 35636 blocks
(gdb) c
```

Now we return to our Python session and load a second experiment:
```
>>> el2 = el.from_file('/data/split_00.expt', check_format=False)
```

And see what has changed on the heap:
```
(gdb) <ctrl-c>
Program received signal SIGINT, Interrupt.
(gdb) heap diff before
Blocks retrieved 10000
Blocks retrieved 20000
Blocks retrieved 30000
Changes from e to current
   +1,452,176 bytes, +44 blocks

  Free-d blocks:
  New blocks:
    0x00005600d3fa3410 -> 0x00005600d3fa344f       64 bytes uncategorized::64 bytes |@....V.......V......|
    0x00005600d3fa3450 -> 0x00005600d3fa348f       64 bytes uncategorized::64 bytes |.4...V.......V......|
[...]
    0x00005600d54a3320 -> 0x00005600d54a336f       80 bytes uncategorized::80 bytes |.....V.......V..le('|
    0x00005600d54a3370 -> 0x00005600d54a339f       48 bytes uncategorized::48 bytes |..(................|
    0x00005600d54a36f0 -> 0x00005600d54fa54f   355936 bytes uncategorized::355936 bytes |....................|
    0x00005600d54fa550 -> 0x00005600d55513af   355936 bytes uncategorized::355936 bytes |....................|
    0x00005600d55513b0 -> 0x00005600d55a820f   355936 bytes uncategorized::355936 bytes |....................|
    0x00005600d55a8210 -> 0x00005600d55ff06f   355936 bytes uncategorized::355936 bytes |....................|
```
The last four chunks on the heap are each 356 kB, accounting for the vast majority of the 1.45 MB
used by the experiment.  I don't know how to directly classify the objects stored there, but we
will try the following trick: (1) set memory watchpoints on the addresses; (2) delete the Python
object we just created; (3) see what object we are destructing when we touch that memory.

Step 1: Set watchpoints. Only four watchpoints are supported at once, so we use the start and end
of the first two chunks.
```
(gdb) awatch *0x00005600d54a36f0
(gdb) awatch *0x00005600d54fa54f
(gdb) awatch *0x00005600d54fa550
(gdb) awatch *0x00005600d55513af
(gdb) c
```
Step 2: Delete the Python object.
```
>>> del el2
```
Step 3: We hit a watchpoint and check the stack trace.

```
Hardware access (read/write) watchpoint 20: *0x00005600d54fa54f

Value = 0
0x00007fbd1daa5f84 in __GI___libc_free (mem=0x5600d54fa550) at malloc.c:3109
3109	  p = mem2chunk (mem);
(gdb) bt
#0  0x00007fbd1daa5f84 in __GI___libc_free (mem=0x5600d54fa550) at malloc.c:3109
#1  0x00007fbd1b4dbde1 in scitbx::af::aligned_free (p=0x5600d54fa550) at /img/modules/cctbx_project/scitbx/array_family/memory.h:71
#2  0x00007fbd1b4dbe8f in scitbx::af::sharing_handle::~sharing_handle (this=0x5600d51d8c80, __in_chrg=<optimized out>)
    at /img/modules/cctbx_project/scitbx/array_family/shared_plain.h:77
#3  0x00007fbd1b4dbeaa in scitbx::af::sharing_handle::~sharing_handle (this=0x5600d51d8c80, __in_chrg=<optimized out>)
    at /img/modules/cctbx_project/scitbx/array_family/shared_plain.h:82
#4  0x00007fbd13e83631 in scitbx::af::shared_plain<boost::shared_ptr<dxtbx::model::Detector> >::m_dispose (this=0x5600d51d8210)
    at /img/modules/cctbx_project/scitbx/array_family/shared_plain.h:374
#5  0x00007fbd13e733c4 in scitbx::af::shared_plain<boost::shared_ptr<dxtbx::model::Detector> >::~shared_plain (this=0x5600d51d8210, __in_chrg=<optimized out>)
    at /img/modules/cctbx_project/scitbx/array_family/shared_plain.h:227
#6  0x00007fbd13e610fa in scitbx::af::shared<boost::shared_ptr<dxtbx::model::Detector> >::~shared (this=0x5600d51d8210, __in_chrg=<optimized out>)
    at /img/modules/cctbx_project/scitbx/array_family/shared.h:11
#7  0x00007fbd13e63bfc in dxtbx::ImageSetData::~ImageSetData (this=0x5600d51d81e8, __in_chrg=<optimized out>) at /img/modules/dxtbx/src/dxtbx/imageset.h:154
#8  0x00007fbd13e63cc8 in dxtbx::ImageSet::~ImageSet (this=0x5600d51d81e0, __in_chrg=<optimized out>) at /img/modules/dxtbx/src/dxtbx/imageset.h:546
#9  0x00007fbd13e63ce4 in dxtbx::ImageSet::~ImageSet (this=0x5600d51d81e0, __in_chrg=<optimized out>) at /img/modules/dxtbx/src/dxtbx/imageset.h:546
#10 0x00007fbd13e9d4a9 in boost::checked_delete<dxtbx::ImageSet> (x=0x5600d51d81e0) at /img/modules/boost/boost/core/checked_delete.hpp:36
#11 0x00007fbd13f387f6 in boost::detail::sp_counted_impl_p<dxtbx::ImageSet>::dispose (this=0x5600d51de960)
    at /img/modules/boost/boost/smart_ptr/detail/sp_counted_impl.hpp:93
[...]
```
Scanning down the stack we see the memory is accessed while destructing a `dxtbx::ImageSetData`
(frame 7). In that process we are destructing a `scitbx::af::shared` array of `boost::shared_ptr`
to `dxtbx::model::Detector`. This must be member data of the `ImageSetData`, so we look at the
bottom of the class definition in `dxtbx/src/dxtbx/imageset.h`:
```
    499   scitbx::af::shared<beam_ptr> beams_;
    500   scitbx::af::shared<detector_ptr> detectors_;
    501   scitbx::af::shared<goniometer_ptr> goniometers_;
    502   scitbx::af::shared<scan_ptr> scans_;
```
where `beam_ptr` et al. have been typedef'd to the boost shared pointers we saw in the stack trace:
```
    156   typedef boost::shared_ptr<BeamBase> beam_ptr;
[etc.]   
```
We move up the stack to the ImageSetData class and confirm that the arrays `beams_`, etc. are the
large objects we saw:
```
(gdb) frame 7
#7  0x00007fbd13e63bfc in dxtbx::ImageSetData::~ImageSetData (this=0x5600d51d81e8, __in_chrg=<optimized out>) at /img/modules/dxtbx/src/dxtbx/imageset.h:154
154	class ImageSetData {
(gdb) p beams_.size()
$1 = 22245
```
We look at the `ImageSetData` constructor to see why the arrays are so big:
```
    169   ImageSetData(boost::python::object reader, masker_ptr masker)
    170       : reader_(reader),
    171         masker_(masker),
    172         beams_(boost::python::len(reader)),
    173         detectors_(boost::python::len(reader)),
    174         goniometers_(boost::python::len(reader)),
    175         scans_(boost::python::len(reader)),
    176         reject_(boost::python::len(reader)) {}
```
where the `reader` argument is the Python class `FormatMultiImage.Reader` with these methods
`__len__` and `__init__`:
```
     26     def __init__(self, format_class, filenames, num_images=None, **kwargs):
     [...]
     31         if num_images is None:
     32             format_instance = self.format_class.get_instance(
     33                 self._filename, **self.kwargs
     34             )
     35             self._num_images = format_instance.get_num_images()
     36         else:
     37             self._num_images = num_images
     [...]
     50     def __len__(self):
     51         return self._num_images
```
By separately setting a trace in the `Reader` constructor, we find the argument `num_images` was
set as follows in `FormatMultiImage.get_imageset`:
```
    148         if cls.get_num_images == FormatMultiImage.get_num_images:
    149             assert single_file_indices is not None
    150             assert min(single_file_indices) >= 0
    151             num_images = max(single_file_indices) + 1
    [...]
    175    	    reader = cls.get_reader()(filenames, num_images=num_images, **format_kwargs)
```
Finally we confirm that `single_file_indices` is a large number for this experiment file:
```
# grep -A2 single_file_indices /data/split_00.expt
      "single_file_indices": [
        22244
      ],
```
We note that 22,244 is not an unusual number of images in a single container file; for example it
corresponds to 3 minutes of data at 120 Hz.

The arrays `beams_`, `detectors_`, etc. each contain 22,245 `shared_ptr`s, which are each 16 B
(i.e. two 64-bit addresses). Therefore the memory footprint of 4x356 kB is fully explained by this
problem.

[DXTBX PR 438](https://github.com/cctbx/dxtbx/pull/438) reduces this memory usage by making the
`ImageSetData` member arrays len(indices) instead of len(max(indices)).
