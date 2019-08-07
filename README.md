# Computational Crystallography Toolbox

Greetings, earthling. We see you intend to embark on the arduous quest of becoming an elite computational crystallographer. Many before you have dared to try. We hope the CCTBX package grants you success along your journey. 

## Contents

* [About](#intro)

* [Installation (Mac OS / Linux)](#install)

  * [Download the Bootstrap installer script](#getboot)
     * [hot / update](#hot)
     * [base](#conda)
     * [build (auto)](#build)
     * [build (manual)](#manual_build)
  * [Bootstrap usage (full)](#boot_doc)
  * [Install without conda](#deprecated)

* [Installation on Windows](#windows)

* [Actually using CCTBX](#using)

* [Examples](#Examples)

* [Additional information] (#additional)
 

<a name="intro"></a>
## About 

The Computational Crystallography Toolbox (cctbx) is being developed as the open source component of the PHENIX system. The goal of the PHENIX project is to advance automation of macromolecular structure determination. PHENIX depends on the cctbx, but not vice versa. This hierarchical approach enforces a clean design as a reusable library. The cctbx is therefore also useful for small-molecule crystallography and even general scientific applications.

The cctbx also provides some of the key component of the Olex 2 software. Olex 2 is dedicated to the workflow of small molecule crystallographic studies. It features a powerful and flexible refinement engine, olex2.refine, which is developed as part of the cctbx,
in the smtbx top-module.

To maximize reusability and, maybe even more importantly, to give individual developers a notion of privacy, the cctbx is organized as a set of smaller modules. This is very much like a village (the cctbx project) with individual houses (modules) for each family (groups of developers, of any size including one).

The cctbx code base is available without restrictions and free of charge to all interested developers, both academic and commercial. The entire community is invited to actively participate in the development of the code base. A sophisticated technical infrastructure that enables community based software development is provided by GitHub. This service is also free of charge and open to the entire world.

The cctbx is designed with an open and flexible architecture to promote extendability and easy incorporation into other software environments. The package is organized as a set of ISO C++ classes with Python bindings. This organization combines the computational efficiency of a strongly typed compiled language with the convenience and flexibility of a dynamically typed scripting language in a strikingly uniform and very maintainable way.

Use of the Python interfaces is highly recommended, but optional. The cctbx can also be used purely as a C++ class library.


<a name="install"></a>
## Installation (Mac OS / Linux)

Current efforts are to incorporate cctbx into a conda environment. Here is a basic installation workflow. Note, cctbx is closely linked to the DIALS diffraction software package, and the DXTBX software package, which provide algorithms for processing experimental diffraction data. Therefore this build will include builds of DIALS and DXTBX, as well as a few others. 


<a name="getboot"></a>
### Download the installer script

The installer script [bootstrap](https://raw.githubusercontent.com/cctbx/cctbx_project/darwin_omp/libtbx/auto_build/bootstrap) downloads dependencies and builds the software. Place it in your main working directory, e.g. ```~/crystal```

```
mkdir ~/crystal  # make a working directory
cd ~/crystal
wget https://raw.githubusercontent.com/cctbx/cctbx_project/darwin_omp/libtbx/auto_build/bootstrap
``` 

Make ```bootstrap.py``` an executable and look at the options (we'll just use a few)

```
chmod +x bootstrap.py
./bootstrap.py --help
```

Note, in the past bootstrap required that the first python binary in your path be python version 2, or else it will fail. Verify this by running
 
```
python -c "import sys;print (sys.version_info[0])"
```

and seeing a 2 printed to the screen. If its a 3, you might need to change your PATH variable and visible python binaries accordingly.

<a name="hot"></a>
### Hot / update: Getting cctbx project sources 
Typically the first step is to download the internal sources created by the CCTBX project developers and collaborators. This is done automatically using the ```hot``` and ```update``` arguments to bootstrap: 
 
```
./bootstrap.py --builder=dials hot update
```

Running Hot and update with the ```--builder=dials``` argument will download the packages that dials depends on, in this case ```cctbx_project``` and all its goodies. This is the *"builder"* that most developers will use and all of the source materials are available to the public. The packages will be stored in the newly created ```~/crystal/modules``` folder.

```./bootstrap.py hot update --builder=dials``` can be run multiple times to bring in the latests updates to the sources in the ```modules``` folder


<a name="conda"></a>
### Base: getting external dependencies
Now you want to download a base python installation as well as any other external dependencies. This is most easily done using the ```conda``` package manager, and bootstrap makes this relatively painless. 

Bringing in a new conda uses about 2 GB of space. If however you already have a conda install that you love, then with a little help you can use that one as well. See towards the end of this section. 

**To bring a new conda install in:** first verify there is not conda in your path by typing ```conda``` in the terminal and verifying its not there. Next, verify you do not have a ```CONDA_PREFIX``` environment variable set (sometimes this can be set in a .bashrc for example, or sourced from somewhere not so obvious). To do this , type ```printenv | grep CONDA``` into the terminal and verify the output is empty. Also, verify you don't have ```miniconda``` folders in home folders or you don't have a ```~/.condarc``` file. If you have them already, you might want to proceed to the below instructions for installing with a preexisting conda However, to proceed with a fresh install, type 

```
./bootstrap.py base --use_conda --builder=dials
```

and it will download a ```mc3``` folder (miniconda3) and create a ```conda_base``` folder in the current directory. The ```conda_base``` is actually a conda environment, and it has all of the cctbx dependencies for your current operating system.

**If you already have conda installed**: Create an empty conda environment

```
source ~/miniconda3/etc/profile.d/conda.sh
conda create -name cctbx
./bootstrap.py base --use_conda ~/miniconda3/envs/cctbx --builder=dials
```

This will also bring in a lot of downloads.


<a name="build"></a>
### Build: Auto-building

One can auto-build with the bootstrap script. Typically, after the conda environment is setup as described above, one can

```
./bootstrap.py build --use_conda ./conda_base --build-dir ompbuild --nproc 8
```

This should compile the code in a build folder with the default configuration for your current OS. It is useful to have multiple build folders, for example a CUDA-compatible build. In such a case one needs the NVCC configuration flags, so one would

```
./bootstrap.py build --use_conda ./conda_base --build-dir cudabuild --config-flags="--enable_cuda" --nproc 8
```

In order to use run cctbx scripts, one should always source the setpaths.sh script in the desired build folder.

The full list of configuration flags is:

```
Usage: libtbx/configure.py [options] module_name[=redirection_path] ...

Options:
  -h, --help            show this help message and exit
  -r DIRECTORY, --repository=DIRECTORY
                        path to source code repository (may be specified
                        multiple times; paths are searched in the order given)
  --current_working_directory=DIRECTORY
                        preferred spelling of current working directory (to
                        resolve ambiguities due to soft links)
  --build=release|max_optimized|quick|debug|debug_optimized|profile
                        build mode (default: release)
  --compiler=STRING     select non-standard compiler (platform dependent)
  --warning_level=WARNING_LEVEL
                        manipulate warning options (platform dependent)
  --static_libraries    build all libraries statically
  --static_exe          link all executables statically (implies
                        --static_libraries)
  --scan_boost          enable implicit dependency scan
  --write_full_flex_fwd_h
                        create full flex_fwd.h files to work around platform-
                        specific problems (see comments in
                        cctbx.source_generators.flex_fwd_h)
  --command_version_suffix=STRING
                        version suffix for commands in bin directory
  --use_environment_flags
                        add compiler flags from environment variables:
                        CXXFLAGS, CFLAGS, CPPFLAGS, LDFLAGS
  --force_32bit         Force 32-bit compilation on Mac OS 10.6 (Snow Leopard)
                        Not compatible with /usr/bin/python: please run
                        configure with /System/Library/Frameworks/Python.frame
                        work/Versions/2.x/bin/python
  --msvc_arch_flag=None|SSE|SSE2
                        choose MSVC CPU architecture instruction set for
                        optimized builds
  --build_boost_python_extensions=True|False
                        build Boost.Python extension modules (default: True)
  --enable_openmp_if_possible=True|False
                        use OpenMP if available and known to work (default:
                        False)
  --enable_boost_threads=ENABLE_BOOST_THREADS
                        make Boost.Threads available
  --enable_cuda         Use optimized CUDA routines for certain calculations.
                        Requires at least one NVIDIA GPU with compute
                        capability of 2.0 or higher, and CUDA Toolkit 4.0 or
                        higher (default: False)
  --use_conda           Use conda as the source for Python and dependencies
                        (default: False)
  --opt_resources=True|False
                        use opt_resources if available (default: False)
  --precompile_headers  Precompile headers, especially Boost Python ones
                        (default: don't)
  --boost_python_no_py_signatures
                        disable Boost.Python docstring Python signatures
  --boost_python_bool_int_strict
                        disable Boost.Python implicit bool<->int conversions
  --clear_scons_memory  remove scons build signatures and config cache
  --no_bin_python       do not create <build directory>/bin/python even in a
                        development environment
  --enable_cxx11        use C++11 standard
  --python3warn=none|warn|fail
                        Python3 migration warnings. 'warn': print warnings
                        when running code that may cause problems in Python 3.
                        'fail': stop execution on warnings. 'none': disable
                        warnings (default)
  --skip_phenix_dispatchers
                        Skip all dispatchers with 'phenix' in the title
```

<a name="boot_doc"></a>
### Bootstrap usage

The full list of bootstrap options is

```
usage: bootstrap.py [-h] [--builder BUILDER] [--cciuser CCIUSER]
                    [--sfuser SFUSER] [--revert REVERT] [--sfmethod SFMETHOD]
                    [--git-ssh] [--git-reference GIT_REFERENCE]
                    [--with-python WITH_PYTHON] [--nproc NPROC]
                    [--download-only] [-v] [--skip-base-packages SKIP_BASE]
                    [--force-base-build] [--enable-shared] [--mpi-build]
                    [--python3] [--wxpython4] [--config-flags CONFIG_FLAGS]
                    [--use-conda [ENV_DIRECTORY]] [--build-dir BUILD_DIR]
                    [action [action ...]]

  You may specify one or more actions:
    hot - Update static sources (scons, etc.)
    update - Update source repositories (cctbx, cbflib, etc.)
    base - Build base dependencies (python, hdf5, wxWidgets, etc.)
    build - Build
    tests - Run tests
    doc - Build documentation

  The default action is to run: hot, update, base, build

  You can specify which package will be downloaded, configured,
  and built with "--builder". Current builders:
    cctbx, cctbxlite, dials, external, labelit, molprobity, phaser, phaser_tng,
    phenix, phenix_tng, qrefine, xfel

  You can provide your SourceForge username with "--sfuser", and
  your CCI SVN username with "--cciuser". These will checkout
  and update repositories with your credentials. Some builders,
  like phenix, require this argument for access to certain
  repositories.

  You can run the compilation step in parallel by providing a
  the number of processes using "--nproc".
  Complete build output is shown with "-v" or "--verbose".

  Finally, you may specify a specific Python interpreter
  using "--with-python".

  Example:

    python bootstrap.py --builder=cctbx --sfuser=metalheadd hot update build tests
  

positional arguments:
  action                Actions for building

optional arguments:
  -h, --help            show this help message and exit
  --builder BUILDER     Builder: phenix_tng,molprobity,dials,xfel,external,pha
                        ser,phenix,cctbx,phaser_tng,qrefine,labelit,cctbxlite
  --cciuser CCIUSER     CCI SVN username.
  --sfuser SFUSER       SourceForge SVN username.
  --revert REVERT       SVN string to revert all SVN trees
  --sfmethod SFMETHOD   SourceForge SVN checkout method.
  --git-ssh             Use ssh connections for git. This allows you to commit
                        changes without changing remotes and use reference
                        repositories.
  --git-reference GIT_REFERENCE
                        Path to a directory containing reference copies of
                        repositories for faster checkouts.
  --with-python WITH_PYTHON
                        Use specified Python interpreter
  --nproc NPROC         number of parallel processes in compile step.
  --download-only       Do not build, only download prerequisites
  -v, --verbose         Verbose output
  --skip-base-packages SKIP_BASE
  --force-base-build
  --enable-shared
  --mpi-build           Builds software with mpi functionality
  --python3             Install a Python3 interpreter. This is unsupported and
                        purely for development purposes.
  --wxpython4           Install wxpython4 instead of wxpython3. This is
                        unsupported and purely for development purposes.
  --config-flags CONFIG_FLAGS, --config_flags CONFIG_FLAGS
                        Pass flags to the configuration step. Flags should be
                        passed separately with quotes to avoid confusion (e.g
                        --config_flags="--build=debug" --config_flags="--
                        enable_cxx11")
  --use-conda [ENV_DIRECTORY], --use_conda [ENV_DIRECTORY]
                        Use conda for dependencies. The directory to an
                        existing conda environment can be provided. The build
                        will use that environment instead of creating a
                        default one for the builder. Also, if a conda
                        environment is currently active, $CONDA_PREFIX will be
                        used if ENV_DIRECTORY is not provided. Specifying an
                        environment or using $CONDA_PREFIX is for developers
                        that maintain their own conda environment.
  --build-dir BUILD_DIR
                        directory where the build will be. Should be at the
                        same level as modules! default is 'build'
```

<a name="manual_build"></a>
### Manual building
This method of building is not as stable so it is highly recommended to use the auto build. The most straightforward way to manually build is to  create a build folder alongside the ```modules``` and ```conda_base``` folders. 

To start, you will want to activate the newly created conda environment:

```
source ./miniconda3/etc/profile.d/conda.sh
conda activate ./conda_base
mkdir build
cd build
```

At this point, it doesnt hurt to verify your python is as you expect it to be: ```which python``` should point to the ```conda_base/bin/python```. 

Now, from inside the build folder, we are going to run the configure script which is in the cctbx_project. This sets up your system for a build step

```
python ../modules/cctbx_project/libtbx/configure.py dials --enable_openmp_if_possible=True --use_conda
```

**I am not sure what happens if you configure once with some options, but pull in a new repository that requires different options**

After configuring, run make from within the build directory, twice, to compile and to link the binaries. 

```
make
make
```

<a name="deprecated"></a>
### Deprecated installation instructions
One can just execute bootstrap and let the default settings work the magic of installation. While seemingly more simple, **do this at your own risk.** Current efforts to use conda have made cross platform builds and post-install development much more stable. Note this installation procedure will take longer, and requires a lot of space:

```
./bootstrap.py 
```

<a name="windows"></a>
## Installing on windows
On Windows follow the instructions detailed [here](https://github.com/cctbx/cctbx_project/wiki/How-to-build-CCTBX-on-Windows)

<a name="using"></a>
## Actually using cctbx
Once you have a successful build, one can set the environment like so

```
source ~/crystal/build/setpaths.sh
```

which brings in all of the relevant binaries. 

One can access the cctbx python base using the binary wrapper script ```libtbx.python```

```
libtbx.python -c "import sys;print (sys.prefix)"
```

which should return something like ```~/crystal/conda_base```.

Work can also be done interactively

```
libtbx.ipython
```


Further, after having built with the ```--builder=dials``` options above, one can now use the dials command line scripts to process or view data, e.g.

```
dials.image_viewer  my.cbf
```

To see more about DIALS checkout the [dials homepage](https://dials.github.io/).


<a name="Examples"></a>
## Examples
Many example scripts can be found at this [link](https://cctbx.github.io/). Run the scripts using

```libtbx.python my_example.py```

<a name="additional"></a>
## Additional information

#### Updating internal dependencies
hot update can be run at anytime to update the modules

#### Adding more python programs

Additional python dependencies can be pulled in with pip by running e.g.

```
libtbx.python -m pip install joblib
```

#### Adding new modules or other github repos
New python modules can be added by putting them in the modules folder

```
cd modules
git clone https://githiub.com/newproject.git
libtbx.configure newproject
```

#### Compiling code during development
If your new modules involve code that needs to be compiled, after configuring , go to the relevant build folder and run make 

```
# after updating C code in one of the modules
cd ~/crystal/build. # or whatever your build dir is named
make
```




