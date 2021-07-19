# Conda Environments

These files specify a ```conda``` environment that contains the dependencies required for CCTBX projects. They are named according to ```<program>_<python version>_<platform>.txt```. So for the file, ```cctbx_py27_linux-64.txt```, the environment is for building the base CCTBX package using Python 2.7 on a 64-bit Linux distribution. The ```<platform>``` follows the ```conda``` pattern for naming platforms, so ```linux-64``` refers to 64-bit linux, ```osx-64``` refers to 64-bit macOS, and ```win-64``` refers to 64-bit Windows.
## Building with ```conda```
Assuming you have already downloaded the relevant source files with ```python bootstrap.py hot update --builder=<builder>``` into the ```<root directory>```, these directions let you build CCTBX packages with dependencies in a```conda``` environment. There is a ```--use-conda``` flag that will use these environments by default. The flag will also accept ```$CONDA_PREFIX``` as an argument if you want a custom environment for building.

1) Building with default environment,
```python bootstrap.py base build --builder=<builder> --use-conda```
2) Building with the default Python 3.8 environment,
```python bootstrap.py base build --builder=<builder> --use-conda --python=38```
3) Building with the currently active environment
```python bootstrap.py base build --builder=<builder> --use-conda=${CONDA_PREFIX}```

## Using the build
Activating the ```conda``` environment is not necessary once everything is built. You can use the installation normally by setting up your paths with ```build/setpaths.sh``` (or the ```csh``` equivalent).

## Updating environments
The general environment specification, `cctbx.devenv.yml`, is consolidated across platforms by using the [`conda-devenv` package](https://conda-devenv.readthedocs.io). Using Jinja2 syntax, platform specific dependencies can be added similar to how `conda-build` specifies dependencies. This consolidated environment file can be used to update the versioned environment files.

A standard conda environment yml file for a platform can be created by using the `--print` flag. This standard environment file can be modified to change the versions of specific dependencies. Additionally, because the channel order will be sorted alphabetically, the channel order should be modified to the desired order.
