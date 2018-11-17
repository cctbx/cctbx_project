# Conda Environments

**TESTING (not for use in production)**
- macOS may fail due to ```DYLD_LIBRARY_PATH``` being set in the dispatcher

These files specify a ```conda``` environment that contains the dependencies required for CCTBX projects. They are named according to ```<program>_<python version>_<platform>.txt```. So for the file, ```cctbx_py27_linux-64.txt```, the environment is for building the base CCTBX package using Python 2.7 on a 64-bit Linux distribution. The ```<platform>``` follows the ```conda``` pattern for naming platforms, so ```linux-64``` refers to 64-bit linux, ```osx-64``` refers to 64-bit macOS, and ```win-64``` refers to 64-bit Windows.
## Building with ```conda```
Assuming you have already downloaded the relevant source files with ```bootstrap.py hot update --builder=<builder>``` into the ```<root directory>```, these directions let you build CCTBX packages with dependencies in a```conda``` environment.

1) Create a new environment,
```conda create -n <environment name> --file cctbx_py27_linux-64.txt```
2) Activate new environment,
```conda activate <environment name>```
3) Build!
```cd <root directory>```
```python modules/cctbx_project/libtbx/auto_build/bootstrap.py build --builder=<builder> --with-python=`which python` ```

The ```--with-python``` flag is used to specify which ```python``` to use. The ```base``` option is not used because the ```conda``` should environment contain the necessary dependencies.
## Using the build
Activating the ```conda``` environment does not seem to be necessary once everything is built. You can use the installation normally by setting up your paths with ```build/setpaths.sh``` (or the ```csh``` equivalent).
