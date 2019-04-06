## Diffraction Experiment Toolbox

This code has now moved into a separate repository.
If you set up your build environment after April 2019 it should already include
the necessary changes.


### Updating developer installations

You can pick up dxtbx from its new location by running the command
```bash
libtbx.install dxtbx
```
which will install dxtbx in the modules/ directory. You should also run
```bash
cd $(libtbx.find_in_repositories cctbx_project)
git clean -diX dxtbx
```
to remove any `.pyc` and other python runtime debris from the previous dxtbx
path, as these files could interfere with python package discovery.
The clean command will ask for confirmation (enter `1`) before deleting files.


### Code migration

If you have any cctbx\_project development branches that touch the dxtbx
directory and need help in transferring commits over, please open an issue on
Github and we are happy to help.
