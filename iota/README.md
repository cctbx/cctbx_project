## IOTA

This code has now moved into a separate repository.
If you set up your build environment after April 2020 it should already include
the necessary changes.


### Updating developer installations

You can pick up iota from its new location by running the command
```bash
libtbx.install iota
```
which will install iota in the modules/ directory. You should also run
```bash
cd $(libtbx.find_in_repositories cctbx_project)
git clean -diX iota
```
to remove any `.pyc` and other python runtime debris from the previous iota
path, as these files could interfere with python package discovery.
The clean command will ask for confirmation (enter `1`) before deleting files.
