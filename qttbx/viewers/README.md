## Phenix viewers 
Integration between cctbx data structures, a QT gui, and external viewers such as ChimeraX and Molstar

### Switch to branch
```bash
git fetch --all
git switch ChimeraXSelectionViewer
```

### Installation of Molstar Viewer
See [here](https://github.com/phenix-project/phenix-molstar)

### Installation of additional dependencies
```
python -m pip install qtconsole
```

### Enable command line tool
Run from Phenix installation directory
```bash
python bootstrap.py build  --builder=phenix --use_conda --nproc=16
```


### Run a minimal molstar viewer to view some file
```bash
phenix.start_molstar 1yjp.cif
```


### Run a demo with all existing features
```bash
phenix.start_viewer_demo 1yjp.cif
```
