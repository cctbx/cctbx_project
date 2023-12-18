## Phenix viewers 
Integration between cctbx data structures, a QT gui, and external viewers such as ChimeraX and Molstar

# Installation

### Switch to branch
```bash
git fetch --all
git switch ChimeraXSelectionViewer
```

### Installation of Molstar Viewer
See [here](https://github.com/phenix-project/phenix-molstar)
```

### Optional: Installation of additional dependencies
```
python -m pip install qtconsole
```

### Enable command line tool
Run from Phenix installation directory
```bash
python bootstrap.py build  --builder=phenix --use_conda --nproc=16
```


### Open a file from the command line
```bash
phenix.start_viewer 1yjp.cif
```
