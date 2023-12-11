## Phenix viewers 
Integration between cctbx data structures, a QT gui, and external viewers such as ChimeraX and Molstar

### Switch to branch
```bash
git pull
git checkout ChimeraXSelectionViewer
```

### Installation of Molstar Viewer
See [here](https://github.com/phenix-project/phenix-molstar)

```bash
cd qttbx/viewers/molstar
git clone https://github.com/phenix-project/phenix-molstar/
cd phenix-molstar

npm install
npm run build
```

### Enable command line tool
```bash
python bootstrap.py build  --builder=phenix --use_conda --nproc=16
```


### Run a minimal molstar viewer
```bash
phenix.start_molstar data/1yjp.cif
```


### Run a demo with all existing features
```bash
phenix.start_viewer_demo data/1yjp.cif
```
