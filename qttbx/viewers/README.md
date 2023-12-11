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


### Running
The entry point is qttbx/command_line/start_gui.py
```bash
python bootstrap.py build  --builder=phenix --use_conda 
phenix.start_gui
```
