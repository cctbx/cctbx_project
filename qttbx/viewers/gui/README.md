## Phenix Viewer
A molecular viewer using Molstar in a QtWebEngineView window.

### Installation
Assuming you are already working on an active Phenix/CCTBX development branch:
1. Get the branch
```sh
git fetch -a
```
2. Checkout the branch
```sh
git checkout molstar-veiwer-base
```
3. Install dependencies from conda (or mamba)
```sh
mamba -c conda-forge install pyside2 nodejs qtawesome ipykernel qt-webengine qtconsole-base
```
4. Change to Phenix modules directory
5. Download the Molstar branch
```sh
git clone git@github.com:phenix-project/phenix-molstar.git
```
6. Change to the `phenix-molstar` directory
7. Install nodejs dependencies
```sh
npm install
```
8. Compile Typescript to Javascript
```sh
npm run build
9. Refresh the build to register the command line files.
```
10. Start the viewer with a model as an argument
```sh
phenix.start_molstar 1yjp.cif
```
Alternatively, run directly if command line files not recognized:
```sh
phenix.python <MODULES_DIRECTORY>/cctbx_project/qttbx/command_line/start_molstar.py 1yjp.cif
```

### Development guide
The code is structured with the intention to follow a strict model-view-controller (MVC) architecture, to the degree that each component has its own subdirectory in the source code. 
```
cctbx_project/
└── qttbx/
    ├── command_line/
    │   └── start_molstar.py
    ├── programs/
    │   └── start_molstar.py
    └── viewers/
        ├── molstar/
        │   └── molstar.py
        └── gui/
            ├── README.md
            ├── apps/
            │   └── molstar_base_app.py
            ├── controller/
            ├── state/
            └── view/
```
