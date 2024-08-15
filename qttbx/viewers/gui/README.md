## Phenix Viewer
A molecular viewer using Molstar in a QtWebEngineView window.

### Installation
Assuming you are already working on an active Phenix/CCTBX development branch:
1. Get the branch
```sh
git fetch -a
```
2. Checkout the branch
```
git checkout molstar-veiwer-base
```
3. Install dependencies from conda (or mamba)
```
mamba -c conda-forge install pyside2 nodejs qtawesome ipykernel qt-webengine qtconsole-base
```
4. Change to Phenix modules directory
5. Download the Molstar branch
```
git clone git@github.com:phenix-project/phenix-molstar.git
```
6. Change to the `phenix-molstar` directory
7. Install nodejs dependencies
```
npm install
```
8. Compile Typescript to Javascript
```
npm run build
9. Refresh the build to register the command line files.
```
10. Start the viewer with a model as an argument
```
phenix.start_molstar 1yjp.cif
```
Alternatively, run directly if command line files not recognized:
```
phenix.python <MODULES_DIRECTORY>/cctbx_project/qttbx/command_line/start_molstar.py 1yjp.cif
```

### Development guide
The code is structured with the intention to follow a strict model-view-controller (MVC) architecture, to the degree that each component has its own subdirectory in the source code in the ```qttbx/viewers/gui``` subdirectory. 
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
            ├── model/
            └── view/
```
The contents of each directory are described below:

1. ```qttbx/command_line/start_molstar.py``` A stub to register the program as a command line tool
2. ```qttbx/programs/start_molstar.py``` The starting point of the program using the ProgramTemplate class. All parameters parsing and file IO should be done here before the GUI is instantiated.
3. ```qttbx/viewers/molstar/``` This directory contains the python code necessary to interface with the Molstar web app running in a QtWebEngineView window. The ```molstar.py``` file defines a Python api that can be used to affect the Molstar graphics by running corresponding Javascript snippets. Anything in the rest of the GUI that affects the Molstar graphics must come through the ```MolstarGraphics``` class in this file.
4. ```qttbx/viewers/gui/``` Everything in this subdirectory is for the Qt GUI.


