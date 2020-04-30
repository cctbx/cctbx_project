## INSTALL

Assume you've built CCTBX in a folder

```
export XTAL=/path/to/cctbx_root
```

Inside ```cctbx_root``` are the usual ```build``` and ```modules``` folders.


You can install this branch by navigating to the modules folder

```
cd $XTAL/modules/cctbx_project
git checkout diffBragg
git pull    # make sure its up to date
```

Then you should run the build script

```
cd $XTAL/build
make
```

After that you can then run the tests by creating an empty folder, and running the tests from that folder

```
mkdir $XTAL/test_folder  # for example
cd $XTAL/test_folder
libtbx.run_tests_parallel module=simtbx nproc=12
```

## About
diffBragg is made up of tools designed to refine a forward model of the scattering experiment, and it is built on the backbone of the nanoBragg forward simulation package.

#### layout of the code:

### ```src/```
Inside ```src/``` are the low level programs for computing the derivative of the expected scattering with respect to the model parameters. Here we use boost python to provide an interface to these programs (```diffBragg_ext.cpp```) so that this low level code is visible to Python.

The class ```diffBragg``` inherits from ```nanoBragg``` defined in (```nanoBragg_ext.cpp```). The function ```add_diffBragg_spots``` is akin to ```add_nanoBragg_spots``` however in addition to simulating the expected scattering, we simulate the derivatives of the expected scattering at each pixel. There is a list of derivative managers that are used to keep track of the "derivative pixels" for each parameter. The derivative manager list order is hard coded, such that 

* 0=rotX
* 1=rotY
* 2=rotZ
* 3-8=unit cell
* 9=Ncells (mosaic parameter m)
* 10=Fhkl
* 11=originZ

### ```sim_data.py```
Forward modeling depends on certain attributes, for example crystal models, beam models, detector models. The class ```SimData``` attempts to be a high level class to help with the handling of these various attributes. A main property of ```SimData``` is the ```D``` property, which points to an instance of ```diffBragg```. It's important to run ```SimDataInstance.D.free_all()``` when you are done with the instance.

### ```nanoBragg_crystal.py, nanoBragg_beam.py```
These define high level crystal and beam objects that are attributes of a ```SimData``` instance. For example a ```nanoBragg_crystal``` has an associated ```miller_array```, and a ```nanoBragg_beam``` has an associated ```spectrum```. 

### ```refiners/global_refiner.py```
This defines the main refiner class. It computes the outer derivative of teh Likelihood (which depends on the inner derivative computed by the low level ```diffBragg``` class  in ```src/diffBragg.cpp```. Other refiners are mainly used for testing and are historical at this point, and will eventually be removed.

### ```refiners/pixel_refinement.py```
This defines the base class for ```GlobalRefiner``` and it manages the LBFGS interface. The ```__init__``` method of this class has A LOT of attributes that can be adjusted for refinement.


### ```tests/```
The most comprehensive test is ```tst_diffBragg_global_refine.py```. It is a script with arguments, the default arguments for runing the test can be seen in the file ```simtbx/run_tests```. 

### ```refiners/crystal_systems```
Depending on the lattice symmetry, we will refine anywhere from 1 (cubic) to 6 (triclinic) unit cell parameters. The classes here keep track of the number of free parameters, as well as provide useful terms needed to compute the derivative of the expected scattering with respect to the  various unit cell parameters. 




