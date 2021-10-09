# DiffBragg Examples
------------------

# Setup

These instructions assume you have a working CCTBX environment, and that the command `libtbx.python` is in your path. Also, this tutorial will assume you are using NERSC (instructions for private linux clusters will be provided in due time).

The easiest way to to learn diffBragg is by using simulated data. You could use the tools provided in CCTBX to simulate data readily, however for this tutorial, there is another github repository you can download and easily get some images. Start by navigating to your CCTBX buidld's modules folder on your computer. If you are unsure where this is, try issuing the following command

```
$ libtbx.python -c 'import cctbx;print(cctbx.__file__.split("/cctbx_project")[0])'
/Users/dermen/CrystalNew/modules
$ export MOD=/Users/dermen/CrystalNew/modules
```

Now, in the modules folder, download the cxid9114 repository and use `git-lfs` to bring in some data files

```
cd $MOD 
git clone https://github.com/dermen/cxid9114.git
# install git-lfs (if on nersc, just load the module
cd ~/Crystal/modules/cxid9114
module load git-lfs
cd ~/Crystal/modules/cxid9114
git lfs install
git lfs fetch
git lfs pull
```

Also, install `mpi4py`.

# Simulating images

The first step in simulating data using `cxid9114` is to create a background image:

```
libtbx.python $MOD/cxid9114/sim/d9114_mpi_sims.py -odir . --bg-name mybackground.h5 --make-background   --sad 
```

To view the background image, install the proper image format provided in the repository.

```
dxtbx.install_format  -u $MOD/cxid9114/format/FormatD9114.py
dials.image_viewer mybackground.h5 color_scheme=heatmap brightness=150
```

should produce the image:

[INSERT]

Now, we can simulate the diffraction using this image as background. The command is


```
srun -n8 -c2 libtbx.python $MOD/cxid9114/sim/d9114_mpi_sims.py  -o test -odir poly_images --add-bg --add-noise --profile gauss --bg-name mybackground.h5 -trials 13  --oversample 0 --Ncells 10 --xtal_size_mm 0.00015 --mos_doms 1 --mos_spread_deg 0.01  --saveh5 --readout  --masterscale 1150 --sad --bs7real --masterscalejitter 115 --overwrite --gpu -g 8
```

the above command will simulate 13 diffraction patterns per rank from randomly oriented crystals. If you are not on NERSC, then drop the `srun -n8 -c2` prefix , and if you do not have GPU support, drop the `--gpu -g8` flags at the end of the command, and instead add the flag `--force-mono` in order to only simuilate a single wavelength per pattern. The images can also be opened with the image viewer:

```
dials.image_viewer  poly_images/job0/test_rank0_data0_fluence0.h5
```

Now that we have patterns, we can process them using diffBragg wrapper script `hopper_process`.

### using `hopper_process`

We have simulated 104 images, and now we shall process them using the command line tool hopper_process, a program akin to dials.stills_process, however if slight modifcations. Crucially, we will be disabling outlier rejection and all forms of refinement that are usually done during `dials.stills_process` analysis, in favor of the pixel refinement tools in diffBragg, which are wrapped in `hopper_process`. 

#### Monochromatic `hopper_process`

If X-ray spectra are available for your data, then they should be encoded in the image format. In the simulated data images, we store the spectra in the hdf5 datasets `"wavelengths"` and `"spectrum"`. We can grab the spectra manually as a demonstration of the API

```python
import h5py
h = h5py.File("poly_images/job0/test_rank0_data0_fluence0.h5", "r")
a,b = h['wavelengths'][()], h['spectrum'][()]. # `spectrum` was a poor choice of name here, these are just the weights
h.close()

import pylab as plt
plt.plot( a, b)
plt.xlabel("wavelength $\AA$",fontsize=12)
plt.gca().tick_params(labelsize=11)

# check the mono wavelength is set correctly
import dxtbx
dx = dxtbx.load("poly_images/job0/test_rank0_data0_fluence0.h5")
ave_wave = dx.get_beam(0).get_wavelength()
plt.subplots_adjust(bottom=0.15)
plt.vlines( ave_wave, b.min(), b.max(),ls='--', color='tomato')
plt.show()
```

<p align="center">
<img src="https://user-images.githubusercontent.com/2335439/136639491-caeb9490-e9f1-492b-8ac2-e531184a544f.png" />
</p>

These are also encoded in the format class `FormatD9114.py` which was installed above using `dxtbx`, so check in that file for examples on how to plug in your own format class.

Lets assume that we do not have spectra - we can still run hopper_process using the mono-wavelength associated with each image. To do that, we execute the following command

```
DIFFBRAGG_USE_CUDA=1 srun -n 8 -c2  simtbx.diffBragg.hopper_process process.phil "poly_images/job*/*.h5"  mp.method=mpi mp.nproc=8 num_devices=8  dispatch.integrate=True spectrum_from_imageset=False  output.output_dir=poly_images/procMono
``` 

The configuration file `process.phil` contains

<details>
  <summary>Example process.phil</summary>
  
```
spotfinder {
  threshold.algorithm=dispersion
  threshold.dispersion.gain=28
  threshold.dispersion.global_threshold=40
  threshold.dispersion.kernel_size=[2,2]
  threshold.dispersion.sigma_strong=1
  threshold.dispersion.sigma_background=6
  filter.min_spot_size=2
}

indexing {
  method=fft1d
  known_symmetry.unit_cell=79.1,79.1,38.4,90,90,90
  known_symmetry.space_group=P43212
  refinement_protocol.mode=None
  stills {
    rmsd_min_px=4 
    refine_all_candidates=False 
  }
}

integration.summation.detector_gain=28

# hopper process config:
save_modelers=True
silence_dials_loggers=True
partial_correct=False
reidx_obs=False  

diffBragg {
  no_Nabc_scale=False
  method="L-BFGS-B"
  use_restraints=False
  space_group=P43212
  spectrum_from_imageset = True 
  downsamp_spec {
    skip = True
  }
  roi {
    shoebox_size=12
    fit_tilt=True
    fit_tilt_using_weights = False
    hotpixel_mask = None 
    reject_edge_reflections = False
    reject_roi_with_hotpix = False
    pad_shoebox_for_background_estimation=10
  }
  refiner {
    adu_per_photon = 28
    sigma_r=3
  }
  simulator {
    oversample=2
    crystal.has_isotropic_ncells = False
    structure_factors.mtz_name = iobs_all.mtz 
    structure_factors.mtz_column = "I(+),SIGI(+),I(-),SIGI(-)"
    init_scale = 1
    beam.size_mm = 0.001
    detector.force_zero_thickness = True
  }
  init {
    Nabc=[13,14,15]
    G=100
  }
  mins {
    Nabc=[3,3,3]
    detz_shift=-1.5
    RotXYZ=[-15,-15,-15]
    G=0
  }
  maxs {
    RotXYZ=[15,15,15]
    Nabc=[1600,1600,1600]
    G=1e12
    detz_shift=1.5
  }
  sigmas {
    RotXYZ=[1e-3,1e-3,1e-3]
  }
  fix.detz_shift=True
}
```
</details>

however any command line parameters supercede whats in the file. Notice the command line argument `spectrum_from_imageset=False`. This tells `hopper_process` to read the wavelength from the beam file, and ignore any X-ray spectra that might be present in the data. By setting this flag, we are refining a monochromatic model for each shot. The command takes 71 seconds to run on a single NERSC compute node utilizing 8 GPUs and 8 processors (1 GPU per process), and optimizing models for 104 shots. We have prepared a simple script called `quick_detresid.py` here for analyzing the results


<details>
  <summary>quick_detresid.py</summary>
  
```python
"""quick_detresid.py"""
import glob
import sys
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
import numpy as np
from dxtbx.model import ExperimentList, Experiment
import os

dirname = sys.argv[1]
fnames = glob.glob(dirname+ "/*indexed.refl")
if not fnames:
    print("no fnames")
    exit()
detname = fnames[0].replace("indexed.refl", "refined.expt")
assert os.path.exists(detname)
El = ExperimentListFactory.from_json_file(detname, check_format=False)
dummie_explist = ExperimentList()

D = El[0].detector
B = El[0].beam
R = flex.reflection_table()
for i,f in enumerate(fnames):
    dummie_exp = Experiment()
    dummie_exp.detector = D
    dummie_exp.beam = B
    dummie_explist.append(dummie_exp)
    r = flex.reflection_table.from_file(f)
    eid = r.experiment_identifiers()
    for k in eid.keys():
        del eid[k]
    r['id'] = flex.int(len(r), i)
    R.extend(r)
    print(i)
    
pid = R['panel']

cents = []
sel = []
for i in range(len(R)):
    p = pid[i]
    x,y,z = R[i]['xyzcal.px']
    if np.isnan(x):
        sel.append(False)
        x = y = z = 0
    else:
        sel.append(True)
    xm,ym = D[p].pixel_to_millimeter((x,y))
    cents.append((xm,ym,0))
    print(i)
    
R['xyzcal.mm'] = flex.vec3_double(cents)
R2 = R.select(flex.bool(sel))
R2['id'] = flex.int(len(R2),0)
print("Good ones %d / %d" % (len(R2), len(R)))
ref = os.path.join(dirname, "comb.refl")
exp = os.path.join(dirname, "comb.expt")
R2['delpsical.rad'] = flex.double(len(R2), 0)
R2.as_file(ref)

dummie_explist.as_file(exp)
print("wrote %s, %s" % (ref, exp))

### Call detresiduals
detresid_phil="""
dot_size=3 
colormap=Oranges
repredict_input_reflections=False 
hierarchy_level=1 
residuals {
  plot_max=0.3 
  exclude_outliers_from_refinement=False 
}

plots {
  include_offset_dots=True 
  include_scale_bar_in_pixels=1 
  positional_displacements=True
  deltaXY_by_deltapsi=False
  per_image_RMSDs_histogram = False
  include_radial_and_transverse = False
  unit_cell_histograms = False
  pos_vs_neg_delta_psi=False
}
"""
o = open("_detresid.phil", "w")
o.write(detresid_phil)
o.close()

os.system("cctbx.xfel.detector_residuals %s/comb.* _detresid.phil" % dirname)
```
</details>

It is a simple wrapper which combines the relevant outputs and calls the command line program `cctbx.xfel.detector_residuals`. We are interested at this point to see how well the refined model predicts the strong spot observations. Issue the command

```
libtbx.python quick_detresid.py poly_images/procMono
```

and you will see an image display, as well as some numbers print to the string, indicating how well the monochromatic model predicts the data:

```
RMSD (microns) 94.76634677545071
Overall radial RMSD (microns) 65.9469202178172
Overall transverse RMSD (microns) 68.05633104237863
```

<p align="center">
<img src="https://user-images.githubusercontent.com/2335439/136640760-7fb111b0-1ff1-48d6-b5d2-0c8296b691cf.png" />
</p>

Not bad for a CSPAD with 110 micron pixels. Can we do better using a polychromatic model ? All we need to do to test this is drop the flag ```spectrum_from_imageset=False``` from the command line (its set to `True` in the ```process.phil```):

```
DIFFBRAGG_USE_CUDA=1 srun -n 8 -c2  simtbx.diffBragg.hopper_process process.phil "poly_images/job*/*.h5"  mp.method=mpi mp.nproc=8 num_devices=8  dispatch.integrate=True   output.output_dir=poly_images/procPoly
libtbx.python quick_detresid.py poly_images/procPoly

#RMSD (microns) 61.79812988365175
#Overall radial RMSD (microns) 46.86256672413951
#Overall transverse RMSD (microns) 40.285340970907185
```

<p align="center">
<img src="https://user-images.githubusercontent.com/2335439/136640765-d2d0b274-cd0c-4613-9ea8-42e2f497a418.png" />
</p>

The numbers that print to the screen now are slightly more optimized, owing to the fact that we used a polychromatic model which is more in-line with the data. In fact these numbers represent the best we can do when we know our detector geometry perfectly. This time, the data took 200 seconds to process.


### Fixing a faulty geometry

Geometry refinement is currently using the `lmfit` module. Install it using pip:

```
libtbx.python -m pip install lmfit
```

We have prepared a faulty experimental geometry with which to process the data, taken directly from real expeirmental errors associated with the CSPAD geometry. This is to simulate the scenario where the geometry is not well known, and one wishes to optimize it using pixel refinement. One must extract it from its raw form and write it to disk, using the following simple script

```python
from cxid9114.geom.multi_panel import CSPAD2
from dxtbx.model import Experiment, ExperimentList
El = ExperimentList()
E = Experiment()
E.detector = CSPAD2
El.append(E)
El.as_json("badGeo.expt")
```

Now we can process the simulated images using the command

```
DIFFBRAGG_USE_CUDA=1 srun -n 8 -c2  simtbx.diffBragg.hopper_process process2.phil "poly_images/job*/*.h5"  mp.method=mpi mp.nproc=8 num_devices=8  dispatch.integrate=True  output.output_dir=poly_images/procBad reference_geometry=badGeo.expt 
```

The flag `reference_geometry=badGeo.expt` forces the geometry file stored in badGeo.expt to override the geometry defined by the image format. In this way, all models we optimize are subjected to the errors in the geometry, and this should be reflected in the prediction offsets. Indeed we find that the numbers are much worse, and the geometry is visibly distorted

```
libtbx.python quick_detresid.py poly_images/procBad

#RMSD (microns) 117.83984508972222
#Overall radial RMSD (microns) 92.76564120886626
#Overall transverse RMSD (microns) 72.66887161555252
``` 

<p align="center">
<img src="https://user-images.githubusercontent.com/2335439/136641295-b2f4f727-0faa-4da2-b3be-df9693f713a6.png" />
</p>

In order to fix the bad geometry, we can use the program `geometry_refiner.py` located at `$MOD/cctbx_project/simtbx/diffBragg/geoemtry_refiner.py`. The script `hopper_process` does not currently support optimization of the detector geometry, other than a per-shot shift in the overall detector position (which we have kept fixed so as not to interfere with the geometry refinement we are about to do). `hopper_process` operates on single images, whereas geometry refinement benefits from jointly optimizing multiple images together. As it is a large refinement, we will regulate it with parameter restraints. To determine the proper restraint targets, examine the summary pandas file located in the `poly_images/procBad` folder using the following script

```python
import pandas
import numpy as np
import pylab as plt

df = pandas.read_pickle("poly_images/procBad/hopper_process_summary.pkl")
df[['a','c']].hist(bins=30)
na,nb,nc = np.vstack(df.ncells).T
df['Nc'] = nc
df['Nb'] = nb
df['Na'] = na
df['logG'] = np.log10(df.spot_scales)
df[['Na','Nb', 'Nc', 'logG']].hist(bins=30)

print(df[['a','c','Na','Nb', 'Nc', 'logG']].median().to_string())

plt.show()

#a       79.313096
#c       38.507866
#Na       7.331476
#Nb       7.521334
#Nc       8.755147
#logG     6.970101
``` 

<p align="center">
<img src="https://user-images.githubusercontent.com/2335439/136643006-b2c5f2ab-fe36-4907-bef0-8aa73825ec57.png" />
</p>

<p align="center">
<img src="https://user-images.githubusercontent.com/2335439/136643034-817e9843-87e4-470b-8600-49af8d820099.png" />
</p>

The previous script produces histograms of the unit cell parameters, as well as the mosaic domain blocksize (Na,Nb,Nc), and the log of the per-shot scale factors (logG). From the histogram, we see suitable values to use for the restraint targets (real experimental data might not look so nice, but restraints still help!). 

This information is passed to `geometry_refiner.py` in the configuation file

<details>
  <summary>geom.phil</summary>

```
method="L-BFGS-B"
use_restraints=True
space_group=P43212
spectrum_from_imageset = True
downsamp_spec {
  skip = True
}
roi {
  shoebox_size=12
  fit_tilt=True
  fit_tilt_using_weights = False
  hotpixel_mask = None 
  reject_edge_reflections = False
  reject_roi_with_hotpix = False
  pad_shoebox_for_background_estimation=10
}
refiner {
  adu_per_photon = 28
  sigma_r=3
  panel_group_file = cspad_32panel.txt
}
simulator {
  oversample=2
  crystal.has_isotropic_ncells = False
  structure_factors.mtz_name = iobs_all.mtz 
  structure_factors.mtz_column = "I(+),SIGI(+),I(-),SIGI(-)"
  init_scale = 1
  beam.size_mm = 0.001
  detector.force_zero_thickness = True
}
mins {
  Nabc=[3,3,3]
  detz_shift=-1.5
  RotXYZ=[-15,-15,-15]
  G=0
}
maxs {
  RotXYZ=[15,15,15]
  Nabc=[1600,1600,1600]
  G=1e12
  detz_shift=1.5
}
betas {
  RotXYZ=1e-5
  Nabc=[1e3,1e3,1e2]
  G=1e5
  ucell=[1e-3,1e-3]
  detz_shift=1e-8
}
centers {
  Nvol=None
  G=1e7
  Nabc=[9,9,9]
  ucell=[79.3, 38.5]
}
sigmas {
  RotXYZ=[1e-2,1e-2,1e-2]
}
fix {
  detz_shift=True
  ucell = False
  Nabc = False
  RotXYZ = False
}

geometry {
  optimize = True
  betas {
    panel_rot=[1,1,1]
    panel_xyz=[1,1,1]
  }
  min {
    panel_rotations=-5,-1,-1
    panel_translations=-1,-1,-1
  }
  max {
    panel_rotations=5,1,1
    panel_translations=1,1,1
  }
}

```
</details>

Note that the phil scopes for `roi, refiner, simulator` are mostly copied from the `process.phil`, however there are more parameters now, namely those to do with parameter restraints. In particular, the `betas` are reciprocal variances for the restraint targets (small `betas` are more restrained), and `centers` are the restraint targets themselves, taken frmo the histograms. We know that the mosaic domain size will tend to lower values in the presence of systematic errors (smaller domain sizes lead to larger spot profiles, and larger spot profiles overlap more with the strong observations in the presence of geometry errors). All of the panel geometry corrections have restraint targets of `0` by default, as well as the crystal misorientations `RotXYZ`.

##### Panel groups
A user-defined 2-column text file defined which panels are treated as single contiguous units. The first column is the panel ID (0-63 for the CSPAD), and the second column is the group ID. The user has full control over the groupings. From examining the detector residuals images above, we can see the natural panel groping of the CSPAD is panels (0,1), (2,3), (4,5), etc. 

```python
with open("single_panel.txt", "w") as o:
    for pid in range(64):
        o.write("%d %d\n" % (pid, 0))

with open("cspad_32panel.txt", "w") as o:
    for pid in range(64):
        groupid = int(pid/2)
        o.write("%d %d\n" % (pid, group_id))

with open("cspad_quads.txt", "w") as o:
    for pid in range(64):
        groupid = int(pid/16)
        o.write("%d %d\n" % (pid, group_id))
```

Here we will use the 32-panel model defined in `cspad_32panel.txt`. One quirk is that the `panel_group_file=cspad_32panel.txt` parameter is defined the the `refiner` phil scope. This is a remnant of an older version of geometry refiner and will likely change, so keep an eye on `cctbx_project/simtbx/diffBragg/phil.py` !

To run geometry refinement, issue the command

```
DIFFBRAGG_USE_CUDA=1 srun -n8 -c2 libtbx.python $DDZ/refiners/geometry_refiner.py --phil geom.phil --cmdlinePhil  optimized_detector_name=optGeo.expt input_pkl=poly_images/procBad/hopper_process_summary.pkl lbfgs_maxiter=2000 num_devices=8 geometry.first_n=80
```

The `input_pkl` is the file written by hopper_process if it terminates successfully. If for some reason hopper_process terminates prematurely, one can create this file following this script

```python
import pandas
import glob
fnames = glob.glob("poly_images/procBad/pandas/rank*/*pkl")
df = pandas.concat([pandas.read_pickle(f) for f in fnames])
df.to_pickle("geomRef_input.pkl") # for example

# *******************************************************************
# One can also run geometry refinement on the integrated reflections.
# This is useful to model spots evenly across all panels, 
# and out to the corners of the detector.
# *******************************************************************
df["integrated_refls"] = [f.replace("_indexed", "_integrated") for f in df.stage1_refls]
import os
# ensure the integrated pickles are present!
df = df.loc[ [os.path.exists(f) for f in df.integrated_refls]]
print("%d of those records had integration results" % len(df))
df.to_pickle("geomRef_integ.pkl")
```

The geometry refiner will load experiments and reflection tables that are stored in the input pickle in the columns `"exp_name"` and `"stage1_refls"`, however there is a PHIL parameter `refls_key` for choosing alternate reflection tables. If you had run the above script, setting `refls_key=integrated_refls` and `input_pickle=geomRef_itneg.pkl` would launch a refinement that models all pixels in the integration shoeboxes. The flag `geometry.first_n=80` specifies the number of diffraction shots to use in the refinement.

Now, we can re-run `hopper_process` using this optimized geometry by adding the flag `reference_geometry=optGeo.expt`:

```
DIFFBRAGG_USE_CUDA=1 srun -n 8 -c2  simtbx.diffBragg.hopper_process process2.phil "poly_images/job*/*.h5"  mp.method=mpi mp.nproc=8 num_devices=8  dispatch.integrate=True  output.output_dir=poly_images/procOpt reference_geometry=optGeo.expt
libtbx.python quick_detresid.py poly_images/procOpt

#RMSD (microns) 66.94415897822252
#Overall radial RMSD (microns) 52.0290261621068
#Overall transverse RMSD (microns) 42.124824414167406
```

<p align="center">
<img src="https://user-images.githubusercontent.com/2335439/136643005-b32b09bb-547b-4360-8f7f-73f14125b243.png" />
</p>

The result shows significant improvement with the optimized geometry, more so than even the monochromatic models that were subjected to zero detector inaccuracy!

This particular dataset has rather fat spots that are dominated by the mosaic domain size as opposed to the spot spectra, but nonetheless we get better results using polychromatic models. Plotting the average prediction offset as a function of resolution shows this readily. The script

<details>
  <summary>pred_offsets.py</summary>

```python
"""pred_offsets.py"""
import glob
import os
from pylab import *
from dials.array_family import flex
from joblib import Parallel, delayed

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("dirnames", type=str, nargs="+", help="hopper_process output folders")
parser.add_argument("-j", type=int, default=1, help="number of procs")
parser.add_argument("-nbins", type=int, default=10, help="number of resolution bins")
args = parser.parse_args()


NJ=args.j


def xy_to_polar(refl,DET, dials=False):
    x, y, _ = refl["xyzobs.px.value"]
    if dials:
        xcal, ycal, _ = refl["dials.xyzcal.px"]
    else:
        xcal, ycal, _ = refl["xyzcal.px"]

    pid = refl['panel']
    panel = DET[pid]
    x,y = panel.pixel_to_millimeter((x,y))
    xcal,ycal = panel.pixel_to_millimeter((xcal,ycal))

    xyz_lab = panel.get_lab_coord((x,y))
    xyz_cal_lab = panel.get_lab_coord((xcal, ycal))

    diff = np.array(xyz_lab) - np.array(xyz_cal_lab)

    xy_lab = np.array((xyz_lab[0], xyz_lab[1]))
    rad = xy_lab / np.linalg.norm(xy_lab)
    tang = np.array([-rad[1], rad[0]])

    rad_component = abs(np.dot(diff[:2], rad))
    tang_component = abs(np.dot(diff[:2], tang))
    pxsize = panel.get_pixel_size()[0]
    return rad_component/pxsize, tang_component/pxsize


def main(jid, njobs, dirname):

    from dxtbx.model import ExperimentList

    fnames = glob.glob("%s/*indexed.refl" % dirname)
    print("%d fnames" % len(fnames)) 
    assert fnames
    detpath = fnames[0].replace("indexed.refl", "refined.expt")
    assert os.path.exists(detpath)
    DET = ExperimentList.from_file(detpath, False)[0].detector

    all_d = []
    all_r = []
    all_t = []
    reso = []
    for i_f, f in enumerate(fnames):
        if i_f % njobs != jid:
            continue
        R = flex.reflection_table.from_file(f)
        if len(R)==0:
            continue
        xyobs = R['xyzobs.px.value'].as_numpy_array()[:,:2]
        xycal = R['xyzcal.px'].as_numpy_array()[:,:2]
        reso += list( 1./np.linalg.norm(R['rlp'], axis=1))
        d = np.sqrt(np.sum( (xyobs -xycal)**2, 1))
        all_d += list(d)
        rad,theta = zip(*[xy_to_polar(R[i_r],DET,dials=False) for i_r in range(len(R))])
        all_r += list(rad)
        all_t += list(theta)
        print(i_f)
    return all_d, all_r, all_t, reso


def results_from_folder(dirname, nbins):
    results = Parallel(n_jobs=NJ)(delayed(main)(j,NJ,dirname) for j in range(NJ))

    all_d, all_r, all_t,  reso = [],[],[],[]
    for d,r,t, dspacing in results:
        all_d += d
        all_r += r
        all_t += t
        reso += dspacing

    bins = [ b[0]-1e-6 for b in np.array_split(np.sort(reso), nbins)] + [max(reso)+1e-6]
    digs = np.digitize(reso, bins)

    all_d = np.array(all_d)
    all_r = np.array(all_r)
    all_t = np.array(all_t)
    reso = np.array(reso)
    ave_d, ave_r,  ave_t,  ave_res =[],[],[],[]
    for i_bin in range(1, nbins+1):
        sel = digs==i_bin
        ave_d.append( np.median(all_d[sel]))
        ave_r.append( np.median(all_r[sel]))
        ave_t.append( np.median(all_t[sel]))

        ave_res.append( np.median(reso[sel]))
    return ave_d, ave_r, ave_t, ave_res

all_d = []
all_r = []
all_t= []
for dirname in args.dirnames:
    d,r,t, ave_res = results_from_folder(dirname, args.nbins)
    all_d.append(d)
    all_r.append(r)
    all_t.append(t)

for vals_series, title in [(all_d, "overall"), (all_r, "radial component"), (all_t, "tangential component")]:

    figure()
    gca().set_title(title)
    from itertools import cycle
    colors = cycle(["tomato", "chartreuse", "plum"])
    markers = cycle(["o", "s", "*", ">"])
    for i_d, dirname in enumerate(args.dirnames):
        vals = vals_series[i_d]
        plot(vals[::-1], color=next(colors), marker=next(markers), mec='k', label=dirname)
    xticks = range(args.nbins)
    xlabels = ["%.2f" % r for r in ave_res]
    gca().set_xticks(xticks)
    gca().set_xticklabels(xlabels[::-1], rotation=90)
    gcf().set_size_inches((5,4))
    subplots_adjust(bottom=0.2, left=0.15, right=0.98, top=0.9)
    gca().tick_params(labelsize=10, length=0) # direction='in')
    grid(1, color="#777777", ls="--", lw=0.5)
    xlabel("resolution ($\AA$)", fontsize=11, labelpad=5)
    ylabel("prediction offset (pixels)", fontsize=11)
    leg = legend(prop={"size":10})
    fr = leg.get_frame()
    fr.set_facecolor("bisque")
    fr.set_alpha(.5)
    gca().set_facecolor("gainsboro")

show()
```
</details>

```
libtbx.python pred_offsets.py  poly_images/procMono/ poly_images/procPoly/ poly_images/procBad/ poly_images/procOpt/
```

<p align="center">
<img src="https://user-images.githubusercontent.com/2335439/136645704-d5958f38-71e2-4cb7-b5ff-4b3085baa896.png" />
</p>

<p align="center">
<img src="https://user-images.githubusercontent.com/2335439/136645710-2a594c09-b790-4c3a-a225-561446387a7e.png" />
</p>

<p align="center">
<img src="https://user-images.githubusercontent.com/2335439/136645711-7ecd1ace-b153-44af-8967-5fbfcfe5a7bf.png" />
</p>



##  more on hopper_process

The script runs dials.stills_process with a few alterations

* diffBragg refinement is performed after indexing, provided that `skip_hopper=False` . 
* `reidx_obs=True` will re-index the strong spot observations after running the normal stills indexing algorithm, *prior* to running refinement. This is useful when obtaining an indexing solution warrants using a high-res cutoff. In such a case, one can detect strong spots out to the corners of the detector, use the indexing_refinement phil param [TODO lookup name] to limit the spots which are fed into indexing, and then re-index the strong spots out to the corners of the camera to grab more spots for diffBragg refinement. 
* After refinement, the diffBragg model is used to predict integration positions on the detector, and then the dials integration program is used to compute integrated spot intensities.

Also, if `save_models=True`, then a modellers folder will be written   
that contains information about the per-shot modellers



