from __future__ import absolute_import, division, print_function
import os
import sys
import glob
import libtbx.option_parser
from xfel.amo.pnccd_ana                 import mpi_fxs_index
from xfel.amo.pnccd_ana                 import mpi_fxs_c2
from xfel.amo.pnccd_ana                 import mpi_fxs_bg
from xfel.amo.pnccd_ana                 import mpi_fxs_calib
from xfel.amo.pnccd_ana                 import mpi_fxs_mask

def launch(argv=None) :

  """Function to launch: Single CPU, Multi-Processor interactive jobs or MPI batch jobs
       for the purpose of:

       1. Indexing                (mpi_fxs_index.py)
       2. Masking                 (mpi_fxs_mask.py)
       3. Background computation  (mpi_fxs_bg.py)
       4. Calibration             (mpi_fxs_calib.py)
       5. C2 computation          (mpi_fxs_c2.py)

     Example tcsh-file for launching Single, Interactive or Batch jobs:

     #################################################################
     #!/bin/tcsh

     # Provide Run number in input (eg ./run.csh 111)
     #echo $1

     # Batch parameters
     set processing_q  = psanaq
     set nr_cpu        = 120
     set job_name      = run_$1

     # Processing parameters
     set file_type     = idx
     set job_type      = correlation
     set experiment    = amok5415
     set detector      = pnccdFront
     set distance      = 369.0
     set mask_thr      = 2.0
     set wavelength    = 7.094
     set x             = 513
     set y             = 662
     set dx            = 5
     set dy            = 5

     set mask          = /reg/neh/home1/malmerbe/develop/amok5415_files/processing/masks/Mask_map_61_2.dat
     set index         = /reg/neh/home1/malmerbe/develop/amok5415_files/index_files/Selected_run111_COM.dat

     set first         = 0
     set last          = 20

     # Launch Batch Job

     bsub -q ${processing_q} -a mympi -n ${nr_cpu} -o log_files/${job_name}_${detector}%J.log -J ${job_name}_${detector} mpi_fxs_launch -j ${job_type} -e ${experiment} -r $1 -a ${detector} -d ${distance} -t ${mask_thr} -f ${file_type} -w ${wavelength}  -x ${x} -y ${y} --dx ${dx} --dy ${dy} --index ${index} -m ${mask}

     # Launch Interactive job

     mpirun -n ${nr_cpu}   mpi_fxs_launch -j ${job_type} -e ${experiment} -r $1 -a ${detector} -d ${distance} -t ${mask_thr}  -f ${file_type} -w ${wavelength} -x ${x} -y ${y} --dx ${dx} --dy ${dy} --index ${index} --first ${first} --last ${last} -m ${mask}

     # Launch single CPU

     mpi_fxs_launch -j ${job_type} -e ${experiment} -r $1 -a ${detector} -d ${distance} -t ${mask_thr} -f ${file_type} -w ${wavelength}  -x ${x} -y ${y} --dx ${dx} --dy ${dy}  --first ${first} --last ${last} -m ${mask} --index ${index}

     ######################################################################

     Instructions for viewing intermediate results using psplot

     1. In terminal type: bjobs -w                              (to get the id of the first analysis node (psanaXXXX))
     2. In terminal type: psplot -s psanaXXXX IND AVE SAXS C2   (to select the type of results to view)

  """

  if argv == None:
     argv = sys.argv[1:]

  try:
     from mpi4py import MPI
  except ImportError:
     raise Sorry("MPI not found")

  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  size = comm.Get_size()

  command_line = (libtbx.option_parser.option_parser(
    usage="""

%s    -j job type -e  experiment  -r  run  -a  address -d  detector distance  [-o  outputdir]
               [-I first frame] [-F last frame] [-H  hit intensity]   [-m  mask path]
               [-B background image path]  [-M  background mask path]  [-i  index path] [-G geometry path]
               [-u  pixel size]  [-w  wavelength override]  [-t  threshold for masking]  [-q  nr of q bins]
               [-p  nr of phi bins]  [-Q  sampled stepsize q]  [-P  sampled stepsize phi]  [-x beam center x] [-y beam center y]
               [-z rmax for beam c]  [-Z  sample step size in r]  [-x  bound in x for beam c ]  [-y  bound in y for beam c],
               [-s  particle radius]  [-S  q-range for size]  [-g plot output] [-E path to external directory] [-f data file type]
""" % libtbx.env.dispatcher_name)
                .option(None, "--job", "-j",
                        type="string",
                        default=None,
                        dest="job",
                        help="job type: index | background |correlation |calibration |mask ")
                .option(None, "--experiment", "-e",
                        type="string",
                        default=None,
                        dest="experiment",
                        help="experiment name: amo86615 ")
                .option(None, "--run", "-r",
                        type="int",
                        default=None,
                        dest="run",
                        help="run number: eg 120")
                .option(None, "--address", "-a",
                        type="string",
                        default="None",
                        dest="address",
                        help="detector address: pnccdBack | pnccdFront ")
                .option(None, "--outputdir", "-o",
                        type="string",
                        default=None,
                        dest="outputdir",
                        metavar="PATH",
                        help="Optional path to output directory for output files")
                .option(None, "--first", "-I",
                        type="int",
                        default=None,
                        dest="first",
                        help="first frame to process")
                .option(None, "--last", "-F",
                        type="int",
                        default=None,
                        dest="last",
                        help="last frame to process")
                .option(None, "--hit", "-H",
                        type="float",
                        default=None,
                        dest="hit",
                        help="hit intensity")
                .option(None, "--mask path", "-m",
                        type="string",
                        default=None,
                        dest="mask_path",
                        metavar="PATH",
                        help="static mask file")
                .option(None, "--background image path", "-B",
                        type="string",
                        default=None,
                        dest="bg_img_path",
                        metavar="PATH",
                        help="background image file")
                .option(None, "--background mask path", "-M",
                        type="string",
                        default=None,
                        dest="bg_msk_path",
                        metavar="PATH",
                        help="background mask file")
                .option(None, "--index", "-i",
                        type="string",
                        default=None,
                        dest="param_path",
                        metavar="PATH",
                        help="index file")
                .option(None, "--geometry", "-G",
                        type="string",
                        default=None,
                        dest="geom_path",
                        metavar="PATH",
                        help="geometry file")
                .option(None, "--distance", "-d",
                        type="float",
                        default=None,
                        dest="det_distance",
                        help="detector distance [mm]")
                .option(None, "--pixel size", "-u",
                        type="float",
                        default=0.075,
                        dest="det_pixel",
                        help="pixels size [mm]")
                .option(None, "--wavelength", "-w",
                        type="float",
                        default=None,
                        dest="lambda_b",
                        help="wavelength overwrite [ang]")
                .option(None, "--mask threshold", "-t",
                        type="float",
                        default=None,
                        dest="thr",
                        help="masking threshold")
                .option(None, "--nq", "-q",
                        type="int",
                        default=None,
                        dest="nQ",
                        help="nr of q bins")
                .option(None, "--nPhi", "-p",
                        type="int",
                        default=None,
                        dest="nPhi",
                        help="nr of phi bins")
                .option(None, "--dq", "-Q",
                        type="int",
                        default=1,
                        dest="dQ",
                        help="sampled step size in q")
                .option(None, "--dphi", "-P",
                        type="int",
                        default=1,
                        dest="dP",
                        help="sampled step size in phi")
                .option(None, "--x", "-x",
                        type="int",
                        default=0,
                        dest="x",
                        help="beam center in x")
                .option(None, "--y", "-y",
                        type="int",
                        default=0,
                        dest="y",
                        help="beam center in y")
                .option(None, "--r max", "-z",
                        type="int",
                        default=None,
                        dest="r_max",
                        help="r_max for beam center refinement")
                .option(None, "--dr", "-Z",
                        type="int",
                        default=1,
                        dest="dr",
                        help="step size for beam center refinement")
                .option(None, "--dx", "-X",
                        type="int",
                        default=5,
                        dest="dx",
                        help="grid size in x for for beam center refinement, x+/-dx")
                .option(None, "--dy", "-Y",
                        type="int",
                        default=5,
                        dest="dy",
                        help="grid size in y for for beam center refinement, y+/-dy")
                .option(None, "--particle radius", "-s",
                        type="float",
                        default=None,
                        dest="r0",
                        help="initial particle radius [ang]")
                .option(None, "--q-range for sizing", "-S",
                        type="float",
                        default=None,
                        dest="q_bound",
                        help="q-max for sizing [ang^-1]")
                .option(None, "--plot", "-g",
                        type="int",
                        default=0,
                        dest="plot",
                        help="plot output [0/1]")
                .option(None, "--datadir", "-E",
                        type="string",
                        default=None,
                        dest="xtc_dir",
                        metavar="PATH",
                        help="Alternative path to xtc or h5 directory")
                .option(None, "--file type", "-f",
                        type="string",
                        default='idx',
                        dest="ftype",
                        help="Type of file (idx, smd, idx_ffb, smd_ffb, h5, xtc)")
                ).process(args=argv)

  # Check mandatory parameters
  if len(command_line.args) > 0 or \
     command_line.options.job is None or \
     command_line.options.experiment is None or \
     command_line.options.run is None or \
     command_line.options.address is None or \
     command_line.options.det_distance is None:
   command_line.parser.show_help()
   return

  # Create output directory
  if command_line.options.outputdir is None :
     if rank == 0 :
        cpath = os.getcwd()
        rpath = os.path.join(cpath,'run_' + str(command_line.options.run) +'/')
        if not os.path.exists(rpath) :
           os.mkdir(rpath)
        ppath = os.path.join(rpath, str(command_line.options.address) + '/' )
        if not os.path.exists(ppath) :
           os.mkdir(ppath)
        jpath = os.path.join(ppath, str(command_line.options.job) + '/' )
        if not os.path.exists(jpath) :
           os.mkdir(jpath)
        d=[]
        for name in glob.glob(jpath + '*'):
            d.append(name)
        tpath = os.path.join(jpath,str(len(d)+1))
        os.mkdir(tpath)

        command_line.options.outputdir = tpath


  cargs = command_line.options

  if command_line.options.job   == 'index' :
     mpi_fxs_index.compute_index(cargs)
  elif command_line.options.job == 'background' :
     mpi_fxs_bg.compute_bg(cargs)
  elif command_line.options.job == 'correlation' :
     mpi_fxs_c2.compute_c2(cargs)
  elif command_line.options.job == 'calibration' :
     mpi_fxs_calib.compute_calib(cargs)
  elif command_line.options.job == 'mask' :
     mpi_fxs_mask.compute_mask(cargs)
  else:
     print("*** No recognizable job type chose: index | background | correlation ***")
     command_line.parser.show_help()
     return


  # Move launch and log file to directory

  MPI.Finalize()

if (__name__ == "__main__"):
  sys.exit(launch(sys.argv[1:]))
