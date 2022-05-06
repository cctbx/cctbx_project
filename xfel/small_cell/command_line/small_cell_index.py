from __future__ import absolute_import, division, print_function
#-*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tabwidth: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cctbx.small_cell_index
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

import xfel.small_cell.small_cell
from xfel.small_cell.small_cell import small_cell_index
import libtbx.load_env
import libtbx.option_parser
import sys,os
from six.moves import zip

small_cell_phil_str = """
small_cell {
  powdercell = None
    .type=unit_cell
    .help = "Specify unit cell for powder rings"
  spacegroup = None
    .type=str
    .help = "Specify spacegroup for the unit cell"
  high_res_limit = 1.5
    .type=float
    .help= "Highest resolution limit to process"
  min_spots_to_integrate = 5
    .type=int
    .help= "At least this many spots needed to have been indexed to integrate the image"
  interspot_distance = 5
    .type=int
    .help= "Minimum distance in pixels between a prediction and a spotfinder spot to be accepted"
  faked_mosaicity = 0.005
    .type=float
    .help= "Non-experimentally determined mosaicity to use for each image"
  spot_connection_epsilon = 2.e-3
    .type=float
    .help= "Epsilon for comparing measured vs. predicted inter-spot distances when building the maximum clique"
  d_ring_overlap_limit = 5
    .type = int
    .help = "Number of d rings a spot can overlap before it is removed from consideration. Set to None to use all spots, but this can be time consuming"
  override_wavelength = None
    .type=float
    .help = "Use to override the wavelength found in the image file"
  write_gnuplot_input = False
    .type = bool
    .help = "Use to produce a series of files as inputs to gnuplot to show the indexing results"
  max_calls_to_bronk = 100000
    .type = int
    .help = "Terminate indexing on this many calls to the maximum clique finder."
            "This eliminates a long tail of slow images with too many spots."
}
"""

dials_phil_str = """
include scope dials.algorithms.spot_finding.factory.phil_scope
"""

def run(argv=None):
  if (argv is None):
    argv = sys.argv

  from iotbx.phil import parse
  small_cell_phil = parse(small_cell_phil_str+dials_phil_str,process_includes=True)

  welcome_message = """
  %s [-s] -t PATH <directory or image paths>

  cctbx.small_cell: software for indexing sparse, still patterns.

  An excellent knowledge of the unit cell, detector distance, wavelength and
  beam center is required. Specify at least the unit cell in the target phil
  file passed in with the -t parameter.

  If the image can be integrated, the integrated intensities will be found in
  a *.int file (plain text) and in a cctbx.xfel integration pickle file.

  See Brewster, A.S., Sawaya, M.R., Rodriguez, J., Hattne, J., Echols, N.,
  McFarlane, H.T., Cascio, D., Adams, P.D., Eisenberg, D.S. & Sauter, N.K.
  (2015). Acta Cryst. D71, doi:10.1107/S1399004714026145.

  Showing phil parameters:

  """ % libtbx.env.dispatcher_name
  welcome_message += small_cell_phil.as_str(attributes_level = 2)

  command_line = (libtbx.option_parser.option_parser(
      usage=welcome_message)
                  .option(None, "--target", "-t",
                          type="string",
                          default=None,
                          dest="target",
                          metavar="PATH",
                          help="Target phil file")
                  .option(None, "--skip_processed_files", "-s",
                          action="store_true",
                          default=False,
                          dest="skip_processed_files",
                          help="Will skip images that have a .int file already created")
                  ).process(args=argv[1:])

  paths = command_line.args

  # Target phil file and at least one file to process are required
  if len(paths) == 0:
    command_line.parser.print_usage()
    return

  # Parse the target
  args = []
  if command_line.options.target is not None:
    args.append(parse(file_name=command_line.options.target,process_includes=True))

  horiz_phil = small_cell_phil.fetch(sources = args).extract()

  for path in paths:
    # process an entire directory
    if os.path.isdir(path):
      files = os.listdir(path)

      try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        # determine which subset of the files in this directory this process will
        # work on
        chunk = len(files) // size
        myfiles = files[rank*chunk:(rank+1)*chunk]
        if rank == 0:
          myfiles += files[len(files)-len(files)%size:len(files)]
      except ImportError as e:
        print("MPI not found, multiprocessing disabled")
        myfiles = files

      counts = []
      processed = []

      for file in myfiles:
        if (os.path.splitext(file)[1] == ".pickle" or os.path.splitext(file)[1] == ".edf") and os.path.basename(file)[0:3].lower() != "int" and file != "spotfinder.pickle":
          if command_line.options.skip_processed_files and os.path.exists(file + ".int"):
            print("Skiping %s as it has already been processed"%file)
            continue
          counts.append(small_cell_index(os.path.join(path,file),horiz_phil))
          if counts[-1] == None: counts[-1] = 0
          processed.append(file)
      for file, count in zip(processed,counts):
        print("%s %4d spots in max clique"%(file,count))

    # process a single file
    elif os.path.isfile(path):
      if os.path.splitext(path)[1] == ".txt":
        # Given a list of a file names in a text file, process each file listed
        f = open(path, "r")
        for line in f.readlines():
          if os.path.isfile(line.strip()):
            count = small_cell_index(line.strip(),horiz_phil)
            if count != None:
              print("%s %4d spots in max clique"%(line.strip(),count))
        f.close()
      elif os.path.splitext(path)[1] == ".int":
        # Summarize a .int file, providing completeness and multiplicity statistics
        f = open(path, "r")
        hkls_all = []
        hkls_unique = []
        files = []
        for line in f.readlines():
          strs = line.strip().split()
          src = strs[0].split(":")[0]
          if not src in files:
            files.append(src)
          hkl = (int(strs[7]), int(strs[8]), int(strs[9]))
          if not hkl in hkls_unique:
            hkls_unique.append(hkl)
          hkls_all.append(hkl)
        print("%d unique hkls from %d orginal files.  Completeness: "%(len(hkls_unique),len(files)))
        from cctbx.crystal import symmetry
        import cctbx.miller
        from cctbx.array_family import flex
        sym = symmetry(unit_cell=horiz_phil.small_cell.powdercell,
                       space_group_symbol=horiz_phil.small_cell.spacegroup)
        millerset = cctbx.miller.set(sym,flex.miller_index(hkls_unique),anomalous_flag=False)
        millerset = millerset.resolution_filter(d_min=horiz_phil.small_cell.high_res_limit)
        millerset.setup_binner(n_bins=10)
        data = millerset.completeness(True)
        data.show()
        data = millerset.completeness(False)
        print("Total completeness: %d%%\n"%(data * 100))

        print("%d measurements total from %d original files. Multiplicity (measurements/expected):"%(len(hkls_all),len(files)))
        millerset = cctbx.miller.set(sym,flex.miller_index(hkls_all),anomalous_flag=False)
        millerset = millerset.resolution_filter(d_min=horiz_phil.small_cell.high_res_limit)
        millerset.setup_binner(n_bins=10)
        data = millerset.completeness(True)
        data.show()
        print("Total multiplicty: %.3f"%(len(hkls_all)/len(millerset.complete_set().indices())))

        f.close()
      else:
        # process a regular image file
        count = small_cell_index(path,horiz_phil)
        if count != None:
          print("%s %4d spots in max clique"%(path,count))
    else:
      print("Not a file or directory: %s"%path)

  if xfel.small_cell.small_cell.app is not None:
    del xfel.small_cell.small_cell.app

if __name__=='__main__':
  sys.exit(run())
