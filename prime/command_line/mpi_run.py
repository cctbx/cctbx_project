from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime.mpi_run

from prime.command_line.solve_indexing_ambiguity import indexing_ambiguity_handler
import sys, os
from subprocess import call

if __name__ == "__main__":
  argv = sys.argv[1:] if len(sys.argv) > 1 else None
  # determine indexing ambiguity and setup iparams
  txt_indexing_ambiguity = "Determine if there is an indexing ambiguity on the dataset"
  print txt_indexing_ambiguity
  idah = indexing_ambiguity_handler()
  sol_fname, iparams = idah.run(argv)
  if sol_fname is None:
    print "No ambiguity."
    txt_indexing_ambiguity += "\nNo ambiguity."
  else:
    print "Ambiguity is solved. Solution file was saved to :"+str(sol_fname)
    txt_indexing_ambiguity += "Ambiguity is solved. Solution file was saved to :"+str(sol_fname)
    iparams.indexing_ambiguity.index_basis_in = sol_fname
  # setup parameters
  iparams.flag_volume_correction = False
  if iparams.partiality_model == "Lognormal":
    iparams.voigt_nu = 0.008 #use voigt_nu as lognpdf zero parameter
  # scaling
  qCmd = 'mpirun prime.mpi_scale '+' '.join(argv)
  cmd = ['bsub','-n', str(iparams.n_processors),'-q', iparams.queue.qname, '-o', os.path.join(iparams.run_no,'log_doMerge.out'), qCmd]
call(cmd)
