# LIBTBX_SET_DISPATCHER_NAME phenix.suitename
# LIBTBX_SET_DISPATCHER_NAME molprobity.suitename

#import libtbx.load_env
#from libtbx.utils import Usage
from iotbx.cli_parser import CCTBXParser
from libtbx.program_template import ProgramTemplate
from libtbx.utils import multi_out, show_total_time

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

from dualparse import parseArgs, parseArgs1   # must come before suitename
from suitename import main  # from suitename.py in PARENT directory


def run(args):
  # create parser
  logger = multi_out()
  logger.register('stderr', sys.stderr)
  logger2 = multi_out()
  logger2.register('stdout', sys.stdout)

  working_phil = parseArgs(Program, logger)
  working_phil.show()


class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  program_name = "suitename"
  description="""
   < insert help text here>
""" % locals()

  # The substructure below has been commented out because it would make it more
  # difficult to integrate with suitename's ability to run independently
  # of CCTBX
  master_phil_str = """
    suitename {
      # output {
        string=False  
          .type=bool
          .help="output in string format, 3 characters per suite"
        kinemage=False
          .type=bool
          .help="output in kinemage format, useful for visualization"
        chart=False
          .type=bool
          .help="modifier to standard report, output without statistical summary"
        nosequence = False
          .type=bool
          .help="modifier to string format, do not include base letters"
        causes=False
          .type=bool
          .help="output extra details concerning the causes of each assignment made"
        test=False
          .type=bool
          .help="display a lat of additional information about program internals"
        # }
      # compute {
        satellites=False
          .type=bool
          .help="use the special satelliteWidths values for satellites" 
        nowannabe=False
          .type=bool
          .help="do not consider 'wannabe' clusters"
        noinc=False
          .type=bool
          .help="do not display incomplete suites"
        etatheta=False
          .type=bool
        altid="A"
          .type=str
          .help="which alternate conformer to use (A, B, etc)"
        altidfield = 3
          .type=int
          .help="which field gives the alternate conformer code"        
        # }
      # input {
        anglefields = 9
          .type=int
          .help="number of angle fields provided, for textual input only"
        pointidfields = 7
          .type=int
          .help="number of point id fields before the angle fields"
        # }
      }
"""
  datatypes = ['model', 'phil']  # also 
  data_manager_options = ['model_skip_expand_with_mtrix']
  known_article_ids = ['molprobity']

  def validate(self):
    pass

  def run(self):
      pass

  # end of class Program


# def rest_of_old_run(args):  -- an intermediate state, outdated technique
#   namespace, others = parser.parse_args(sys.argv[1:])
#   # whatever the old fashioned parser won't use becomes part of <others>
#   # and will be given to the phil parser
# 
#   parser = CCTBXParser(
#     program_class=Program,
#     logger=logger)
#   namespace = parser.parse_args(others)
# 
#   # start program
#   print('Starting job', file=logger)
#   print('='*79, file=logger)
#   phil1 = parser.working_phil.extract()
#   args2 = parseCommandLine()
#   task = Program(
#     parser.data_manager, phil1, logger=logger2)  
#   main()

#=============================================================================
# def run(args):
# 
#   # create parser
#   logger = multi_out()
#   logger.register('stderr', sys.stderr)
#   logger2 = multi_out()
#   logger2.register('stdout', sys.stdout)
# 
#   parser = CCTBXParser(
#     program_class=cablam.Program,
#     logger=logger)
#   namespace = parser.parse_args(sys.argv[1:])
# 
#   # start program
#   print('Starting job', file=logger)
#   print('='*79, file=logger)
#   task = cablam.Program(
#     parser.data_manager, parser.working_phil.extract(), logger=logger2)
# 
#   # validate inputs
#   task.validate()
# 
#   # run program
#   task.run()
# 
#   # stop timer
#   print('', file=logger)
#   print('='*79, file=logger)
#   print('Job complete', file=logger)
#   show_total_time(out=logger)

# =============================================================================
if __name__ == '__main__':
  run(sys.argv[1:])
