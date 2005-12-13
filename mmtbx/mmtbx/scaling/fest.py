from cctbx import maptbx
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx import adptbx
import cctbx.sgtbx.lattice_symmetry
import cctbx.sgtbx.cosets
from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
import mmtbx.scaling
from mmtbx.scaling import absolute_scaling, relative_scaling
from mmtbx.scaling import matthews, twin_analyses
from mmtbx.scaling import basic_analyses, data_statistics,pair_analyses
import libtbx.phil.command_line
from cStringIO import StringIO
from scitbx.python_utils import easy_pickle
from mmtbx.scaling import sad_scale, sir_scale, rip_scale, siras_scale
from mmtbx.scaling import twmad_scale
import sys, os


from mmtbx.scaling import fa_estimation

def print_banner():
  print "#########################################"
  print "####              Fest               ####"
  print "####    Delta F and FA estimation    ####"
  print "#########################################"


def run(args):
  print_banner()

  scenarios = [ "SAD",#0

                "SIR",#1
                "RIP",#2

                "SIRAS",#3
                "RIPAS",#4

                "MIR",#5
                "MIRAS",#6

                "2WMAD",#7
                "3WMAD",#8
                "4WMAD",#9

                "2WMAD+NAT",#10
                "3WMAD+NAT",#11
                "4WMAD+NAT",#12

                "2WSAD",#13
                "3WSAD",#14
                "4WSAD" ]#15
  if args[0] in scenarios:
    print "experiment: " , args[0]
  else:
    print "Unknown experiment type"
    print " Choose from the following list:"
    for method in scenarios:
      print method

    raise Sorry("Unmknown experimenttype")

  ## this is a rather ugly enumeration of all cases
  if args[0]==scenarios[0]:
    sad_scale.run( args[1:] )

  elif args[0]==scenarios[1]:
    sir_scale.run( args[1:] )

  elif args[0]==scenarios[2]:
    rip_scale.run( args[1:] )

  elif args[0]==scenarios[3]:
    siras_scale.run( args[1:] )

  elif args[0]==scenarios[4]:
    siras_scale.run( args[1:] )

  elif args[0]==scenarios[7]:
    twmad_scale.run( args[1:] )
  
  else:
    print 
    print "Sorry, no time. This method has not yet been implemented."
    print "Currently, only SAD,SIR,RIP and 2WMAD are supported."
    print 

if (__name__ == "__main__"):
  run(sys.argv[1:])
