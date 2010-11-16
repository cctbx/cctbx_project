from libtbx.utils import Sorry
from mmtbx.scaling import sad_scale, sir_scale, rip_scale, siras_scale
from mmtbx.scaling import twmad_scale
import sys


def print_banner(command_name):
  hashes   = "#########################################"
  subtitle = "####    Delta F and FA estimation    ####"
  n = len(hashes) - 10 - len(command_name)
  left = " " * (n//2)
  right = left + " " * (n%2)
  print hashes
  print "#### %s%s%s ####" % (left, command_name, right)
  print subtitle
  print hashes

def run(args, command_name="phenix.fest"):
  print_banner(command_name=command_name)

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
  if len(args)==0:
    print
    print "usage: %s <EXPERIMENT TYPE> <FLAGS and/or PARAMETER FILE>" \
      % command_name
    print
    raise Sorry("No instructions received")

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
