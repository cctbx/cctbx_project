from cctbx import sgtbx
from scitbx.python_utils.command_line import parse_options
import sys, os, random

def get_test_space_group_symbols(flag_AllSpaceGroups, flag_AllSettings):
  if (flag_AllSettings):
    return [symbols.extended_hermann_mauguin()
            for symbols in sgtbx.space_group_symbol_iterator()]
  if (flag_AllSpaceGroups):
    sg_numbers = xrange(1, 231)
  else:
    sg_numbers = (1,2,3,15,16,74,75,76,142,143,144,157,167,168,194,195,230)
  return [sgtbx.space_group_symbols(n).extended_hermann_mauguin()
          for n in sg_numbers]

def report_cpu_times():
  t = os.times()
  print "u+s,u,s: %.2f %.2f %.2f" % (t[0] + t[1], t[0], t[1])

def loop_space_groups(argv, flags, call_back, symbols_to_stdout=0):
  if (not flags.RandomSeed): random.seed(0)
  if (len(argv) > 0 + flags.n):
    symbols = argv
  else:
    symbols = get_test_space_group_symbols(
      flags.AllSpaceGroups,
      flags.AllSettings)
  for symbol in symbols:
    if (symbol.startswith("--")): continue
    space_group_info = sgtbx.space_group_info(symbol)
    sys.stdout.flush()
    print >> sys.stderr, space_group_info
    sys.stderr.flush()
    if (symbols_to_stdout):
      print space_group_info
      sys.stdout.flush()
    continue_flag = call_back(flags, space_group_info)
    sys.stdout.flush()
    if (continue_flag == 00000): break
  report_cpu_times()

def parse_options_loop_space_groups(argv, call_back,
                                    keywords=(),
                                    symbols_to_stdout=0):
  flags = parse_options(argv, (
    "Verbose",
    "RandomSeed",
    "AllSpaceGroups",
    "AllSettings") + keywords
  )
  loop_space_groups(argv, flags, call_back, symbols_to_stdout)
