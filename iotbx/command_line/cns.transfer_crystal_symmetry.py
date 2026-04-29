"""Transfer crystal symmetry to a cns file"""
from __future__ import absolute_import, division, print_function
from iotbx import cns
import iotbx.cns.space_group_symbols
from iotbx import crystal_symmetry_from_any
from libtbx.str_utils import show_string
from libtbx.utils import Sorry, Usage, plural_s, detect_binary_file
import sys, os
from six.moves import zip

def run(args):
  if (len(args) != 2):
    raise Usage("""\
iotbx.cns.transfer_crystal_symmetry any_symmetry_source_file cns_input_file
  *********************************************
  NOTE: the cns_input_file is changed in place.
  *********************************************""")
  #
  for file_name in args:
    if (not os.path.exists(file_name)):
      raise Sorry("No such file: %s" % show_string(file_name))
  source, target = args
  crystal_symmetry = crystal_symmetry_from_any.extract_from(source)
  if (crystal_symmetry is None):
    raise Sorry(
      "Unknown file format or unit cell and/or space group"
      " missing from file: " + show_string(source))
  cns_space_group_symbol = cns.space_group_symbols.cns_format(
    space_group_info=crystal_symmetry.space_group_info())
  if (cns_space_group_symbol is None):
    raise Sorry("Space group not available in CNS: %s" %
      show_string(str(crystal_symmetry.space_group_info())))
  sg = '"%s"' % cns_space_group_symbol
  a,b,c,alpha,beta,gamma = ["%.6g" % p
    for p in crystal_symmetry.unit_cell().parameters()]
  parameter_names = ["sg", "a", "b", "c", "alpha", "beta", "gamma"]
  parameters_found = dict(zip(parameter_names, [0]*len(parameter_names)))
  parameters_changed = []
  lines_out = []
  detect_binary = detect_binary_file(monitor_initial=100)
  try: cns_inp = open(target).read().splitlines()
  except IOError as e:
    raise Sorry("Error reading file %s (%s)" % (show_string(target), str(e)))
  end_block_parameter_definition = False
  for line in cns_inp:
    if (detect_binary is not None):
      is_binary = detect_binary.is_binary_file(block=line)
      if (is_binary is not None):
        if (is_binary):
          raise Sorry("%s appears to be a binary file." % show_string(target))
        detect_binary = None
    if (end_block_parameter_definition):
      lines_out.append(line)
    else:
      l = line.strip().replace(" ","")
      if (l == "){-endblockparameterdefinition-}"):
        lines_out.append(line)
        end_block_parameter_definition = True
      else:
        line_out = line
        for p in parameter_names:
          if (l.startswith("{===>}%s=" % p) and l.endswith(";")):
            parameters_found[p] += 1
            line_out = '{===>} %s=%s;' % (p, vars()[p])
            if (line_out != line): parameters_changed.append(p)
            break
        lines_out.append(line_out)
  if (list(parameters_found.values()).count(1) != 7):
    raise Sorry("Unexpected set of variable names in %s:\n  counts: %s" % (
      show_string(target), str(parameters_found)))
  elif (len(parameters_changed) == 0):
    print("Info: no changes, %s was not modified." % show_string(target))
  else:
    string_out = "\n".join(lines_out)
    print("Info: %d change%s" % plural_s(len(parameters_changed)), \
      "(%s)," % ", ".join(parameters_changed), \
      "writing modified file %s." % show_string(target))
    try: print(string_out, file=open(target, "w"))
    except IOError as e:
      raise Sorry("Error writing file %s (%s)" % (show_string(target), str(e)))

if (__name__ == "__main__"):
  run(sys.argv[1:])
