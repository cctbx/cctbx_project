# LIBTBX_SET_DISPATCHER_NAME phenix.csv_as_mtz

import sys, os
from cctbx import miller
from iotbx.option_parser import iotbx_option_parser
from iotbx.pdb import crystal_symmetry_from_pdb
from libtbx.utils import Sorry
from cctbx.array_family import flex

def run(args, command_name = "mmtbx.csv_to_mtz"):
  if (len(args) == 0): args = ["--help"]
  try:
    command_line = (iotbx_option_parser(
      usage="%s [reflection_csv_file] [options]" % command_name,
      description='Example: %s 1m5u-sf.csv --use_model=1m5u.pdb'%command_name)
      .enable_symmetry_comprehensive()
      .option(None, "--output_file",
        action="store",
        default=False,
        type="string",
        help="Output mtz file name.")
      .option(None, "--use_model",
        action="store",
        default=False,
        type="string",
        help="Use PDB model to make better guess about reflection data type.")
    ).process(args=args)
  except Exception, e:
    if(str(e) != "0"): print str(e)
    sys.exit(0)
  crystal_symmetry = command_line.symmetry
  if(command_line.symmetry.unit_cell() is None or
     command_line.symmetry.space_group_info() is None):
    if(command_line.options.use_model):
      crystal_symmetry = crystal_symmetry_from_pdb.extract_from(
         file_name=command_line.options.use_model)
  if(crystal_symmetry.unit_cell() is None or
     crystal_symmetry.space_group_info() is None):
    raise Sorry(
      "Crystal symmetry is not defined. Please use the --symmetry option.\n"
      "Type %s without arguments to see more options."%command_name)
  if(len(command_line.args) > 1):
    print "%d arguments are given from the command line:"% \
      len(command_line.args), command_line.args
    raise Sorry("Please specify one reflection csv file.")
  file_name = command_line.args[0]
  if(not os.path.isfile(file_name)):
    raise Sorry("File is not found: %s"%file_name)
  data = flex.double()
  sigmas = flex.double()
  flags = flex.int()
  column_ids, columns = parse_csv_file(file_name=file_name)
  data_label_root = column_ids[3]
  ms = miller.set(crystal_symmetry, flex.miller_index(columns[0]))
  for d in columns[1]:
    data.append(float(d))
  for sig in columns[2]:
    sigmas.append(float(sig))
  for flag in columns[3]:
    flags.append(int(flag))
  assert len(data) == len(sigmas)
  assert len(flags) == len(data)
  ma = miller.array(ms, data, sigmas)
  if data_label_root.startswith('F'):
    ma.set_observation_type_xray_amplitude()
  elif data_label_root.startswith('I'):
    ma.set_observation_type_xray_intensity()
  else:
    ma.set_observation_type_xray_amplitude()
  flags_ma = miller.set(
      crystal_symmetry = crystal_symmetry,
      indices          = ma.indices()).array(data = flags)
  mtz_dataset = ma.as_mtz_dataset(
    column_root_label = data_label_root)
  mtz_dataset.add_miller_array(
      miller_array      = flags_ma,
      column_root_label = "R-free-flags")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = command_line.options.output_file)

def parse_csv_file(file_name):
  f = open(file_name)
  column_ids = f.readline().strip().split(',')
  data = []
  for id in column_ids:
    data.append([])
  for line in f.readlines():
    temp = line.strip().split(',')
    for i, value in enumerate(temp):
      if i == 0:
        h = int(value)
      elif i == 1:
        k = int(value)
      elif i == 2:
        l = int(value)
      elif i == 3:
        hkl = (h, k, l)
        data[0].append(hkl)
        data[i-2].append(value)
      else:
        data[i-2].append(value)
  return column_ids, data


if(__name__ == "__main__"):
   run(sys.argv[1:])
