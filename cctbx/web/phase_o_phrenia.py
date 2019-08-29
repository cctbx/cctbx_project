from __future__ import absolute_import, division, print_function
from cctbx.examples import phase_o_phrenia
from iotbx.cns import sdb_reader
from cctbx import uctbx
from cctbx import sgtbx
from cctbx.array_family import flex
from cctbx.web import io_utils
from cctbx.web import cgi_utils
from six.moves import range

def interpret_form_data(form):
  inp = cgi_utils.inp_from_form(form,
    (("ucparams", None),
     ("sgsymbol", None),
     ("convention", ""),
     ("format", None),
     ("coor_type", None),
     ("coor_file", None),
     ("coordinates", None),
     ("min_distance_sym_equiv", "0.5"),
     ("d_min", "3"),
     ("min_peak_distance", "4.0"),
     ("max_reduced_peaks", "60")))
  inp.coordinates = cgi_utils.coordinates_from_form(form)
  return inp

def sdb_files_as_xray_structures(lines, unit_cell, space_group_info, min_distance_sym_equiv):
  sdb_files = sdb_reader.multi_sdb_parser(lines)
  xray_structures = []
  for sdb in sdb_files:
    if (unit_cell is not None): sdb.unit_cell = unit_cell
    if (space_group_info is not None): sdb.space_group_info = space_group_info
    xray_structure = sdb.as_xray_structure(min_distance_sym_equiv=min_distance_sym_equiv)
    xray_structures.append(xray_structure)
  return xray_structures

def inp_as_xray_structures(inp):
  unit_cell = inp.ucparams
  if (unit_cell is not None): unit_cell = uctbx.unit_cell(unit_cell)
  space_group_info = inp.sgsymbol
  if (space_group_info is not None):
    space_group_info = sgtbx.space_group_info(symbol=space_group_info, table_id=inp.convention)
  min_distance_sym_equiv=float(inp.min_distance_sym_equiv)
  return sdb_files_as_xray_structures(
    inp.coordinates,
    unit_cell,
    space_group_info,
    min_distance_sym_equiv)

def run(server_info, inp, status):
  print("<pre>")

  if (inp.format == "cns_sdb"):
    print("Minimum distance between symmetrically equivalent sites:", end=' ')
    print(float(inp.min_distance_sym_equiv))
    print()
    structures = inp_as_xray_structures(inp)
    if (len(structures) == 0):
      print("No CNS sdb files found!")
      print()
      print("Note that each file must start with {+ file: some_file_name +}")
      print("in order to be recognized.")
      print()
  else:
    if (inp.ucparams is None): inp.ucparams = ""
    if (inp.sgsymbol is None): inp.sgsymbol = "P1"
    special_position_settings = io_utils.special_position_settings_from_inp(inp)
    special_position_settings.show_summary()
    print("Minimum distance between symmetrically equivalent sites:", end=' ')
    print(special_position_settings.min_distance_sym_equiv())
    print()
    structures = [io_utils.structure_from_inp(inp, status, special_position_settings)]

  d_min = float(inp.d_min)
  print("Minimum d-spacing:", d_min)
  if (d_min <= 0.):
    raise ValueError("d-spacing must be greater than zero.")
  print()

  min_peak_distance = float(inp.min_peak_distance)
  print("Minimum peak distance:", min_peak_distance)
  if (min_peak_distance <= 0.):
    raise ValueError("min_peak_distance must be greater than zero.")
  print()

  max_reduced_peaks = int(inp.max_reduced_peaks)
  print("Maximum number of peaks:", max_reduced_peaks)
  if (max_reduced_peaks <= 0):
    raise ValueError("max_reduced_peaks must be greater than zero.")
  print()

  for structure in structures:
    if (inp.format == "cns_sdb"):
      structure.show_summary().show_scatterers()
      print()
    if (structure.scatterers().size() == 0): continue
    reduced_peaks = phase_o_phrenia.calculate_exp_i_two_phi_peaks(
      xray_structure=structure,
      d_min=d_min,
      min_peak_distance=min_peak_distance,
      max_reduced_peaks=max_reduced_peaks)

    print("Actual number of peaks:", len(reduced_peaks))
    print()

    plot_nx = min(len(reduced_peaks), 60)
    if (plot_nx > 0):
      plot_ny = max(10, plot_nx//3)
      if (plot_nx != max_reduced_peaks):
        print("Number of peaks used for plot:", plot_nx)
        print()
      print("Plot of relative peak heights:")
      print()
      plot = flex.bool(flex.grid(plot_nx, plot_ny))
      for i in range(plot_nx):
        height = reduced_peaks[i].height
        h = int(round(height * plot_ny))
        h = max(0, min(plot_ny, h))
        for j in range(h): plot[(i,j)] = True
      for j in range(plot_ny-1,-1,-1):
        line = ""
        for i in range(plot_nx):
          if (plot[(i,j)]): line += "*"
          else:                  line += " "
        print("    |" + line.rstrip())
      print("    -" + "-" * plot_nx)
      print()

      print("Peak list:")
      print("  Relative")
      print("   height   Fractional coordinates")
      for peak in reduced_peaks:
        print("    %5.1f" % (peak.height*100), " %8.5f %8.5f %8.5f" % peak.site)
      print()

  print("</pre>")
