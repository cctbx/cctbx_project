def run(args):
  import iotbx.phil
  pcl = iotbx.phil.process_command_line(args=args, master_string="""\
symmetry_search {
  structure_factor_d_min = 2
    .type = float
}
""")
  if (len(pcl.remaining_args) != 1):
    from libtbx.utils import Usage
    import libtbx.load_env
    print
    pcl.master.show()
    print
    raise Usage("%s pdb_file [parameters]" % libtbx.env.dispatcher_name)
  pdb_file = pcl.remaining_args[0]
  from libtbx.str_utils import show_string
  print "PDB file:", show_string(pdb_file)
  print
  pcl.work.show()
  print
  params = pcl.work.extract().symmetry_search
  import iotbx.pdb
  pdb_inp = iotbx.pdb.input(file_name=pdb_file)
  xs = pdb_inp.xray_structure_simple()
  print "Unit cell:", xs.unit_cell()
  print
  fc_p1 = xs.structure_factors(
    d_min=params.structure_factor_d_min).f_calc().expand_to_p1()
  from cctbx import symmetry_search
  sf_symm = symmetry_search.structure_factor_symmetry(f_in_p1=fc_p1)
  print sf_symm
  print
  print sf_symm.space_group_info
  print

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
