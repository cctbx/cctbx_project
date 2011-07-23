# LIBTBX_SET_DISPATCHER_NAME phenix.secondary_structure_restraints

from mmtbx.secondary_structure import *
from mmtbx.geometry_restraints import hbond
import iotbx.pdb
import iotbx.phil
from scitbx.array_family import flex
from libtbx.utils import Sorry
import cStringIO
import sys

def run (args, out=sys.stdout, log=sys.stderr) :
  pdb_files = []
  sources = []
  force_new_annotation = False
  master_phil_str = """
show_all_params = False
  .type = bool
filter_outliers = True
  .type = bool
format = *phenix phenix_bonds pymol refmac kinemage
  .type = choice
quiet = False
  .type = bool
refinement {
  secondary_structure {
    %s
  }
}""" % sec_str_master_phil_str
  pcl = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_phil_str,
    pdb_file_def="refinement.secondary_structure.input.file_name")
  working_phil = pcl.work
  params = working_phil.extract()
  ss_params = params.refinement.secondary_structure
  if params.quiet :
    out = cStringIO.StringIO()
  if len(pdb_files) > 0 :
    ss_params.input.file_name.extend(pdb_files)
  pdb_files = ss_params.input.file_name
  if len(pdb_files) == 0 :
    raise Sorry("No PDB files specified.")
  pdb_combined = iotbx.pdb.combine_unique_pdb_files(file_names=pdb_files)
  pdb_structure = iotbx.pdb.input(source_info=None,
    lines=flex.std_string(pdb_combined.raw_records))
  pdb_hierarchy = pdb_structure.construct_hierarchy()
  pdb_hierarchy.atoms().reset_i_seq()
  xray_structure = pdb_structure.xray_structure_simple()
  if len(pdb_hierarchy.models()) != 1 :
    raise Sorry("Multiple models not supported.")
  m = manager(pdb_hierarchy=pdb_hierarchy,
    xray_structure=xray_structure,
    sec_str_from_pdb_file=pdb_structure.extract_secondary_structure(),
    params=ss_params)
  m.find_automatically(log=log)
  prefix_scope="refinement.secondary_structure"
  if (params.format != "phenix") :
    prefix_scope = ""
  ss_phil = None
  working_phil = m.as_phil_str(master_phil=sec_str_master_phil)
  phil_diff = sec_str_master_phil.fetch_diff(source=working_phil)
  if params.format == "phenix_bonds" :
    raise Sorry("Not yet implemented.")
  elif params.format in ["pymol", "refmac", "kinemage"] :
    m.show_summary(out=log)
    build_proxies = m.create_hbond_proxies(
      log=log,
      restraint_type=None,
      as_python_objects=True)
    if (len(build_proxies.proxies) == 0) :
      pass
    elif params.format == "pymol" :
      hbond.as_pymol_dashes(
        proxies=build_proxies.proxies,
        pdb_hierarchy=pdb_hierarchy,
        filter=params.filter_outliers,
        out=out)
    elif params.format == "kinemage" :
      hbond.as_kinemage(
        proxies=build_proxies.proxies,
        pdb_hierarchy=pdb_hierarchy,
        filter=params.filter_outliers,
        out=out)
    else :
      hbond.as_refmac_restraints(
        proxies=build_proxies.proxies,
        pdb_hierarchy=pdb_hierarchy,
        filter=params.filter_outliers,
        out=out)
  else :
    #working_phil.show(out=out)
    print >> out, "# These parameters are suitable for use in phenix.refine."
    if (prefix_scope != "") :
      print >> out, "%s {" % prefix_scope
    if params.show_all_params :
      working_phil.show(prefix="  ", out=out)
    else :
      phil_diff.show(prefix="  ", out=out)
    if (prefix_scope != "") :
      print >> out, "}"
    #print >> out, ss_params_str
    return working_phil.as_str()

if __name__ == "__main__" :
  run(sys.argv[1:])
