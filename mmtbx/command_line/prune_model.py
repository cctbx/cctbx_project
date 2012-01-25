
import libtbx.phil
from libtbx.str_utils import make_header
from libtbx.utils import Usage, multi_out
import os
import sys

master_phil = libtbx.phil.parse("""
include scope mmtbx.utils.cmdline_input_phil_str
prune {
  resolution_factor = 1/4.
    .type = float
  sidechains = True
    .type = bool
  residues = True
    .type = bool
  min_c_alpha_2fofc = 1.0
    .type = float
  min_c_alpha_fofc = -3.0
    .type = float
  min_sidechain_2fofc = 0.5
    .type = float
  max_sidechain_fofc = 2.8
    .type = float
  min_cc = 0.75
    .type = float
  min_cc_sidechain = 0.7
    .type = float
}
output {
  file_name = None
    .type = path
}
""", process_includes=True)

def id_str (chain, residue_group, atom_group) :
  return "%3s %s%4s%s" % (atom_group.resname, chain.id, residue_group.resseq,
    residue_group.icode)

def prune_model (
    f_map_coeffs,
    diff_map_coeffs,
    model_map_coeffs,
    pdb_hierarchy,
    params,
    out=None) :
  if (out is None) :
    out = sys.stdout
  assert (len(pdb_hierarchy.models()) == 1)
  from mmtbx.real_space_correlation import set_details_level_and_radius
  from cctbx import maptbx
  from scitbx.array_family import flex
  unit_cell = f_map_coeffs.unit_cell()
  f_map_fft = f_map_coeffs.fft_map(resolution_factor=params.resolution_factor)
  f_map = f_map_fft.apply_sigma_scaling().real_map()
  diff_map_fft = diff_map_coeffs.fft_map(
    resolution_factor=params.resolution_factor)
  diff_map = diff_map_fft.apply_sigma_scaling().real_map()
  model_map_fft = model_map_coeffs.fft_map(
    resolution_factor=params.resolution_factor)
  model_map = model_map_fft.apply_sigma_scaling().real_map()
  atom_detail, residue_detail, atom_radius = set_details_level_and_radius(
    details_level="automatic",
    d_min=f_map_coeffs.d_min(),
    atom_radius=None)
  n_res_removed = 0
  n_sc_removed = 0
  make_header("Pruning residues and sidechains", out=out)
  for chain in pdb_hierarchy.models()[0].chains() :
    residue_id_hash = {}
    removed_resseqs = []
    if (len(chain.conformers()) > 1) :
      print >> out, "WARNING: chain '%s' has multiple conformers" % chain.id
    first_conf = chain.conformers()[0]
    if (not first_conf.is_protein()) :
      continue
    for j_seq, residue_group in enumerate(chain.residue_groups()) :
      residue_id_hash[residue_group.resid()] = j_seq
      for atom_group in residue_group.atom_groups() :
        ag_id_str = id_str(chain, residue_group, atom_group)
        remove_atom_group = False
        sidechain_atoms = []
        for atom in atom_group.atoms() :
          if (atom.name == " CA ") :
            if (not params.residues) : # sidechains only
              continue
            # first check the absolute map values at the C-alpha position
            site_frac = unit_cell.fractionalize(atom.xyz)
            map_value = f_map.tricubic_interpolation(site_frac)
            diff_map_value = diff_map.tricubic_interpolation(site_frac)
            if (map_value < params.min_c_alpha_2fofc) :
              print >> out, "%s: C-alpha 2mFo-DFc = %5.2f" % (ag_id_str,
                map_value)
              remove_atom_group = True
            elif (diff_map_value < params.min_c_alpha_fofc) :
              print >> out, "%s: C-alpha  mFo-DFc = %5.2f" % (ag_id_str,
                diff_map_value)
              remove_atom_group = True
            if (remove_atom_group) :
              break
          elif (not atom.name.strip() in ["N", "O", "C", "H", "CA", "CB"]) :
            sidechain_atoms.append(atom)
        if (not remove_atom_group) :
          # map values look okay - now check overall CC
          atoms = atom_group.atoms()
          sites_cart = atoms.extract_xyz()
          radii = flex.double([atom_radius] * sites_cart.size())
          sel = maptbx.grid_indices_around_sites(
            unit_cell=unit_cell,
            fft_n_real=f_map.focus(),
            fft_m_real=f_map.all(),
            sites_cart=sites_cart,
            site_radii=get_atom_radii(atoms, atom_radius))
          f_map_sel = f_map.select(sel)
          model_map_sel = model_map.select(sel)
          cc = flex.linear_correlation(x=f_map_sel,
            y=model_map_sel).coefficient()
          if (cc < params.min_cc) and (params.residues) :
            print >> out, "%s: overall CC = %4.2f" % (ag_id_str, cc)
            remove_atom_group = True
          elif (len(sidechain_atoms) > 0) and (params.sidechains) :
            # overall CC is acceptable - now look at sidechain alone
            remove_sidechain = False
            sites_cart = flex.vec3_double()
            sites_cart_nonH = flex.vec3_double()
            for atom in sidechain_atoms :
              sites_cart.append(atom.xyz)
              if (not atom.element.strip() in ["H","D"]) :
                sites_cart_nonH.append(atom.xyz)
            if (len(sites_cart_nonH) == 0) :
              continue # ALA, GLY?
            sel = maptbx.grid_indices_around_sites(
              unit_cell=unit_cell,
              fft_n_real=f_map.focus(),
              fft_m_real=f_map.all(),
              sites_cart=sites_cart,
              site_radii=get_atom_radii(sidechain_atoms, atom_radius))
            f_map_sel = f_map.select(sel)
            model_map_sel = model_map.select(sel)
            cc = flex.linear_correlation(x=f_map_sel,
              y=model_map_sel).coefficient()
            if (cc < params.min_cc_sidechain) :
              print >> out, "%s: sidechain CC = %4.2f" % (ag_id_str, cc)
              remove_sidechain = True
            else :
              # sidechain CC is okay - finally, check mean map values
              sel = maptbx.grid_indices_around_sites(
                unit_cell=unit_cell,
                fft_n_real=f_map.focus(),
                fft_m_real=f_map.all(),
                sites_cart=sites_cart_nonH,
                site_radii=flex.double([atom_radius] * len(sites_cart_nonH)))
              f_map_sel = f_map.select(sel)
              diff_map_sel = diff_map.select(sel)
              mean_f_value = flex.mean(f_map_sel.as_1d())
              mean_diff_value = flex.mean(diff_map_sel.as_1d())
              if (mean_f_value < params.min_sidechain_2fofc) :
                print >> out, "%s: sidechain 2mFo-DFc = %4.2f" % \
                  (ag_id_str, mean_f_value)
                remove_sidechain = True
              elif (mean_diff_value > params.max_sidechain_fofc) :
                print >> out, "%s: sidechain  mFo-DFc = %4.2f" % \
                  (ag_id_str, mean_diff_value)
                remove_sidechain = True
            if (remove_sidechain) :
              assert (params.sidechains)
              for atom in sidechain_atoms :
                atom_group.remove_atom(atom)
              n_sc_removed += 1
        if (remove_atom_group) :
          assert (params.residues)
          residue_group.remove_atom_group(atom_group)
      if (len(residue_group.atom_groups()) == 0) :
        chain.remove_residue_group(residue_group)
        n_res_removed += 1
        removed_resseqs.append(residue_group.resseq_as_int())
    # Final pass: remove lone single/pair residues
    if (params.residues) :
      for residue_group in chain.residue_groups() :
        if (residue_group.icode.strip() != "") :
          continue
        resseq = residue_group.resseq_as_int()
        if ((resseq - 1 in removed_resseqs) and
            ((resseq + 1 in removed_resseqs) or
             (resseq + 2 in removed_resseqs))) :
          print >> out, "Residue %s %s is not part of a continuous chain" % \
            (chain.id, residue_group.resseq)
          chain.remove_residue_group(residue_group)
          removed_resseqs.append(resseq)
  print >> out, "Removed %d residues and %d sidechains" % (n_res_removed,
    n_sc_removed)
  return pdb_hierarchy

def get_atom_radii (atoms, atom_radius) :
  from scitbx.array_family import flex
  radii = flex.double([atom_radius] * len(atoms))
  for i_seq, atom in enumerate(atoms) :
    if (atom.element.strip().upper() in ["H", "D"]) :
      radii[i_seq] = 1.0
  return radii

def run (args, out=None) :
  if (out is None) : out = sys.stdout
  if (len(args) == 0) :
    raise Usage("""\
mmtbx.prune_model model.pdb data.mtz [options...]

Filters protein residues based on CC to 2mFo-DFc map and absolute
(sigma-scaled) values in 2mFo-DFc and mFo-DFc maps.  For fast automatic
correction of MR solutions after initial refinement (ideally with rotamer
correction) to remove spurious loops and sidechains.
""")
  from mmtbx.utils import cmdline_load_pdb_and_data
  import iotbx.pdb
  cmdline = cmdline_load_pdb_and_data(
    args=args,
    master_phil=master_phil,
    out=out,
    process_pdb_file=False,
    create_fmodel=True)
  params = cmdline.params
  fmodel = cmdline.fmodel
  if (params.output.file_name is None) :
    base_name = os.path.basename(params.input.pdb.file_name[0])
    params.output.file_name = os.path.splitext(base_name)[0] + "_pruned.pdb"
  log_file = os.path.splitext(os.path.basename(params.output.file_name))[0] + \
    ".log"
  log = open(log_file, "w")
  out2 = multi_out()
  out2.register("out", out)
  out2.register("log", log)
  map_helper = fmodel.electron_density_map()
  f_map_coeffs = map_helper.map_coefficients(map_type="2mFo-DFc")
  diff_map_coeffs = map_helper.map_coefficients(map_type="mFo-DFc")
  model_map_coeffs = map_helper.map_coefficients(map_type="Fc")
  new_hierarchy = prune_model(
    f_map_coeffs=f_map_coeffs,
    diff_map_coeffs=diff_map_coeffs,
    model_map_coeffs=model_map_coeffs,
    pdb_hierarchy=cmdline.pdb_hierarchy,
    params=params.prune,
    out=out2)
  f = open(params.output.file_name, "w")
  cryst1 = iotbx.pdb.format_cryst1_record(fmodel.xray_structure)
  f.write("REMARK edited by mmtbx.prune_model\n")
  f.write("%s\n" % cryst1)
  f.write(new_hierarchy.as_pdb_string())
  f.close()
  log.close()
  print >> out, "Wrote %s" % params.output.file_name

if (__name__ == "__main__") :
  run(sys.argv[1:])
