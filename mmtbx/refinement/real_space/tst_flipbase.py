from __future__ import division
from scitbx.array_family import flex
from libtbx import group_args
from libtbx.utils import user_plus_sys_time
from mmtbx.refinement.real_space import individual_sites, flipbase
import mmtbx

pdb_str_answer = """\
CRYST1   23.136   23.980   28.180  90.00  90.00  90.00 P 1
ATOM   1240  P    DG B  40       0.649   2.161 -36.081  0.50 17.21           P
ATOM   1241  OP1  DG B  40       1.754   3.055 -35.664  0.50 19.03           O
ATOM   1242  OP2  DG B  40      -0.231   2.567 -37.203  0.50 13.93           O
ATOM   1243  O5'  DG B  40       1.316   0.741 -36.417  0.50 17.66           O
ATOM   1244  C5'  DG B  40       2.224   0.095 -35.507  0.50 17.29           C
ATOM   1245  C4'  DG B  40       2.673  -1.267 -36.015  0.50 18.14           C
ATOM   1246  O4'  DG B  40       1.657  -2.274 -35.767  0.50 18.54           O
ATOM   1247  C3'  DG B  40       2.981  -1.325 -37.503  0.50 18.42           C
ATOM   1248  O3'  DG B  40       4.381  -1.282 -37.640  0.50 18.86           O
ATOM   1249  C2'  DG B  40       2.435  -2.669 -37.974  0.50 19.35           C
ATOM   1250  C1'  DG B  40       1.386  -3.034 -36.928  0.50 20.27           C
ATOM   1251  N9   DG B  40       0.023  -2.677 -37.316  0.50 21.30           N
ATOM   1252  C8   DG B  40      -0.406  -1.445 -37.759  0.50 21.71           C
ATOM   1253  N7   DG B  40      -1.676  -1.380 -38.029  0.50 22.28           N
ATOM   1254  C5   DG B  40      -2.128  -2.655 -37.739  0.50 21.94           C
ATOM   1255  C6   DG B  40      -3.437  -3.178 -37.849  0.50 20.86           C
ATOM   1256  O6   DG B  40      -4.476  -2.610 -38.224  0.50 22.01           O
ATOM   1257  N1   DG B  40      -3.477  -4.511 -37.466  0.50 20.98           N
ATOM   1258  C2   DG B  40      -2.385  -5.232 -37.045  0.50 21.41           C
ATOM   1259  N2   DG B  40      -2.658  -6.493 -36.734  0.50 22.19           N
ATOM   1260  N3   DG B  40      -1.143  -4.774 -36.928  0.50 20.66           N
ATOM   1261  C4   DG B  40      -1.093  -3.469 -37.297  0.50 20.70           C
TER
END
"""

pdb_str_poor = """\
CRYST1   23.136   23.980   28.180  90.00  90.00  90.00 P 1
ATOM   1240  P    DG B  40       0.612   2.104 -36.129  0.50 17.21           P
ATOM   1241  OP1  DG B  40       1.718   3.006 -35.787  0.50 19.03           O
ATOM   1242  OP2  DG B  40      -0.331   2.457 -37.202  0.50 13.93           O
ATOM   1243  O5'  DG B  40       1.193   0.664 -36.476  0.50 17.66           O
ATOM   1244  C5'  DG B  40       2.139   0.076 -35.623  0.50 17.29           C
ATOM   1245  C4'  DG B  40       2.467  -1.338 -36.039  0.50 18.14           C
ATOM   1246  O4'  DG B  40       1.323  -2.180 -35.893  0.50 18.54           O
ATOM   1247  C3'  DG B  40       2.890  -1.504 -37.491  0.50 18.42           C
ATOM   1248  O3'  DG B  40       4.284  -1.545 -37.594  0.50 18.86           O
ATOM   1249  C2'  DG B  40       2.227  -2.815 -37.926  0.50 19.35           C
ATOM   1250  C1'  DG B  40       1.572  -3.325 -36.650  0.50 20.27           C
ATOM   1251  N9   DG B  40       0.318  -4.068 -36.831  0.50 21.30           N
ATOM   1252  C8   DG B  40       0.099  -5.410 -36.618  0.50 21.71           C
ATOM   1253  N7   DG B  40      -1.128  -5.780 -36.836  0.50 22.28           N
ATOM   1254  C5   DG B  40      -1.761  -4.618 -37.199  0.50 21.94           C
ATOM   1255  C6   DG B  40      -3.099  -4.399 -37.545  0.50 20.86           C
ATOM   1256  O6   DG B  40      -4.016  -5.216 -37.606  0.50 22.01           O
ATOM   1257  N1   DG B  40      -3.326  -3.075 -37.841  0.50 20.98           N
ATOM   1258  C2   DG B  40      -2.364  -2.095 -37.807  0.50 21.41           C
ATOM   1259  N2   DG B  40      -2.755  -0.882 -38.105  0.50 22.19           N
ATOM   1260  N3   DG B  40      -1.117  -2.287 -37.484  0.50 20.66           N
ATOM   1261  C4   DG B  40      -0.884  -3.558 -37.193  0.50 20.70           C
TER
END
"""

def get_pdb_inputs(pdb_str, restraints):
  ppf = mmtbx.utils.process_pdb_file_srv(log=False).process_pdb_files(
    raw_records=pdb_str.splitlines())[0]
  xrs = ppf.xray_structure(show_summary = False)
  restraints_manager=None
  if(restraints):
    restraints_manager = mmtbx.restraints.manager(
      geometry      = ppf.geometry_restraints_manager(show_energies = False),
      normalization = True)
  return group_args(
    ph  = ppf.all_chain_proxies.pdb_hierarchy,
    grm = restraints_manager,
    xrs = xrs)

def get_map(xrs):
  f_calc = xrs.structure_factors(d_min = 1.2).f_calc()
  fft_map = f_calc.fft_map(resolution_factor=0.25)
  fft_map.apply_sigma_scaling()
  return fft_map.real_map_unpadded()

def exercise():
  """
  Exercise refine "easy" with DNA/RNA.
  """
  pi_good = get_pdb_inputs(pdb_str=pdb_str_answer, restraints=False)
  map_data = get_map(xrs=pi_good.xrs)
  xrs_good = pi_good.xrs.deep_copy_scatterers()
  pi_good.ph.write_pdb_file(file_name="answer.pdb",
    crystal_symmetry=xrs_good.crystal_symmetry())
  #
  pi_poor = get_pdb_inputs(pdb_str=pdb_str_poor, restraints=True)
  pi_poor.ph.write_pdb_file(file_name="poor.pdb")
  xrs_poor = pi_poor.xrs.deep_copy_scatterers()
  #
  d = xrs_good.distances(other=xrs_poor)
  print d.min_max_mean().as_tuple()
  assert flex.max(d)>3
  assert flex.mean(d)>1

  chain = 'B'
  res_num = 40
  alt_loc = None
  ero = False
  n_refine_cycles = 3
  target_map = get_map(xrs=pi_good.xrs)
  sites_cart = xrs_poor.sites_cart()
  for ch in pi_poor.ph.chains():
    if ch.id.strip() != chain : continue
    print 'chain'
    for rg in ch.residue_groups():
      if rg.resseq_as_int() != res_num : continue
      print "res#"
      if rg.have_conformers() and not alt_loc :
        s = 'Specified residue has alternate conformations. Please specify '
        raise RuntimeError(s + 'alt_loc on the command line')
      for residue in rg.atom_groups():
        if alt_loc and alt_loc != residue.altloc.strip():
          continue
        flipbase.flip_base(residue, angle=180)

        sites_cart.set_selected(residue.atoms().extract_i_seq(),
          residue.atoms().extract_xyz())
        xray_structure = xrs_poor.replace_sites_cart(sites_cart)
        sele = residue.atoms().extract_i_seq()
        print 'real-space refinement BEGIN'.center(79,'*')
        for i in range(n_refine_cycles):
          print 'real-space refinement cycle %i...' % (i + 1)
          ero = individual_sites.easy(
            map_data                    = target_map,
            xray_structure              = xray_structure,
            pdb_hierarchy               = pi_poor.ph,
            geometry_restraints_manager = pi_poor.grm,
            selection                   = sele)
        print 'real-space refinement FINISHED'.center(79,'*')
        xrs_refined = ero.xray_structure
  if not ero : raise RuntimeError('Specified residue not found')
  d = xrs_good.distances(other=xrs_refined)
  print d.min_max_mean().as_tuple()
  assert flex.max(d)<0.07
  assert flex.mean(d)<0.03
  ero.pdb_hierarchy.write_pdb_file(file_name="refined.pdb",
    crystal_symmetry=xrs_good.crystal_symmetry())



if(__name__ == "__main__"):
  timer = user_plus_sys_time()
  exercise()
  print "Time: %6.2f" % timer.elapsed()
  print "OK"
