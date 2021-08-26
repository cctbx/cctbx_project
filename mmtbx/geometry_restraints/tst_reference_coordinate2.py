from __future__ import absolute_import, division, print_function

import iotbx.pdb
from cctbx.array_family import flex
from cctbx import adp_restraints # import dependency
import random
from libtbx import group_args
from mmtbx.refinement.real_space import individual_sites
from mmtbx.rotamer.rotamer_eval import RotamerEval
from cctbx import miller
from cctbx import maptbx
import mmtbx.model
from six.moves import zip

if(1):
  random.seed(0)
  flex.set_random_seed(0)

pdb_str = """\
CRYST1   16.960   19.455   19.841  90.00  90.00  90.00 P 1
SCALE1      0.058962  0.000000  0.000000        0.00000
SCALE2      0.000000  0.051401  0.000000        0.00000
SCALE3      0.000000  0.000000  0.050401        0.00000
ATOM      1  N   ASP A  18      14.998  12.244  11.645  1.00 10.00           N
ATOM      2  CA  ASP A  18      13.705  12.852  11.934  1.00 10.00           C
ATOM      3  C   ASP A  18      12.883  13.036  10.662  1.00 10.00           C
ATOM      4  O   ASP A  18      13.419  13.406   9.616  1.00 10.00           O
ATOM      5  CB  ASP A  18      13.893  14.198  12.639  1.00 10.00           C
ATOM      6  CG  ASP A  18      12.576  14.863  12.988  1.00 10.00           C
ATOM      7  OD1 ASP A  18      12.042  14.589  14.083  1.00 10.00           O
ATOM      8  OD2 ASP A  18      12.076  15.662  12.168  1.00 10.00           O
ATOM      9  N   ASN A  19      11.583  12.770  10.766  1.00 10.00           N
ATOM     10  CA  ASN A  19      10.649  12.917   9.651  1.00 10.00           C
ATOM     11  C   ASN A  19      11.021  12.082   8.426  1.00 10.00           C
ATOM     12  O   ASN A  19      11.667  12.571   7.500  1.00 10.00           O
ATOM     13  CB  ASN A  19      10.489  14.392   9.265  1.00 10.00           C
ATOM     14  CG  ASN A  19       9.422  14.608   8.210  1.00 10.00           C
ATOM     15  OD1 ASN A  19       8.499  13.805   8.068  1.00 10.00           O
ATOM     16  ND2 ASN A  19       9.542  15.699   7.461  1.00 10.00           N
ATOM     17  N   TYR A  20      10.605  10.819   8.431  1.00 10.00           N
ATOM     18  CA  TYR A  20      10.853   9.923   7.306  1.00 10.00           C
ATOM     19  C   TYR A  20       9.622   9.083   6.984  1.00 10.00           C
ATOM     20  O   TYR A  20       9.039   8.461   7.872  1.00 10.00           O
ATOM     21  CB  TYR A  20      12.051   9.014   7.591  1.00 10.00           C
ATOM     22  CG  TYR A  20      13.384   9.727   7.571  1.00 10.00           C
ATOM     23  CD1 TYR A  20      14.071   9.919   6.379  1.00 10.00           C
ATOM     24  CD2 TYR A  20      13.957  10.205   8.742  1.00 10.00           C
ATOM     25  CE1 TYR A  20      15.291  10.569   6.354  1.00 10.00           C
ATOM     26  CE2 TYR A  20      15.176  10.856   8.726  1.00 10.00           C
ATOM     27  CZ  TYR A  20      15.838  11.036   7.530  1.00 10.00           C
ATOM     28  OH  TYR A  20      17.052  11.683   7.510  1.00 10.00           O
ATOM     29  N   ARG A  21       9.247   9.073   5.704  1.00 10.00           N
ATOM     30  CA  ARG A  21       8.081   8.343   5.187  1.00 10.00           C
ATOM     31  C   ARG A  21       6.836   8.363   6.082  1.00 10.00           C
ATOM     32  O   ARG A  21       6.119   7.368   6.189  1.00 10.00           O
ATOM     33  CB  ARG A  21       8.450   6.904   4.781  1.00 10.00           C
ATOM     34  CG  ARG A  21       8.957   6.011   5.906  1.00 10.00           C
ATOM     35  CD  ARG A  21       9.297   4.620   5.397  1.00 10.00           C
ATOM     36  NE  ARG A  21       9.792   3.753   6.462  1.00 10.00           N
ATOM     37  CZ  ARG A  21       9.018   2.985   7.221  1.00 10.00           C
ATOM     38  NH1 ARG A  21       7.705   2.974   7.034  1.00 10.00           N
ATOM     39  NH2 ARG A  21       9.555   2.228   8.168  1.00 10.00           N
ATOM     40  N   GLY A  22       6.583   9.505   6.713  1.00 10.00           N
ATOM     41  CA  GLY A  22       5.417   9.664   7.563  1.00 10.00           C
ATOM     42  C   GLY A  22       5.733   9.547   9.041  1.00 10.00           C
ATOM     43  O   GLY A  22       5.411  10.440   9.825  1.00 10.00           O
ATOM     44  N   TYR A  23       6.364   8.441   9.423  1.00 10.00           N
ATOM     45  CA  TYR A  23       6.697   8.193  10.821  1.00 10.00           C
ATOM     46  C   TYR A  23       7.809   9.120  11.303  1.00 10.00           C
ATOM     47  O   TYR A  23       8.750   9.413  10.567  1.00 10.00           O
ATOM     48  CB  TYR A  23       7.103   6.730  11.024  1.00 10.00           C
ATOM     49  CG  TYR A  23       7.345   6.352  12.468  1.00 10.00           C
ATOM     50  CD1 TYR A  23       6.289   6.007  13.302  1.00 10.00           C
ATOM     51  CD2 TYR A  23       8.629   6.335  12.997  1.00 10.00           C
ATOM     52  CE1 TYR A  23       6.504   5.660  14.623  1.00 10.00           C
ATOM     53  CE2 TYR A  23       8.854   5.989  14.317  1.00 10.00           C
ATOM     54  CZ  TYR A  23       7.789   5.652  15.124  1.00 10.00           C
ATOM     55  OH  TYR A  23       8.008   5.307  16.438  1.00 10.00           O
ATOM     56  N   SER A  24       7.692   9.578  12.545  1.00 10.00           N
ATOM     57  CA  SER A  24       8.685  10.471  13.130  1.00 10.00           C
ATOM     58  C   SER A  24       9.375   9.821  14.325  1.00 10.00           C
ATOM     59  O   SER A  24       8.725   9.448  15.302  1.00 10.00           O
ATOM     60  CB  SER A  24       8.037  11.791  13.551  1.00 10.00           C
ATOM     61  OG  SER A  24       7.455  12.450  12.440  1.00 10.00           O
ATOM     62  N   LEU A  25      10.694   9.688  14.240  1.00 10.00           N
ATOM     63  CA  LEU A  25      11.474   9.083  15.315  1.00 10.00           C
ATOM     64  C   LEU A  25      11.794  10.103  16.403  1.00 10.00           C
ATOM     65  O   LEU A  25      10.928  10.472  17.196  1.00 10.00           O
ATOM     66  CB  LEU A  25      12.769   8.468  14.772  1.00 10.00           C
ATOM     67  CG  LEU A  25      12.681   7.150  13.995  1.00 10.00           C
ATOM     68  CD1 LEU A  25      12.147   7.356  12.583  1.00 10.00           C
ATOM     69  CD2 LEU A  25      14.037   6.460  13.961  1.00 10.00           C
TER
END
"""

def get_pdb_inputs(pdb_str):
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str.split('\n'))
  model = mmtbx.model.manager(model_input = pdb_inp)
  model.process(make_restraints=True)
  return group_args(
    ph  = model.get_hierarchy(),
    grm = model.get_restraints_manager(),
    xrs = model.get_xray_structure())

def get_tmo(inp, d_min):
  sel = inp.ph.atom_selection_cache().selection(
    string = "name CA or name CB or name N or name O or name C")
  f_calc = inp.xrs.select(
    selection=sel).structure_factors(d_min=d_min).f_calc()
  fft_map = f_calc.fft_map(resolution_factor = 0.25)
  fft_map.apply_sigma_scaling()
  target_map = fft_map.real_map_unpadded()
  return group_args(
    data             = target_map,
    miller_array     = f_calc,
    crystal_gridding = fft_map)

def shake_sites(xrs, random, shift, grm=None):
  from mmtbx.dynamics import cartesian_dynamics
  if(random):
    xrs.shake_sites_in_place(mean_distance = shift)
  else:
    grad_calc = cartesian_dynamics.gradients_calculator_geometry_restraints(
      restraints_manager = grm)
    cartesian_dynamics.run(
      xray_structure       = xrs,
      gradients_calculator = grad_calc,
      temperature          = 1000,
      n_steps              = 100000,
      time_step            = 0.0005,
      stop_cm_motion       = True,
      stop_at_diff         = shift)
  return xrs

def compute_map(target_map, xray_structure):
  mc = target_map.miller_array.structure_factors_from_scatterers(
    xray_structure = xray_structure).f_calc()
  fft_map = miller.fft_map(
    crystal_gridding     = target_map.crystal_gridding,
    fourier_coefficients = mc)
  fft_map.apply_sigma_scaling()
  return fft_map.real_map_unpadded()

def show(pdb_hierarchy, tm, xrs, grm, prefix):
  map = compute_map(target_map=tm, xray_structure=xrs)
  cc = flex.linear_correlation(
    x=map.as_1d(),
    y=tm.data.as_1d()).coefficient()
  es = grm.energies_sites(sites_cart = xrs.sites_cart())
  rmsd_a = es.angle_deviations()[2]
  rmsd_b = es.bond_deviations()[2]
  print("%s: overall CC: %6.4f rmsd_bonds=%6.3f rmsd_angles=%6.3f"%(
    prefix, cc, rmsd_b, rmsd_a))
  pdb_hierarchy.adopt_xray_structure(xrs)
  rotamer_manager = RotamerEval()
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue in chain.residues():
        sites_cart = residue.atoms().extract_xyz()
        sel = maptbx.grid_indices_around_sites(
          unit_cell  = xrs.unit_cell(),
          fft_n_real = map.focus(),
          fft_m_real = map.all(),
          sites_cart = sites_cart,
          site_radii = flex.double(sites_cart.size(), 2))
        ccr = flex.linear_correlation(
          x=map.select(sel).as_1d(),
          y=tm.data.select(sel).as_1d()).coefficient()
        fmt = "%s: %4s %10s CC: %6.4f"
        print(fmt%(prefix, residue.resname, rotamer_manager.evaluate_residue(residue),ccr))

def exercise(d_min=5, random_seed=1111111):
  inp = get_pdb_inputs(pdb_str=pdb_str)
  xrs_good = inp.xrs.deep_copy_scatterers()
  target_map = get_tmo(inp=inp, d_min = d_min)
  inp.ph.write_pdb_file(file_name="start.pdb")
  show(prefix="GOOD",
      pdb_hierarchy = inp.ph,
      tm=target_map,
      xrs=xrs_good,
      grm=inp.grm.geometry)
  #
  sites_cart_reference = []
  selections_reference = []
  pdb_hierarchy_reference = inp.ph.deep_copy()
  pdb_hierarchy_reference.reset_i_seq_if_necessary()
  for model in inp.ph.models():
    for chain in model.chains():
      for residue in chain.residues():
        sites_cart_reference.append(residue.atoms().extract_xyz())
        selections_reference.append(residue.atoms().extract_i_seq())
  #
  sites_cart_reference_for_chi_only = []
  selections_reference_for_chi_only = []
  for model in inp.ph.models():
    for chain in model.chains():
      for residue in chain.residues():
        s1 = flex.vec3_double()
        s2 = flex.size_t()
        for atom in residue.atoms():
          if(not atom.name.strip().upper() in ["O"]):
            s1.append(atom.xyz)
            s2.append(atom.i_seq)
        sites_cart_reference_for_chi_only.append(s1)
        selections_reference_for_chi_only.append(s2)
  #
  xrs_poor = shake_sites(xrs=xrs_good.deep_copy_scatterers(), random=False,
    shift=2.0, grm=inp.grm)
  inp.ph.adopt_xray_structure(xrs_poor)
  inp.ph.write_pdb_file(file_name="poor.pdb")
  #
  for use_reference_torsion in ["no", "yes_add_once", "yes_add_per_residue",
                                "yes_manual"]:
    es = inp.grm.energies_sites(sites_cart = xrs_good.sites_cart()) # it's essential to update grm
    inp.ph.adopt_xray_structure(xrs_poor)
    random.seed(random_seed)
    flex.set_random_seed(random_seed)
    print("*"*79)
    print("use_reference_torsion:", use_reference_torsion)
    print("*"*79)
    show(prefix="START",pdb_hierarchy = inp.ph, tm=target_map, xrs=xrs_poor, grm=inp.grm.geometry)
    #
    if(use_reference_torsion == "yes_add_per_residue"):
      inp.grm.geometry.remove_chi_torsion_restraints_in_place()
      for sites_cart, selection in zip(sites_cart_reference, selections_reference):
        inp.grm.geometry.add_chi_torsion_restraints_in_place(
          pdb_hierarchy   = pdb_hierarchy_reference,
          sites_cart      = sites_cart,
          selection       = selection,
          chi_angles_only = True,
          sigma           = 1)
    if(use_reference_torsion == "yes_add_once"):
      inp.grm.geometry.remove_chi_torsion_restraints_in_place()
      inp.grm.geometry.add_chi_torsion_restraints_in_place(
        pdb_hierarchy   = pdb_hierarchy_reference,
        sites_cart      = xrs_good.sites_cart(),
        chi_angles_only = True,
        sigma           = 1)
    if(use_reference_torsion == "yes_manual"):
      inp.grm.geometry.remove_chi_torsion_restraints_in_place()
      for sites_cart, selection in zip(sites_cart_reference_for_chi_only,
                                       selections_reference_for_chi_only):
        inp.grm.geometry.add_chi_torsion_restraints_in_place(
          pdb_hierarchy   = pdb_hierarchy_reference,
          sites_cart      = sites_cart,
          selection       = selection,
          chi_angles_only = True,
          sigma           = 1)
    #
    tmp = xrs_poor.deep_copy_scatterers()
    rsr_simple_refiner = individual_sites.simple(
      target_map                  = target_map.data,
      selection                   = flex.bool(tmp.scatterers().size(), True),
      real_space_gradients_delta  = d_min/4,
      max_iterations              = 500,
      geometry_restraints_manager = inp.grm.geometry)
    refined = individual_sites.refinery(
      refiner          = rsr_simple_refiner,
      optimize_weight  = True,
      xray_structure   = tmp,
      start_trial_weight_value = 50,
      rms_bonds_limit  = 0.02,
      rms_angles_limit = 2.0)
    assert refined.sites_cart_result is not None
    tmp = tmp.replace_sites_cart(refined.sites_cart_result)
    inp.ph.adopt_xray_structure(tmp)
    show(prefix="FINAL",pdb_hierarchy = inp.ph, tm=target_map, xrs=tmp, grm=inp.grm.geometry)
    inp.ph.write_pdb_file(file_name="final_%s.pdb"%str(use_reference_torsion))

if (__name__ == "__main__"):
  exercise()
  print("OK")
