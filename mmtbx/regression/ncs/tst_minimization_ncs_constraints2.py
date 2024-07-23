from __future__ import absolute_import, division, print_function
import mmtbx.refinement.minimization_ncs_constraints
from scitbx.array_family import flex
from libtbx import group_args
import iotbx.ncs as ncs
import iotbx.pdb
import time
from iotbx.ncs import ncs_group_master_phil
import iotbx.phil
from six.moves import range
import mmtbx.model

pdb_answer_0 = """\
CRYST1   18.415   14.419   12.493  90.00  90.00  90.00 P 1
ATOM      5  O   HOH S   3       5.648   5.054   4.347  1.00 10.00           O
ATOM      6  O   HOH S   4       2.000   5.459   5.854  1.00 10.00           O
ATOM      7  O   HOH S   5       3.498   4.230   8.986  1.00 10.00           O
TER
ATOM      1  N   THR A   1      11.782  12.419   4.645  1.00 10.00           N
ATOM      2  CA  THR A   1      11.671  11.061   4.125  1.00 10.00           C
ATOM      3  C   THR A   1      11.746  10.033   5.249  1.00 10.00           C
ATOM      4  O   THR A   1      12.561  10.157   6.163  1.00 10.00           O
ATOM      5  CB  THR A   1      12.772  10.760   3.092  1.00 10.00           C
ATOM      6  OG1 THR A   1      12.672  11.682   2.000  1.00 10.00           O
ATOM      7  CG2 THR A   1      12.635   9.340   2.565  1.00 10.00           C
TER
ATOM      1  N   THR D   1      13.010   5.595  10.010  1.00 10.00           N
ATOM      2  CA  THR D   1      14.035   5.945   9.034  1.00 10.00           C
ATOM      3  C   THR D   1      15.310   6.423   9.720  1.00 10.00           C
ATOM      4  O   THR D   1      16.415   6.054   9.323  1.00 10.00           O
ATOM      5  CB  THR D   1      13.543   7.038   8.066  1.00 10.00           C
ATOM      6  OG1 THR D   1      13.237   8.229   8.802  1.00 10.00           O
ATOM      7  CG2 THR D   1      12.300   6.572   7.323  1.00 10.00           C
TER
ATOM      1  N   THR B   1       6.768   9.093   9.237  1.00 10.00           N
ATOM      2  CA  THR B   1       7.284   8.654   7.945  1.00 10.00           C
ATOM      3  C   THR B   1       8.638   7.968   8.097  1.00 10.00           C
ATOM      4  O   THR B   1       9.495   8.426   8.852  1.00 10.00           O
ATOM      5  CB  THR B   1       7.423   9.832   6.963  1.00 10.00           C
ATOM      6  OG1 THR B   1       6.144  10.446   6.765  1.00 10.00           O
ATOM      7  CG2 THR B   1       7.962   9.350   5.625  1.00 10.00           C
TER
ATOM      1  N   THR C   1       9.093   2.000  10.493  1.00 10.00           N
ATOM      2  CA  THR C   1       8.879   2.702   9.233  1.00 10.00           C
ATOM      3  C   THR C   1      10.081   3.570   8.875  1.00 10.00           C
ATOM      4  O   THR C   1      10.652   4.241   9.734  1.00 10.00           O
ATOM      5  CB  THR C   1       7.618   3.584   9.284  1.00 10.00           C
ATOM      6  OG1 THR C   1       6.472   2.770   9.559  1.00 10.00           O
ATOM      7  CG2 THR C   1       7.417   4.305   7.960  1.00 10.00           C
TER
ATOM      3  O   HOH E   1       4.021   6.412   7.953  1.00 10.00           O
ATOM      4  O   HOH E   2       4.186   2.319   6.140  1.00 10.00           O
TER
END
"""
pdb_poor_0 = """\
CRYST1   18.415   14.419   12.493  90.00  90.00  90.00 P 1
ATOM      5  O   HOH S   3       5.648   5.054   4.347  1.00 10.00           O
ATOM      6  O   HOH S   4       2.000   5.459   5.854  1.00 10.00           O
ATOM      7  O   HOH S   5       3.498   4.230   8.986  1.00 10.00           O
TER
ATOM      1  N   THR A   1      11.782  12.419   4.645  1.00 10.00           N
ATOM      2  CA  THR A   1      11.671  11.061   4.125  1.00 10.00           C
ATOM      3  C   THR A   1      11.746  10.033   5.249  1.00 10.00           C
ATOM      4  O   THR A   1      12.561  10.157   6.163  1.00 10.00           O
ATOM      5  CB  THR A   1      12.772  10.760   3.092  1.00 10.00           C
ATOM      6  OG1 THR A   1      12.672  11.682   2.000  1.00 10.00           O
ATOM      7  CG2 THR A   1      12.635   9.340   2.565  1.00 10.00           C
TER
ATOM      1  N   THR D   1      13.010   5.595  10.010  1.00 10.00           N
ATOM      2  CA  THR D   1      14.035   5.945   9.034  1.00 10.00           C
ATOM      3  C   THR D   1      15.310   6.423   9.720  1.00 10.00           C
ATOM      4  O   THR D   1      16.415   6.054   9.323  1.00 10.00           O
ATOM      5  CB  THR D   1      13.543   7.038   8.066  1.00 10.00           C
ATOM      6  OG1 THR D   1      13.237   8.229   8.802  1.00 10.00           O
ATOM      7  CG2 THR D   1      12.300   6.572   7.323  1.00 10.00           C
TER
ATOM      1  N   THR B   1       6.279   9.382   8.893  1.00 10.00           N
ATOM      2  CA  THR B   1       7.095   8.715   7.884  1.00 10.00           C
ATOM      3  C   THR B   1       8.190   7.873   8.531  1.00 10.00           C
ATOM      4  O   THR B   1       8.828   8.299   9.493  1.00 10.00           O
ATOM      5  CB  THR B   1       7.741   9.729   6.922  1.00 10.00           C
ATOM      6  OG1 THR B   1       6.718  10.485   6.262  1.00 10.00           O
ATOM      7  CG2 THR B   1       8.586   9.011   5.880  1.00 10.00           C
TER
ATOM      1  N   THR C   1       8.953   1.712   9.905  1.00 10.00           N
ATOM      2  CA  THR C   1       8.573   2.732   8.935  1.00 10.00           C
ATOM      3  C   THR C   1       9.744   3.658   8.623  1.00 10.00           C
ATOM      4  O   THR C   1      10.481   4.067   9.520  1.00 10.00           O
ATOM      5  CB  THR C   1       7.385   3.573   9.434  1.00 10.00           C
ATOM      6  OG1 THR C   1       6.258   2.721   9.672  1.00 10.00           O
ATOM      7  CG2 THR C   1       7.008   4.629   8.406  1.00 10.00           C
TER
ATOM      3  O   HOH E   1       4.021   6.412   7.953  1.00 10.00           O
ATOM      4  O   HOH E   2       4.186   2.319   6.140  1.00 10.00           O
TER      39      HOH E   2
END
"""

pdb_poor_1 = """\
CRYST1   18.415   14.419   12.493  90.00  90.00  90.00 P 1
ATOM      5  O   HOH S   3       5.648   5.054   4.347  1.00 10.00           O
ATOM      6  O   HOH S   4       2.000   5.459   5.854  1.00 10.00           O
ATOM      7  O   HOH S   5       3.498   4.230   8.986  1.00 10.00           O
TER       4      HOH S   5
ATOM      1  N   THR A   1      11.797  12.521   4.849  1.00 10.00           N
ATOM      2  CA  THR A   1      11.607  11.311   4.057  1.00 10.00           C
ATOM      3  C   THR A   1      11.119  10.155   4.925  1.00 10.00           C
ATOM      4  O   THR A   1      11.603   9.957   6.039  1.00 10.00           O
ATOM      5  CB  THR A   1      12.906  10.893   3.344  1.00 10.00           C
ATOM      6  OG1 THR A   1      13.340  11.949   2.478  1.00 10.00           O
ATOM      7  CG2 THR A   1      12.683   9.631   2.525  1.00 10.00           C
TER      12      THR A   1
ATOM      1  N   THR D   1      13.010   5.595  10.010  1.00 10.00           N
ATOM      2  CA  THR D   1      14.035   5.945   9.034  1.00 10.00           C
ATOM      3  C   THR D   1      15.310   6.423   9.720  1.00 10.00           C
ATOM      4  O   THR D   1      16.415   6.054   9.323  1.00 10.00           O
ATOM      5  CB  THR D   1      13.543   7.038   8.066  1.00 10.00           C
ATOM      6  OG1 THR D   1      13.237   8.229   8.802  1.00 10.00           O
ATOM      7  CG2 THR D   1      12.300   6.572   7.323  1.00 10.00           C
TER      20      THR D   1
ATOM      1  N   THR B   1       6.768   9.093   9.237  1.00 10.00           N
ATOM      2  CA  THR B   1       7.284   8.654   7.945  1.00 10.00           C
ATOM      3  C   THR B   1       8.638   7.968   8.097  1.00 10.00           C
ATOM      4  O   THR B   1       9.495   8.426   8.852  1.00 10.00           O
ATOM      5  CB  THR B   1       7.423   9.832   6.963  1.00 10.00           C
ATOM      6  OG1 THR B   1       6.144  10.446   6.765  1.00 10.00           O
ATOM      7  CG2 THR B   1       7.962   9.350   5.625  1.00 10.00           C
TER      28      THR B   1
ATOM      1  N   THR C   1       9.175   1.430  10.253  1.00 10.00           N
ATOM      2  CA  THR C   1       9.097   2.543   9.315  1.00 10.00           C
ATOM      3  C   THR C   1      10.328   3.438   9.417  1.00 10.00           C
ATOM      4  O   THR C   1      10.801   3.736  10.513  1.00 10.00           O
ATOM      5  CB  THR C   1       7.835   3.394   9.551  1.00 10.00           C
ATOM      6  OG1 THR C   1       6.668   2.578   9.396  1.00 10.00           O
ATOM      7  CG2 THR C   1       7.777   4.547   8.561  1.00 10.00           C
TER      36      THR C   1
ATOM      3  O   HOH E   1       4.021   6.412   7.953  1.00 10.00           O
ATOM      4  O   HOH E   2       4.186   2.319   6.140  1.00 10.00           O
TER      39      HOH E   2
END
"""

pdb_poor_2 = """\
CRYST1   18.415   14.419   12.493  90.00  90.00  90.00 P 1
ATOM      1  O   HOH S   3       5.648   5.054   4.347  1.00 10.00           O
ATOM      2  O   HOH S   4       2.000   5.459   5.854  1.00 10.00           O
ATOM      3  O   HOH S   5       3.498   4.230   8.986  1.00 10.00           O
TER
ATOM      4  N   THR A   1      11.360  12.428   4.988  1.00 10.00           N
ATOM      5  CA  THR A   1      11.640  10.782   4.445  1.00 10.00           C
ATOM      6  C   THR A   1      11.410  10.474   5.073  1.00 10.00           C
ATOM      7  O   THR A   1      12.394  10.079   6.303  1.00 10.00           O
ATOM      8  CB  THR A   1      13.109  11.143   2.643  1.00 10.00           C
ATOM      9  OG1 THR A   1      12.760  11.290   2.263  1.00 10.00           O
ATOM     10  CG2 THR A   1      12.701   9.443   2.402  1.00 10.00           C
TER
ATOM     11  N   THR D   1      13.010   5.595  10.010  1.00 10.00           N
ATOM     12  CA  THR D   1      14.035   5.945   9.034  1.00 10.00           C
ATOM     13  C   THR D   1      15.310   6.423   9.720  1.00 10.00           C
ATOM     14  O   THR D   1      16.415   6.054   9.323  1.00 10.00           O
ATOM     15  CB  THR D   1      13.543   7.038   8.066  1.00 10.00           C
ATOM     16  OG1 THR D   1      13.237   8.229   8.802  1.00 10.00           O
ATOM     17  CG2 THR D   1      12.300   6.572   7.323  1.00 10.00           C
TER
ATOM     18  N   THR B   1       6.970   8.670   9.020  1.00 10.00           N
ATOM     19  CA  THR B   1       7.512   9.099   7.575  1.00 10.00           C
ATOM     20  C   THR B   1       9.003   7.644   7.782  1.00 10.00           C
ATOM     21  O   THR B   1       9.316   8.109   8.577  1.00 10.00           O
ATOM     22  CB  THR B   1       7.065  10.156   6.729  1.00 10.00           C
ATOM     23  OG1 THR B   1       6.149  10.590   6.419  1.00 10.00           O
ATOM     24  CG2 THR B   1       7.542   9.150   6.053  1.00 10.00           C
TER
ATOM     25  N   THR C   1       9.497   2.326  10.790  1.00 10.00           N
ATOM     26  CA  THR C   1       8.804   2.690   9.478  1.00 10.00           C
ATOM     27  C   THR C   1      10.315   3.465   9.108  1.00 10.00           C
ATOM     28  O   THR C   1      10.242   3.818   9.832  1.00 10.00           O
ATOM     29  CB  THR C   1       7.933   3.341   9.276  1.00 10.00           C
ATOM     30  OG1 THR C   1       6.866   2.882   9.784  1.00 10.00           O
ATOM     31  CG2 THR C   1       7.545   4.707   8.366  1.00 10.00           C
TER
ATOM     32  O   HOH E   1       4.021   6.412   7.953  1.00 10.00           O
ATOM     33  O   HOH E   2       4.186   2.319   6.140  1.00 10.00           O
TER
END
"""

ncs_params_str_0 = """
ncs_group {
  reference = chain A
  selection = chain B
  selection = chain C
}
"""

ncs_params_str_1 = """
ncs_group {
  reference = chain B
  selection = chain C
  selection = chain A
}
"""

def set_scattering_dictionary(xray_structure, d_min,
  scattering_table="n_gaussian"):
  xray_structure.scattering_type_registry(
    table = scattering_table,
    d_min = d_min,
    types_without_a_scattering_contribution=["?"])

def get_inputs(prefix, pdb_answer, pdb_poor, ncs_params_str, real_space, d_min):
  pdb_file_name_answer = "answer_%s.pdb"%prefix
  of=open(pdb_file_name_answer, "w")
  print(pdb_answer, file=of)
  of.close()
  #
  pdb_file_name_poor = "poor_%s.pdb"%prefix
  of=open(pdb_file_name_poor, "w")
  print(pdb_poor, file=of)
  of.close()
  #
  pdb_inp_answer = iotbx.pdb.input(file_name=pdb_file_name_answer)
  ph_answer = pdb_inp_answer.construct_hierarchy()
  ph_answer.atoms().reset_i_seq()
  xrs_answer = pdb_inp_answer.xray_structure_simple()
  sites_cart_answer = xrs_answer.sites_cart()
  #
  pdb_inp_poor = iotbx.pdb.input(file_name=pdb_file_name_poor)
  ph_poor = pdb_inp_poor.construct_hierarchy()
  ph_poor.atoms().reset_i_seq()
  xrs_poor = pdb_inp_poor.xray_structure_simple()
  #
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_poor)
  model = mmtbx.model.manager(model_input=pdb_inp)
  p = model.get_default_pdb_interpretation_params()
  p.pdb_interpretation.const_shrink_donor_acceptor=0.6
  model.process(pdb_interpretation_params=p, make_restraints=True)
  restraints_manager = model.get_restraints_manager()
  restraints_manager.geometry.remove_c_beta_torsion_restraints_in_place()
  #
  phil_groups = ncs_group_master_phil.fetch(
      iotbx.phil.parse(ncs_params_str)).extract()

  pdb_inp = iotbx.pdb.input(lines=pdb_answer,source_info=None)
  ncs_inp = ncs.input(
      hierarchy=pdb_inp.construct_hierarchy(),
      ncs_phil_groups=phil_groups.ncs_group)
  ncs_groups = ncs_inp.get_ncs_restraints_group_list()
  # print "ncs_groups:", len(ncs_groups)
  # print "master isel:", list(ncs_groups[0].master_iselection)
  # for c in ncs_groups[0].copies:
  #   print "copy isel:", list(c.iselection)
  #
  set_scattering_dictionary(xray_structure = xrs_answer, d_min = d_min)
  set_scattering_dictionary(xray_structure = xrs_poor,   d_min = d_min)
  #
  #
  map_data, fmodel = None, None
  if(real_space):
    fc = xrs_answer.structure_factors(d_min=d_min, algorithm="direct").f_calc()
    fft_map = fc.fft_map(resolution_factor = 0.25)
    fft_map.apply_sigma_scaling()
    map_data = fft_map.real_map_unpadded()
  else:
    f_obs = abs(xrs_answer.structure_factors(d_min=d_min,
      algorithm="direct").f_calc())
    params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
    params.algorithm = "direct"
    fmodel = mmtbx.f_model.manager(
      f_obs                        = f_obs,
      xray_structure               = xrs_poor,
      sf_and_grads_accuracy_params = params,
      target_name                  = "ls_wunit_k1")
    if(1): print("d_min:", d_min)
    if(1): print("r_work(start):", fmodel.r_work())
  return group_args(
    fmodel             = fmodel,
    map_data           = map_data,
    xrs_poor           = xrs_poor,
    d_min              = d_min,
    ncs_groups         = ncs_groups,
    restraints_manager = restraints_manager,
    ph                 = ph_poor)

def macro_cycle(
      restraints_manager,
      fmodel,
      map_data,
      d_min,
      xrs_poor,
      pdb_poor,
      ncs_groups,
      refine_selection,
      actions):
  assert [fmodel, map_data].count(None) == 1
  for cycle in range(100):
    #for action in actions:
    #  refine_sites, refine_transformations = action
    data_weight = 1
    if(fmodel is not None):
      tg_object = mmtbx.refinement.minimization_ncs_constraints.\
        target_function_and_grads_reciprocal_space(
          fmodel                    = fmodel,
          ncs_restraints_group_list = ncs_groups,
          refine_selection          = refine_selection,
          restraints_manager        = restraints_manager,
          data_weight               = data_weight,
          refine_sites              = True)
    else:
      tg_object = mmtbx.refinement.minimization_ncs_constraints.\
        target_function_and_grads_real_space(
          map_data                   = map_data,
          xray_structure             = xrs_poor,
          ncs_restraints_group_list  = ncs_groups,
          refine_selection           = refine_selection,
          real_space_gradients_delta = d_min/4,
          restraints_manager         = restraints_manager,
          data_weight                = data_weight,
          refine_sites               = True)
    minimized = mmtbx.refinement.minimization_ncs_constraints.lbfgs(
      target_and_grads_object      = tg_object,
      xray_structure               = xrs_poor,
      ncs_restraints_group_list    = ncs_groups,
      refine_selection             = refine_selection,
      finite_grad_differences_test = False,
      max_iterations               = 100,
      refine_sites                 = True)
    if(fmodel is not None):
      fmodel.update_xray_structure(
        xray_structure = xrs_poor, update_f_calc=True)
      if(0): print(cycle, fmodel.r_work())
  if(fmodel is not None):
    rf = fmodel.r_work()
    print("R(final):", rf)
    return rf
  else:
    return None

def call(prefix, refine_selection, pdb_poor, pdb_answer, ncs_params_str,
         actions, real_space, d_min):
  inp = get_inputs(
    prefix         = prefix,
    pdb_poor       = pdb_poor,
    pdb_answer     = pdb_answer,
    ncs_params_str = ncs_params_str,
    real_space     = real_space,
    d_min          = d_min)
  ph_poor = inp.ph.deep_copy()
  ph_poor.adopt_xray_structure(inp.xrs_poor)
  result = macro_cycle(
    restraints_manager = inp.restraints_manager,
    fmodel             = inp.fmodel,
    map_data           = inp.map_data,
    d_min              = inp.d_min,
    xrs_poor           = inp.xrs_poor,
    pdb_poor           = ph_poor,
    ncs_groups         = inp.ncs_groups,
    refine_selection   = refine_selection,
    actions            = actions)
  inp.ph.adopt_xray_structure(inp.xrs_poor)
  inp.ph.write_pdb_file(file_name="%s.pdb"%prefix)
  return result

def check_result(result_file_name):
  s1=iotbx.pdb.input(file_name=result_file_name).atoms().extract_xyz()
  s2=iotbx.pdb.input(file_name="answer_"+result_file_name).atoms().extract_xyz()
  r = flex.sqrt((s1 - s2).dot()).min_max_mean().as_tuple()
  print(r)
  return r

def run():
  """
  NCS constrained refinement of coordinates. Chains A, B and C are
  NCS related, chain D is refined without NCS constraints, chain E is fixed.
  Real- and reciprocal-space refinement.
  THIS TEST NOW DOES NOT EXERCISE REFINEMENT OF NCS OPERATORS !!!
  """
  refine_selection = flex.size_t(range(3,31))
  for real_space in [True, False]:
    print("real_space:", real_space)
    if(real_space):
      suffix = "Real"
      d_min  = 1.0
    else:
      suffix = "Reciprocal"
      d_min  = 1.5
    if(1):
      #
      # Refine NCS constrained sites and NCS operators.
      #
      prefix = "tst_%s_CoordinatesOnly"%suffix
      rf = call(
        prefix           = prefix,
        refine_selection = refine_selection,
        pdb_poor         = pdb_poor_2,
        pdb_answer       = pdb_answer_0,
        ncs_params_str   = ncs_params_str_0,
        actions          = [[False,True],[True,False]],
        real_space       = real_space,
        d_min            = d_min)
      r = check_result(result_file_name=prefix+".pdb")
      if(real_space):
        assert r[1]<0.1  , r[1]
        assert r[2]<0.05 , r[2]
        assert rf is None
      else:
        assert r[1]<0.05    , r[1]
        assert r[2]<0.0055  , r[2]
        assert rf  < 0.005 , rf

if (__name__ == "__main__"):
  t0=time.time()
  run()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
