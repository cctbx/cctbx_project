from __future__ import absolute_import, division, print_function
from scitbx import matrix
from rstbx.dps_core import dps_core
from rstbx.dps_core.sampling import hemisphere_shortcut
from rstbx.dps_core.basis_choice import SelectBasisMetaprocedure

class algorithm_parameters:
  max_cell_edge_for_dps_fft = 100.0 # in Angstroms.
                                    # 100 Angstrom should be suitable for small molecular work
                                    # use Bragg spot spacing to estimate the maximum
                                    # cell edge for macromolecular cases
  directional_sampling_granularity = 0.029 # in radians; 0.029 radian granularity
                                           # is both fast enough and granular enough
                                           # to sample the hemisphere for small molecular work.
                                           # Tradeoff between performance and coverage occurs
                                           # for macromolecular work; DPS paper uses 0.029
  max_cell_edge_basis_choice = 8.0  # in Angstroms. This input parameter is important.
                                    # choice of basis vector is extremely sensitive to
                                    # this parameter--for silicon, must choose a value less
                                    # than twice the cell edge
def do_index(reciprocal_space_vectors, verbose=True):
  D = dps_core()
  D.setMaxcell(algorithm_parameters.max_cell_edge_for_dps_fft)
  D.setXyzData(reciprocal_space_vectors)
  hemisphere_shortcut(ai = D,
    characteristic_sampling = algorithm_parameters.directional_sampling_granularity,
    max_cell = algorithm_parameters.max_cell_edge_basis_choice)
  M = SelectBasisMetaprocedure(D)

  from rstbx.dps_core.lepage import iotbx_converter
  L = iotbx_converter(D.getOrientation().unit_cell().minimum_cell(),5.0)
  supergroup = L[0]

  triclinic = D.getOrientation().unit_cell()

  cb_op = supergroup['cb_op_inp_best'].c().as_double_array()[0:9]
  orient = D.getOrientation()
  orient_best = orient.change_basis(matrix.sqr(cb_op).transpose())
  constrain_orient = orient_best.constrain(supergroup['system'])
  D.setOrientation(constrain_orient)

  if verbose:
    for subgroup in L:
      print(subgroup.short_digest())
    print("\ntriclinic cell=%s volume(A^3)=%.3f"%(triclinic,triclinic.volume()))
    print("\nafter symmetrizing to %s:"%supergroup.reference_lookup_symbol())
    M.show_rms()
  return D,L
