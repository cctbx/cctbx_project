from __future__ import absolute_import, division, print_function
import os
from io import StringIO
import time

from libtbx.utils import Sorry
from scitbx.array_family import flex

# from cctbx.geometry_restraints.manager import manager as standard_manager

from mmtbx.geometry_restraints import base_qm_manager, mopac_manager, orca_manager

harkcal = 627.50946900
bohrang = 0.52918
#
# QM runner
#
def qm_runner(qmm,
              cleanup=True,
              file_read=False,
              log=StringIO(),
              ):
  def get_func(manager, attr):
    return getattr(manager, 'get_%s' % attr, None)
  redirect_output=True
  if qmm.program=='test':
    func = get_func(base_qm_manager.base_manager, qmm.program_goal)
  elif qmm.program=='orca':
    func = get_func(orca_manager.orca_manager, qmm.program_goal)
    coordinate_filename_ext='.xyz'
    log_filename_ext='.log'
    # raise Sorry('Orca temporarily unsupported. Consider using MOPAC.')
  elif qmm.program=='mopac':
    func = get_func(mopac_manager.mopac_manager, qmm.program_goal)
    coordinate_filename_ext='.arc'
    log_filename_ext='.out'
    redirect_output=False
  else:
    raise Sorry('QM program not found or set "%s"' % qmm.program)
  if func is None:
    raise Sorry('QM manager does not have get_%s' % qmm.program_goal)
  ligand_xyz, buffer_xyz = func(qmm,
                                cleanup=cleanup,
                                file_read=file_read,
                                coordinate_filename_ext=coordinate_filename_ext,
                                log_filename_ext=log_filename_ext,
                                redirect_output=redirect_output,
                                log=log)
  return ligand_xyz, buffer_xyz
#
# ORCA
#
'''
                                .--------------------.
          ----------------------|Geometry convergence|-------------------------
          Item                value                   Tolerance       Converged
          ---------------------------------------------------------------------
          Energy change      -0.2772205978            0.0000050000      NO
          RMS gradient        0.0273786248            0.0001000000      NO
          MAX gradient        0.1448471259            0.0003000000      NO
          RMS step            0.0097797205            0.0020000000      NO
          MAX step            0.0482825340            0.0040000000      NO
          ........................................................
          Max(Bonds)      0.0256      Max(Angles)    0.97
          Max(Dihed)        0.59      Max(Improp)    0.00
          ---------------------------------------------------------------------

          ----------------------------------------------------------------------------
                                  WARNING !!!
       The optimization did not converge but reached the maximum number of
       optimization cycles.
       Please check your results very carefully.
    ----------------------------------------------------------------------------
          '''
# def process_orca_convergence(lines):
#   s = ''
#   for line in lines:
#     tmp = line.split()
#     if tmp[-1] in ['YES', 'NO']:
#       # rc[tmp[0]]=tmp[-1]
#       s+= '%s ' % tmp[-1]
#   return s


# class manager(standard_manager):
#   def __init__(self,
#                params,
#                log=StringIO()):
#     # self.gradients_factory = gradients_factory
#     adopt_init_args(self, locals(), exclude=["log"])
#     self.validate()
#     assert 0

#   def validate(self):
#     qi = self.params.qi
#     assert qi.use_quantum_interface
#     assert qi.selection
#     if qi.orca.use_orca:
#       print('Orca')
#     assert 0

#   def get_engrad(self, sites_cart):
#     self.execution_manager.set_sites_cart(sites_cart)
#     return self.execution_manager.get_engrad()

#   def get_opt(self, sites_cart):
#     assert 0
#     self.execution_manager.set_sites_cart(sites_cart)
#     return self.execution_manager.get_opt()

#   def set_qm_info(self,
#                   method,
#                   basis_set,
#                   solvent_model,
#                   charge,
#                   multiplicity,
#                   ):
#     adopt_init_args(self, locals())
#     if self.basis_set is None:
#       self.basis_set = ''
#     if self.solvent_model is None:
#       self.solvent_model = ''
#     self.execution_manager = orca_manager( self.qm_atoms,
#                                            self.method,
#                                            self.basis_set,
#                                            self.solvent_model,
#                                            self.charge,
#                                            self.multiplicity
#                                            )

#   def set_qm_atoms(self, qm_atoms):
#     self.qm_atoms = qm_atoms
#     self.qm_iseqs = []
#     for atom in self.qm_atoms:
#       self.qm_iseqs.append(atom.i_seq)

#   def energies_sites(self,
#                      sites_cart,
#                      flags=None,
#                      custom_nonbonded_function=None,
#                      compute_gradients=False,
#                      gradients=None,
#                      disable_asu_cache=False,
#                      normalization=False,
#                      external_energy_function=None,
#                      extension_objects=[],
#                      site_labels=None,
#                      log=None):
#     result = standard_manager.energies_sites(
#       self,
#       sites_cart,
#       flags=flags,
#       custom_nonbonded_function=custom_nonbonded_function,
#       compute_gradients=compute_gradients,
    #   gradients=gradients,
    #   disable_asu_cache=disable_asu_cache,
    #   normalization=normalization,
    #   external_energy_function=external_energy_function,
    #   extension_objects=extension_objects,
    #   site_labels=site_labels,
    #   )
    # if compute_gradients:
    #   qm_sites_cart = []
    #   for i_seq in self.qm_iseqs:
    #     qm_sites_cart.append(sites_cart[i_seq])
    #   # coordinates = self.get_opt(qm_sites_cart)
    #   # print(list(coordinates))
    #   # assert 0
    #   energy, gradients = self.get_engrad(qm_sites_cart)
    #   for i_seq, gradient in zip(self.qm_iseqs, gradients):
    #     result.gradients[i_seq]=gradient
    # return result

def main():
  from iotbx import pdb
  pdb_lines = '''
HETATM   97  S   SO4 A  13      31.477  38.950  15.821  0.50 25.00           S
HETATM   98  O1  SO4 A  13      31.243  38.502  17.238  0.50 25.00           O
HETATM   99  O2  SO4 A  13      30.616  40.133  15.527  0.50 25.00           O
HETATM  100  O3  SO4 A  13      31.158  37.816  14.905  0.50 25.00           O
HETATM  101  O4  SO4 A  13      32.916  39.343  15.640  0.50 25.00           O
'''
  pdb_inp = pdb.input(lines=pdb_lines, source_info='lines')
  qi_grm = orca_manager(pdb_inp.atoms(),
                        'PM3',
                        '',
                        '',
                        -2,
                        1,
                        preamble='test',
                        )
  print(qi_grm)
  energy, gradients = qi_grm.get_engrad()
  print(energy, list(gradients))
  coordinates = qi_grm.get_opt()
  print(list(coordinates))

if __name__ == '__main__':
  main()
