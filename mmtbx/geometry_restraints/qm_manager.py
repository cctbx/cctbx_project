from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry
from mmtbx.geometry_restraints import base_qm_manager, mopac_manager, orca_manager

harkcal = 627.50946900
bohrang = 0.52918
#
# QM runner
#
def qm_runner(qmm,
              cleanup=True,
              file_read=False,
              check_file_read_safe=True,
              log=None,
              ):
  def get_func(manager, attr):
    return getattr(manager, 'get_%s' % attr, None)
  redirect_output=True
  # fake
  coordinate_filename_ext='.nwm'
  log_filename_ext='.nwm'
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
                                check_file_read_safe=check_file_read_safe,
                                coordinate_filename_ext=coordinate_filename_ext,
                                log_filename_ext=log_filename_ext,
                                redirect_output=redirect_output,
                                log=log)
  return ligand_xyz, buffer_xyz

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
