from __future__ import division

import iotbx
from mmtbx.ligands.hierarchy_utils import _new_atom

def get_as_hierarchy(filename):
  model = iotbx.cif.reader(filename).model()
  for key,block in model.items():
    if key=='comp_list':
      code=block.get_loop_or_row('_chem_comp')
      code=code.get('_chem_comp.id')[0]
      continue
    loop = block.get_loop_or_row('_chem_comp_atom')
    ag = iotbx.pdb.hierarchy.atom_group()
    ag.resname=code

    for j, tmp in enumerate(loop.iterrows()):
      xyz = (float(tmp.get('_chem_comp_atom.x')),
             float(tmp.get('_chem_comp_atom.y')),
             float(tmp.get('_chem_comp_atom.z')),
             )
      atom = _new_atom(tmp['_chem_comp_atom.atom_id'],
                       tmp['_chem_comp_atom.type_symbol'],
                       xyz,
                       1.,
                       20.,
                       True,
                       )
      # atom.set_serial(j+1)
      ag.append_atom(atom)
  rg = iotbx.pdb.hierarchy.residue_group()
  rg.resseq='   1'
  rg.append_atom_group(ag)
  chain = iotbx.pdb.hierarchy.chain()
  chain.id='A'
  chain.append_residue_group(rg)
  model = iotbx.pdb.hierarchy.model()
  model.append_chain(chain)
  ph = iotbx.pdb.hierarchy.root()
  ph.append_model(model)
  ph.reset_atom_i_seqs()
  return ph

def as_cif_object(filename):
  model = iotbx.cif.reader(filename).model()
  return model
