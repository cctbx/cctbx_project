from __future__ import division

import os
import libtbx.load_env

import iotbx
from mmtbx.ligands.hierarchy_utils import _new_atom

phenix_repository_dir = os.path.dirname(libtbx.env.dist_path("phenix"))
geostd_directory = os.path.join(phenix_repository_dir,
                                "chem_data",
                                "geostd",
                                )

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

def remove_atoms(cif_object, names):
  for s, b in cif_object.items():
    if s=='comp_list': continue
    for loop in b.iterloops():
      remove=[]
      for i, row in enumerate(loop.iterrows()):
        for key, item in row.items():
          if key.find('atom_id')>-1:
            if item in names:
              remove.append(i)
              break
      if remove:
        remove.reverse()
        for r in remove:
          loop.delete_row(r)

def remove_atoms_for_reduce(cif_object):
  for s, b in cif_object.items():
    for key, item in b.items():
      if key=='_chem_comp.group':
        assert len(item)==1
        if item[0] in ['RNA', 'DNA']:
          remove_atoms(cif_object, ["HO3'", 'HOP2'])
          return True
    assert 0

def get_geostd_file(code, pH=None, file_format='cif'):
  assert code
  if pH: assert file_format=='cif', 'pH=%s only for cif' % pH
  basename = "%s.%s" % (code.upper(), file_format)
  if file_format=='cif':
    basename='data_%s' % basename
  if pH:
    basename = basename.replace('.cif', '_pH_%s.cif' % pH)
  filename = os.path.join(geostd_directory,
                          code[0].lower(),
                          basename,
                          )
  if os.path.exists(filename):
    return filename
  return None

def get_geostd_cif_file(code, pH=None):
  return get_geostd_file(code, pH=pH, file_format='cif')

def get_cif_list(filename, verbose=False):
  import iotbx.cif
  list_cif = os.path.join(geostd_directory,
                          "list",
                          filename,
                          )
  cif = iotbx.cif.reader(list_cif, strict=False).model()
  if verbose: print(cif.keys())
  return cif
