from __future__ import absolute_import, division, print_function

from libtbx.utils import format_cpu_times
from libtbx.test_utils import show_diff
from iotbx import pdb

def awl(a):
  return a.fetch_labels()

def exercise_awl_id_str():
  # The same as atom id_str test, but converting to atom_with_labels for all checks.
  a = pdb.hierarchy.atom()
  a.set_name(new_name="NaMe")
  a.set_serial(new_serial="B1234")
  assert awl(a).id_str() == 'pdb="NaMe           "'
  assert awl(a).id_str(pdbres=True) == 'pdbres="          "'
  ag = pdb.hierarchy.atom_group(altloc="A", resname="longGLY")
  ag.append_atom(a)
  assert awl(a).id_str() == 'pdb="NaMeAlongGLY       "', "'%s'" % awl(a).id_str()
  assert awl(a).id_str(True) == 'pdbres="longGLY       "'
  rg = pdb.hierarchy.residue_group(resseq="1234", icode="J")
  rg.append_atom_group(ag)
  assert awl(a).id_str() == 'pdb="NaMeAlongGLY  1234J"'
  assert awl(a).id_str(pdbres=True) == 'pdbres="longGLY  1234J"'
  ch = pdb.hierarchy.chain(id="dlinCh")
  ch.append_residue_group(rg)
  assert awl(a).id_str() == 'pdb="NaMeAlongGLYdlinCh1234J"'
  assert awl(a).id_str(pdbres=True) == 'pdbres="longGLYdlinCh1234J"'
  md = pdb.hierarchy.model()
  md.append_chain(ch)
  assert awl(a).id_str() == 'pdb="NaMeAlongGLYdlinCh1234J"'
  assert awl(a).id_str(pdbres=True) == 'pdbres="longGLYdlinCh1234J"'
  md.id = ""
  assert awl(a).id_str() == 'pdb="NaMeAlongGLYdlinCh1234J"'
  assert awl(a).id_str(pdbres=True) == 'pdbres="longGLYdlinCh1234J"'
  md.id = "1"
  assert awl(a).id_str() == 'model="   1" pdb="NaMeAlongGLYdlinCh1234J"'
  assert awl(a).id_str(pdbres=True) == 'model="   1" pdbres="longGLYdlinCh1234J"'
  md.id = "12345678"
  assert awl(a).id_str() == 'model="12345678" pdb="NaMeAlongGLYdlinCh1234J"'
  assert awl(a).id_str(pdbres=True) == 'model="12345678" pdbres="longGLYdlinCh1234J"'

  md.id = "12345678"
  a.segid = "1234"

  assert awl(a).id_str(suppress_segid=False) \
      == 'model="12345678" pdb="NaMeAlongGLYdlinCh1234J" segid="1234"', a.id_str(suppress_segid=False)
  assert awl(a).id_str(suppress_segid=True) \
      == 'model="12345678" pdb="NaMeAlongGLYdlinCh1234J"'
  assert awl(a).id_str(pdbres=True) \
      == 'model="12345678" pdbres="longGLYdlinCh1234J" segid="1234"'
  assert awl(a).id_str(pdbres=True, suppress_segid=True) \
      == 'model="12345678" pdbres="longGLYdlinCh1234J"'
  md.id = ""
  assert awl(a).id_str() == 'pdb="NaMeAlongGLYdlinCh1234J" segid="1234"'
  assert awl(a).id_str(pdbres=True) \
      == 'pdbres="longGLYdlinCh1234J" segid="1234"'
  assert awl(a).id_str(pdbres=True, suppress_segid=True) \
      == 'pdbres="longGLYdlinCh1234J"'
  rt = pdb.hierarchy.root()
  rt.append_model(md)
  assert awl(a).id_str() == 'pdb="NaMeAlongGLYdlinCh1234J" segid="1234"'
  assert awl(a).id_str(pdbres=True) == 'pdbres="longGLYdlinCh1234J" segid="1234"'
  md.id = "    "
  assert awl(a).id_str() == 'model="    " pdb="NaMeAlongGLYdlinCh1234J" segid="1234"'
  assert awl(a).id_str(pdbres=True) \
      == 'model="    " pdbres="longGLYdlinCh1234J" segid="1234"'


def exercise_atom_residue_id_str():
  a = pdb.hierarchy.atom()
  a.set_name(new_name="NaMe")
  a.set_serial(new_serial="B1234")
  assert a.id_str() == 'pdb="NaMe           "'
  assert a.id_str(pdbres=True) == 'pdbres="          "'
  ag = pdb.hierarchy.atom_group(altloc="A", resname="longGLY")
  ag.append_atom(a)
  assert a.id_str() == 'pdb="NaMeAlongGLY       "'
  assert a.id_str(True) == 'pdbres="longGLY       "'

  rg = pdb.hierarchy.residue_group(resseq="1234", icode="J")
  rg.append_atom_group(ag)
  assert a.id_str() == 'pdb="NaMeAlongGLY  1234J"'
  assert a.id_str(pdbres=True) == 'pdbres="longGLY  1234J"'
  ch = pdb.hierarchy.chain(id="dlinCh")
  ch.append_residue_group(rg)
  assert a.id_str() == 'pdb="NaMeAlongGLYdlinCh1234J"'
  assert a.id_str(pdbres=True) == 'pdbres="longGLYdlinCh1234J"'
  md = pdb.hierarchy.model()
  md.append_chain(ch)
  assert a.id_str() == 'pdb="NaMeAlongGLYdlinCh1234J"'
  assert a.id_str(pdbres=True) == 'pdbres="longGLYdlinCh1234J"'
  md.id = ""
  assert a.id_str() == 'pdb="NaMeAlongGLYdlinCh1234J"'
  assert a.id_str(pdbres=True) == 'pdbres="longGLYdlinCh1234J"'
  md.id = "1"
  assert a.id_str() == 'model="   1" pdb="NaMeAlongGLYdlinCh1234J"'
  assert a.id_str(pdbres=True) == 'model="   1" pdbres="longGLYdlinCh1234J"'
  md.id = "12345678"
  assert a.id_str() == 'model="12345678" pdb="NaMeAlongGLYdlinCh1234J"'
  assert a.id_str(pdbres=True) == 'model="12345678" pdbres="longGLYdlinCh1234J"'

  md.id = "12345678"
  a.segid = "1234"

  assert a.id_str(suppress_segid=False) \
      == 'model="12345678" pdb="NaMeAlongGLYdlinCh1234J" segid="1234"', a.id_str(suppress_segid=False)
  assert a.id_str(suppress_segid=True) \
      == 'model="12345678" pdb="NaMeAlongGLYdlinCh1234J"'
  assert a.id_str(pdbres=True) \
      == 'model="12345678" pdbres="longGLYdlinCh1234J" segid="1234"'
  assert a.id_str(pdbres=True, suppress_segid=True) \
      == 'model="12345678" pdbres="longGLYdlinCh1234J"'
  md.id = ""
  assert a.id_str() == 'pdb="NaMeAlongGLYdlinCh1234J" segid="1234"'
  assert a.id_str(pdbres=True) \
      == 'pdbres="longGLYdlinCh1234J" segid="1234"'
  assert a.id_str(pdbres=True, suppress_segid=True) \
      == 'pdbres="longGLYdlinCh1234J"'
  rt = pdb.hierarchy.root()
  rt.append_model(md)
  assert a.id_str() == 'pdb="NaMeAlongGLYdlinCh1234J" segid="1234"'
  assert a.id_str(pdbres=True) == 'pdbres="longGLYdlinCh1234J" segid="1234"'
  md.id = "    "
  assert a.id_str() == 'model="    " pdb="NaMeAlongGLYdlinCh1234J" segid="1234"'
  assert a.id_str(pdbres=True) \
      == 'model="    " pdbres="longGLYdlinCh1234J" segid="1234"'
  #
  cf = ch.only_conformer()
  rd = cf.only_residue()


  # residue
  assert rd.id_str() == 'model="    " pdbres="longGLYdlinCh1234J" segid="1234"', rd.id_str()
  assert rd.id_str(suppress_segid=1) == 'model="    " pdbres="longGLYdlinCh1234J"'
  md.id = "12345678"
  assert rd.id_str() == 'model="12345678" pdbres="longGLYdlinCh1234J" segid="1234"'
  assert rd.id_str(suppress_segid=1) == 'model="12345678" pdbres="longGLYdlinCh1234J"'
  del cf
  assert rd.id_str(suppress_segid=-1) == 'pdbres="longGLY  1234J" segid="1234"', rd.id_str(suppress_segid=-1)
  assert rd.id_str(suppress_segid=1) == 'pdbres="longGLY  1234J"'
  #
  a2 = pdb.hierarchy.atom().set_segid(new_segid="abcd")
  ag.append_atom(atom=a2)
  cf = ch.only_conformer()
  rd = cf.only_residue()
  assert rd.id_str(suppress_segid=1) == 'model="12345678" pdbres="longGLYdlinCh1234J"'
  assert rd.id_str(suppress_segid=-1) == 'model="12345678" pdbres="longGLYdlinCh1234J"'
  try: rd.id_str()
  except ValueError as e:
    assert not show_diff(str(e), '''\
residue.id_str(suppress_segid=false): segid is not unique:
  model="12345678" pdbres="longGLYdlinCh1234J" segid="1234"''')
  else: raise Exception_expected

def exercise_ag_id_str():
  # The function defined in python
  ag = pdb.hierarchy.atom_group()
  ag.altloc='A'
  ag.resname = "ALA"
  assert ag.id_str() == 'AALA       ', "'%s'" % ag.id_str()
  ag.altloc="longAlt"
  ag.resname = "ALAnine"
  assert ag.id_str() == 'longAltALAnine       ', "'%s'" % ag.id_str()

  h = pdb.input(source_info=None, lines="""\
ATOM      6  CA  ASN B   2      -6.522   2.038   2.831  1.00 14.10           C""").construct_hierarchy()
  ag = h.only_model().only_chain().residue_groups()[0].atom_groups()[0]
  assert ag.id_str() == ' ASN B   2 ', "'%s'" % ag.id_str()
  ag.resname='longASN'
  assert ag.id_str() == ' longASN B   2 '
  h.only_model().only_chain().id='dlinB'
  assert ag.id_str() == ' longASNdlinB   2 '
  ag.altloc='A'
  assert ag.id_str() == 'AlongASNdlinB   2 '

def exercise_rg_id_str():
  # The function defined in python
  h = pdb.input(source_info=None, lines="""\
ATOM      6  CA  ASN B   2      -6.522   2.038   2.831  1.00 14.10           C""").construct_hierarchy()
  rg = h.only_model().only_chain().residue_groups()[0]
  assert rg.id_str() == ' B   2 ', "'%s'" % rg.id_str()
  rg.icode='F'
  rg.resseq="A001"
  h.only_model().only_chain().id="dlinB"
  assert rg.id_str() == 'dlinBA001F', "'%s'" % rg.id_str()



if (__name__ == "__main__"):
  exercise_atom_residue_id_str()
  exercise_awl_id_str()
  exercise_ag_id_str()
  exercise_rg_id_str()

  print(format_cpu_times())
