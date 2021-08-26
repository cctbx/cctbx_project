from __future__ import absolute_import, division, print_function
from mmtbx.secondary_structure.build import ss_idealization as ssb
import iotbx.pdb
from libtbx.utils import Sorry
import mmtbx.model

def exercise_presence_of_h():
  tst_pdb_1 = """\
HELIX    2   2 GLY A  100  TRP A  101  1                                   2
ATOM    866  N   GLY A 100    -128.903-139.939-104.460  1.00 27.58           N
ATOM    867  CA  GLY A 100    -128.223-140.890-103.603  1.00 27.58           C
ATOM    868  C   GLY A 100    -128.381-142.356-103.959  1.00 27.58           C
ATOM    869  O   GLY A 100    -129.283-142.755-104.686  1.00 27.58           O
ATOM    870  H   GLY A 100    -129.862-139.777-104.337  1.00 36.09           H
ATOM    871  N   TRP A 101    -127.491-143.163-103.397  1.00 25.83           N
ATOM    872  CA  TRP A 101    -127.455-144.601-103.619  1.00 25.83           C
ATOM    873  C   TRP A 101    -128.824-145.274-103.752  1.00 25.83           C
ATOM    874  O   TRP A 101    -129.091-145.962-104.734  1.00 25.83           O
ATOM    875  CB  TRP A 101    -126.659-145.240-102.480  1.00 42.04           C
ATOM    876  CG  TRP A 101    -126.207-146.635-102.729  1.00 42.04           C
ATOM    877  CD1 TRP A 101    -126.180-147.300-103.922  1.00 42.04           C
ATOM    878  CD2 TRP A 101    -125.677-147.528-101.760  1.00 42.04           C
ATOM    879  NE1 TRP A 101    -125.663-148.557-103.751  1.00 42.04           N
ATOM    880  CE2 TRP A 101    -125.347-148.723-102.429  1.00 42.04           C
ATOM    881  CE3 TRP A 101    -125.446-147.437-100.384  1.00 42.04           C
ATOM    882  CZ2 TRP A 101    -124.799-149.820-101.767  1.00 42.04           C
ATOM    883  CZ3 TRP A 101    -124.899-148.528 -99.725  1.00 42.04           C
ATOM    884  CH2 TRP A 101    -124.581-149.703-100.415  1.00 42.04           C
ATOM    885  H   TRP A 101    -126.863-142.797-102.745  1.00 42.04           H
"""
  pdb_inp = iotbx.pdb.input(source_info=None, lines=tst_pdb_1)
  ann = pdb_inp.extract_secondary_structure()
  model = mmtbx.model.manager(
    model_input = pdb_inp)
  model.process(make_restraints=True)
  model.set_ss_annotation(ann)
  rm = ssb.substitute_ss(
      model,
      )

def exercise_ca_absent():
  # should raise Sorry
  """ CA from 101 residue are absent. We don't handle it at the moment."""
  tst_pdb_1 = """\
HELIX    2   2 GLY A  100  TRP A  101  1                                   2
ATOM    866  N   GLY A 100    -128.903-139.939-104.460  1.00 27.58           N
ATOM    867  CA  GLY A 100    -128.223-140.890-103.603  1.00 27.58           C
ATOM    868  C   GLY A 100    -128.381-142.356-103.959  1.00 27.58           C
ATOM    869  O   GLY A 100    -129.283-142.755-104.686  1.00 27.58           O
ATOM    870  H   GLY A 100    -129.862-139.777-104.337  1.00 36.09           H
ATOM    871  N   TRP A 101    -127.491-143.163-103.397  1.00 25.83           N
ATOM    873  C   TRP A 101    -128.824-145.274-103.752  1.00 25.83           C
ATOM    874  O   TRP A 101    -129.091-145.962-104.734  1.00 25.83           O
ATOM    875  CB  TRP A 101    -126.659-145.240-102.480  1.00 42.04           C
ATOM    876  CG  TRP A 101    -126.207-146.635-102.729  1.00 42.04           C
ATOM    877  CD1 TRP A 101    -126.180-147.300-103.922  1.00 42.04           C
ATOM    878  CD2 TRP A 101    -125.677-147.528-101.760  1.00 42.04           C
ATOM    879  NE1 TRP A 101    -125.663-148.557-103.751  1.00 42.04           N
ATOM    880  CE2 TRP A 101    -125.347-148.723-102.429  1.00 42.04           C
ATOM    881  CE3 TRP A 101    -125.446-147.437-100.384  1.00 42.04           C
ATOM    882  CZ2 TRP A 101    -124.799-149.820-101.767  1.00 42.04           C
ATOM    883  CZ3 TRP A 101    -124.899-148.528 -99.725  1.00 42.04           C
ATOM    884  CH2 TRP A 101    -124.581-149.703-100.415  1.00 42.04           C
ATOM    885  H   TRP A 101    -126.863-142.797-102.745  1.00 42.04           H
"""
  pdb_inp = iotbx.pdb.input(source_info=None, lines=tst_pdb_1)
  ann = pdb_inp.extract_secondary_structure()
  model = mmtbx.model.manager(
    model_input = pdb_inp)
  model.process(make_restraints=True)
  model.set_ss_annotation(ann)
  try:
    rm = ssb.substitute_ss(
        model)
    rm.run()
  except Sorry as e:
    assert e.args[0].startswith("C, CA or N")
  except Exception:
    assert 0, "Error: This test failed"
  else:
    assert 0, "Error: This test failed"

def exercise_some_atoms_are_absent():
  """ CD1, CD2 from 101 residue are absent"""
  tst_pdb_2 = """\
HELIX    2   2 GLY A  100  TRP A  101  1                                   2
ATOM    866  N   GLY A 100    -128.903-139.939-104.460  1.00 27.58           N
ATOM    867  CA  GLY A 100    -128.223-140.890-103.603  1.00 27.58           C
ATOM    868  C   GLY A 100    -128.381-142.356-103.959  1.00 27.58           C
ATOM    869  O   GLY A 100    -129.283-142.755-104.686  1.00 27.58           O
ATOM    871  N   TRP A 101    -127.491-143.163-103.397  1.00 25.83           N
ATOM    872  CA  TRP A 101    -127.455-144.601-103.619  1.00 25.83           C
ATOM    873  C   TRP A 101    -128.824-145.274-103.752  1.00 25.83           C
ATOM    874  O   TRP A 101    -129.091-145.962-104.734  1.00 25.83           O
ATOM    875  CB  TRP A 101    -126.659-145.240-102.480  1.00 42.04           C
ATOM    876  CG  TRP A 101    -126.207-146.635-102.729  1.00 42.04           C
ATOM    879  NE1 TRP A 101    -125.663-148.557-103.751  1.00 42.04           N
ATOM    880  CE2 TRP A 101    -125.347-148.723-102.429  1.00 42.04           C
ATOM    881  CE3 TRP A 101    -125.446-147.437-100.384  1.00 42.04           C
"""
  pdb_inp = iotbx.pdb.input(source_info=None, lines=tst_pdb_2)
  ann = pdb_inp.extract_secondary_structure()
  model = mmtbx.model.manager(
    model_input = pdb_inp)
  model.process(make_restraints=True)
  model.set_ss_annotation(ann)
  rm = ssb.substitute_ss(
      model)
  rm.run()


def exercise():
  exercise_presence_of_h()
  exercise_ca_absent()
  exercise_some_atoms_are_absent()
  print("OK")

if (__name__ == "__main__"):
  exercise()
