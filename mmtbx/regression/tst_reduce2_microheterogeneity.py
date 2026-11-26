"""H placement at a microheterogeneity site where an atom without altloc is
shared by residues of different names: 1ejg Pro A 22 / Ser B, C 22 with one
blank N, labelled PRO. No non-H atom is added, removed, moved or relabelled
(qmready=True gives every conformer its own N instead,
tst_reduce2_altloc_completion.py).

- The H of Pro go into its altloc A group. As read from a PDB file the Pro is
  a residue group of its own (blank N, A), and its H went into the blank
  group: H of Pro in the Ser conformers.
- Each Ser conformer gets its amide H on the shared N. pdb_interpretation
  does not give the Ser the N (tst_pdb_interpretation_microheterogeneity.py),
  so the H had no bond and was dropped.

Both layouts of the site: two residue groups (PDB file) and one (mmCIF).
"""
from __future__ import absolute_import, division, print_function
import os
from six.moves import cStringIO as StringIO
import iotbx.pdb
import mmtbx.model
from iotbx.cli_parser import run_program
from libtbx.utils import null_out
from scitbx.matrix import col
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.programs import reduce2
from mmtbx.regression.tst_reduce2_pass_environment import m_1ejg_22

def merged(hierarchy):
  """The mmCIF layout: the residue groups of one position merged"""
  for chain in hierarchy.chains():
    primary = None
    for rg in chain.residue_groups():
      if primary is not None and primary.resid() == rg.resid():
        chain.merge_residue_groups(primary=primary, secondary=rg)
      else:
        primary = rg
  hierarchy.atoms().reset_serial()
  return hierarchy

def heavy(hierarchy):
  """The non-H atoms: label (with altloc and residue name), site, occupancy
  and B, in an order of their own (placement sorts the atoms)"""
  return sorted((a.pdb_label_columns(), tuple(round(x, 3) for x in a.xyz),
                 round(a.occ, 2), round(a.b, 2))
                for a in hierarchy.atoms() if not a.element_is_hydrogen())

def site_22(hierarchy):
  """{(altloc, resname): {atom name: atom}} of residue 22"""
  result = {}
  for rg in hierarchy.residue_groups():
    if rg.resseq.strip() != "22": continue
    for ag in rg.atom_groups():
      d = result.setdefault((ag.altloc, ag.resname), {})
      for a in ag.atoms(): d[a.name.strip()] = a
  return result

def angle(a, b, c):
  return (col(a.xyz) - col(b.xyz)).angle(col(c.xyz) - col(b.xyz), deg=True)

def check_site(hierarchy, heavy_before):
  assert heavy(hierarchy) == heavy_before
  site = site_22(hierarchy)
  assert sorted(site.keys()) == [
    ("", "PRO"), ("A", "PRO"), ("B", "SER"), ("C", "SER")], sorted(site.keys())
  # the blank group holds the shared N only: no H of Pro in the Ser conformers
  assert sorted(site[("", "PRO")].keys()) == ["N"], site[("", "PRO")].keys()
  n = site[("", "PRO")]["N"]
  pro_h = sorted(k for k, a in site[("A", "PRO")].items()
                 if a.element_is_hydrogen())
  assert pro_h == ["HA", "HB2", "HB3", "HD2", "HD3", "HG2", "HG3"], pro_h
  for a in site[("A", "PRO")].values():
    if a.element_is_hydrogen(): assert abs(a.occ - 0.57) < 1.e-6, a.occ
  c_prev = [a for a in hierarchy.atoms() if a.name.strip() == "C"
            and a.parent().parent().resseq.strip() == "21"]
  assert len(c_prev) == 1
  for altloc, occ in (("B", 0.21), ("C", 0.22)):
    ser = site[(altloc, "SER")]
    ser_h = sorted(k for k, a in ser.items() if a.element_is_hydrogen())
    assert ser_h == ["H", "HA", "HB2", "HB3", "HG"], ser_h
    h = ser["H"]
    assert abs(h.occ - occ) < 1.e-6, h.occ
    # on the shared N, in the peptide plane, between C of 21 and CA
    assert 0.8 < h.distance(n) < 1.1, h.distance(n)
    a_c, a_ca = angle(c_prev[0], n, h), angle(ser["CA"], n, h)
    assert 110 < a_c < 130, a_c
    assert 105 < a_ca < 130, a_ca
    total = a_c + a_ca + angle(c_prev[0], n, ser["CA"])
    assert total > 350, total

def place(hierarchy):
  cs = iotbx.pdb.input(source_info=None,
    lines=m_1ejg_22.splitlines()).crystal_symmetry()
  model = mmtbx.model.manager(model_input=None, pdb_hierarchy=hierarchy,
    crystal_symmetry=cs, log=null_out())
  o = reduce_hydrogen.place_hydrogens(model=model)
  o.run()
  return o

def input_hierarchy():
  return iotbx.pdb.input(source_info=None,
    lines=m_1ejg_22.splitlines()).construct_hierarchy()

def test_placement(merge):
  h = input_hierarchy()
  assert len([rg for rg in h.residue_groups()
              if rg.resseq.strip() == "22"]) == 2
  if merge: h = merged(h)
  assert len([rg for rg in h.residue_groups()
              if rg.resseq.strip() == "22"]) == (1 if merge else 2)
  before = heavy(h)
  o = place(h)
  check_site(o.get_model().get_hierarchy(), before)
  # nothing of residue 22 is reported as not placed (was: SER A 22 conformer
  # B and C, H on N: parent atom not in the model / no bond to N)
  assert not [r for r in o.not_placed if r.resseq.strip() == "22"], \
    [r.__dict__ for r in o.not_placed if r.resseq.strip() == "22"]

def test_program(suffix):
  """the whole program: the H survive optimization and are written"""
  prefix = "tst_reduce2_microheterogeneity"
  fn_in, fn_out = prefix + suffix, prefix + "_out" + suffix
  h = input_hierarchy()
  cs = iotbx.pdb.input(source_info=None,
    lines=m_1ejg_22.splitlines()).crystal_symmetry()
  with open(fn_in, "w") as f:
    if suffix == ".pdb":
      f.write(m_1ejg_22)
    else:
      f.write(h.as_mmcif_string(crystal_symmetry=cs))
  before = heavy(iotbx.pdb.input(fn_in).construct_hierarchy())
  if os.path.exists(fn_out): os.remove(fn_out)
  logger = StringIO()
  run_program(program_class=reduce2.Program, logger=logger,
    args=[fn_in, "overwrite=True", "output.filename=" + fn_out,
          "output.description_file_name=" + prefix + "_out.txt"])
  out = iotbx.pdb.input(fn_out).construct_hierarchy()
  if suffix == ".cif":
    assert len([rg for rg in out.residue_groups()
                if rg.resseq.strip() == "22"]) == 1
  check_site(out, before)

def run_all():
  test_placement(merge=False)
  test_placement(merge=True)
  test_program(".pdb")
  test_program(".cif")

if __name__ == "__main__":
  run_all()
  print("OK")
