"""reduce2 neutron tests: a returned model passed straight into a second run"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
from scitbx import matrix
from mmtbx.regression import reduce2_neutron_utils
from mmtbx.regression.reduce2_neutron_utils import (
  FLOATING, bgm_atoms, check_water, heavy_pdb, pdb_1yjp, run_prog,
  run_prog_on_model)

reduce2_neutron_utils.set_prefix("tst_r2n_in_memory")

def test_in_memory_rerun():
  # get_results().model passed straight into a second run gives the same
  # result as the written and reread file: its monomer library server
  # already knows BGM (CCD only), its riding-H manager has the first run's
  # X-H lengths, a scattering table set on it would set the water O-H,
  # refinement flags would restrict the idealization; the second run starts
  # from none of them
  bgm = heavy_pdb("BGM", bgm_atoms(), offset=(15.0, 10.0, 15.0),
    cryst1=False).replace(" A   1 ", " B   1 ")
  systems = {"protein": pdb_1yjp.replace("END\n", "") + bgm + "END\n",
             "alone": heavy_pdb("BGM", bgm_atoms())}
  def atom_lines(file_name):
    with open(file_name) as f:
      return [l for l in f if l.startswith(("ATOM", "HETATM"))]
  for system, first, second in (("protein", "none", "hd_and_h"),
      ("alone", "none", "hd_and_h"), ("protein", "none", "none"),
      ("protein", "hd_and_h", "hd_and_h")):
    prefix = "tst_r2n_in_memory_inmem_%s_%s_%s" % (system, first, second)
    prog = run_prog(prefix + "_1", systems[system], neutron_option=first)
    with open(prog.params.output.filename) as f: written = f.read()
    model = prog.get_results().model
    prog = run_prog_on_model(prefix + "_2", model, approach="optimize",
      neutron_option=second)
    with open(prog.params.output.description_file_name) as f:
      description = f.read()
    assert FLOATING not in description, (system, first, second)
    hd = [l for l in atom_lines(prog.params.output.filename)
          if l[17:20] == "BGM" and l[76:78].strip() in ("H", "D")]
    n_d = len([l for l in hd if l[76:78].strip() == "D"])
    if second == "none":
      assert (len(hd), n_d) == (13, 0), len(hd)
    else:
      assert (len(hd), n_d) == (19, 6), (len(hd), n_d)
    reread = run_prog(prefix + "_file", written, approach="optimize",
      neutron_option=second)
    assert atom_lines(prog.params.output.filename) == \
      atom_lines(reread.params.output.filename), (system, first, second)
  # no crystal symmetry in: the returned model has none (not reduce2's box),
  # so the rerun writes the same file as the file rerun, without CRYST1; the
  # cryo-EM placeholder in: both write it
  no_symmetry = "".join(l for l in pdb_1yjp.splitlines(True)
                        if not l.startswith("CRYST1"))
  placeholder = "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 " \
    "P 1           1\n" + no_symmetry
  for case, option in (("nosym", "none"), ("nosym", "hd_and_h"),
      ("placeholder", "hd_and_h")):
    prefix = "tst_r2n_in_memory_inmem_%s_%s" % (case, option)
    prog = run_prog(prefix + "_1", no_symmetry if case == "nosym" else
      placeholder, neutron_option=option)
    with open(prog.params.output.filename) as f: written = f.read()
    prog = run_prog_on_model(prefix + "_2", prog.get_results().model,
      approach="optimize", neutron_option=option)
    reread = run_prog(prefix + "_file", written, approach="optimize",
      neutron_option=option)
    with open(prog.params.output.filename) as f: in_memory = f.read()
    with open(reread.params.output.filename) as f: from_file = f.read()
    cryst1 = [l.split()[1:4] for l in in_memory.splitlines()
              if l.startswith("CRYST1")]
    assert cryst1 == ([] if case == "nosym" else
      [["1.000", "1.000", "1.000"]]), (case, option, cryst1)
    assert in_memory == from_file, (case, option)
  # X-H lengths follow the second run: X-ray first, neutron distances second
  # (the first run's riding-H manager is not reused); refinement flags on the
  # model (here: no site refined) do not restrict the riding-H idealization
  # (either would leave the X-ray C-H)
  from mmtbx.refinement import refinement_flags
  from scitbx.array_family import flex
  prefix = "tst_r2n_in_memory_inmem_neutron_distances"
  prog = run_prog(prefix + "_1", systems["protein"])
  with open(prog.params.output.filename) as f: written = f.read()
  def c_h(file_name):
    h = iotbx.pdb.input(file_name).construct_hierarchy()
    return sorted(set(round((matrix.col(a.xyz) -
      matrix.col(ag.get_atom("CA").xyz)).length(), 2) for ag in h.atom_groups()
      for a in ag.atoms() if a.name.strip() in ("HA", "HA2", "HA3")))
  assert c_h(prog.params.output.filename) == [0.97]
  model = prog.get_results().model
  model.set_refinement_flags(refinement_flags.manager(
    sites_individual=flex.bool(model.get_number_of_atoms(), False)))
  prog = run_prog_on_model(prefix + "_2", model,
    approach="optimize", use_neutron_distances=True)
  assert c_h(prog.params.output.filename) == [1.09]
  reread = run_prog(prefix + "_file", written, approach="optimize",
    use_neutron_distances=True)
  assert atom_lines(prog.params.output.filename) == \
    atom_lines(reread.params.output.filename)
  # a scattering table set on the model (neutron: water O-H 0.98 from the
  # table) is not used: X-ray run, water O-H 0.85 as from the file
  prefix = "tst_r2n_in_memory_inmem_table"
  prog = run_prog(prefix + "_1", pdb_1yjp, exclude_water=False)
  with open(prog.params.output.filename) as f: written = f.read()
  model = prog.get_results().model
  model.setup_scattering_dictionaries(scattering_table="neutron")
  prog = run_prog_on_model(prefix + "_2", model, approach="optimize",
    exclude_water=False)
  check_water(prog.params.output.filename, 0.85, 7)
  reread = run_prog(prefix + "_file", written, approach="optimize",
    exclude_water=False)
  assert atom_lines(prog.params.output.filename) == \
    atom_lines(reread.params.output.filename)

if __name__ == "__main__":
  test_in_memory_rerun()
  print("OK")
