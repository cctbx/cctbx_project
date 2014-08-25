mmtbx.command_line package
==========================

The ``command_line`` directory contains many common utilities as well as
several relatively complex applications, including **Xtriage**, the core
of **MolProbity**, **phenix.maps**, and **phenix.table_one**.  Many of these
tools have in common the same basic inputs and setup procedures.  These are
encapsulated in the top-level mmtbx.command_line module, and simplify the
process of loading models and data to just a few lines of code for each new
application.

API documentation
-----------------

.. automodule:: mmtbx.command_line
    :members:
    :undoc-members:
    :show-inheritance:

A more complex example, which identifies suspect ligands based on electron
density:

  >>> master_phil = mmtbx.command_line.generate_master_phil_with_inputs(
  ...   enable_twin_law=True,
  ...   phil_string="""
  ...     hetatms_only = True
  ...       .type = bool
  ...     skip_single_atoms = True
  ...       .type = bool
  ...     min_acceptable_cc = 0.8
  ...       .type = float""")
  ...
  >>> cmdline = mmtbx.command_line.load_model_and_data(
  ...  args=sys.argv[1:],
  ...  master_phil=master_phil,
  ...  out=sys.stdout,
  ...  process_pdb_file=True,
  ...  create_fmodel=True,
  ...  prefer_anomalous=False)
  ...
  >>> log = cmdline.start_log_file("bad_ligands.log")
  >>> params = cmdline.params
  >>> bad_ligands = mmtbx.real_space_correlation.find_suspicious_residues(
  ...   fmodel=fmodel,
  ...   pdb_hierarchy=pdb_hierarchy,
  ...   hetatms_only=params.hetatms_only,
  ...   skip_single_atoms=params.skip_single_atoms,
  ...   min_acceptable_cc=params.min_acceptable_cc,
  ...   log=log)
  ...
