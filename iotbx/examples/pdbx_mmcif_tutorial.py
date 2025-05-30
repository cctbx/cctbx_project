"""Examples of how to work with PDB and mmCIF files"""
from __future__ import absolute_import, division, print_function
#
# Command to run this example:
#   iotbx.python pdbx_mmcif_tutorial.py
#
# See also:
#   http://cctbx.sourceforge.net/iotbx_cif
#
import os
from six.moves import range
from libtbx.utils import Sorry
from iotbx.pdb.fetch import valid_pdb_id, fetch_and_write

def run(args):
  if len(args) == 0:
    args = ["1hbb"]

  for arg in args:
    if os.path.isfile(arg):
      mmcif_file = arg
      pdb_id = os.path.splitext(os.path.basename(mmcif_file))[0]
      if not valid_pdb_id(pdb_id):
        raise Sorry("Not valid pdb id")
    else:
      # download pdbx/mmcif file from the PDB
      pdb_id = arg
      mirror = "pdbe"
      mmcif_file = fetch_and_write(
        pdb_id, entity="model_cif", mirror=mirror, log=sys.stdout)

    # read the cif file and get an iotbx.cif object
    import iotbx.cif
    cif_reader = iotbx.cif.reader(file_path=mmcif_file)
    cif_object = cif_reader.model()
    cif_block = cif_object[pdb_id]
    # get single items from cif_block
    print("PDB id:", cif_block["_entry.id"])
    # get a looped item from cif_block
    print("Authors:")
    for author in cif_block.get_looped_item("_citation_author.name"):
      print(author)
    print()
    print("Molecular Entities:")
    for pdbx_entity in cif_block.get_looped_item("_entity.pdbx_description"):
      print(pdbx_entity)
    print()

    # extract crystal symmetry information
    import iotbx.cif.builders
    builder = iotbx.cif.builders.crystal_symmetry_builder(cif_block)
    builder.crystal_symmetry.show_summary()

    # 1) this works also for .pdb files, but re-reads the file
    import iotbx.pdb
    pdb_input = iotbx.pdb.input(file_name=mmcif_file)
    hierarchy = pdb_input.construct_hierarchy()

    # 2) This only works for mmcif files, but re-uses the cif_object from above:
    import iotbx.pdb.mmcif
    pdb_input = iotbx.pdb.mmcif.cif_input(cif_object=cif_object)
    hierarchy = pdb_input.construct_hierarchy()

    # some convenience methods of pdb_input object
    print("Software:", pdb_input.get_program_name())
    print("Experiment type:", pdb_input.get_experiment_type())
    print("Solvent content:", pdb_input.get_solvent_content())
    print("Deposition date:", pdb_input.deposition_date())
    r_rfree_sigma = pdb_input.get_r_rfree_sigma(mmcif_file)
    print("R-work/R-free: %s/%s" %(r_rfree_sigma.r_work, r_rfree_sigma.r_free))
    # can also get crystal_symmetry from pdb_input object
    crystal_symmetry = pdb_input.crystal_symmetry()

    print()
    hierarchy.overall_counts().show()
    # level_id can be "model", "chain", "residue_group", "atom_group" or "atom"
    hierarchy.show(level_id="chain")
    # for a more detailed example of interacting with a pdb.hierarchy object,
    # see iotbx/examples/pdb_hierarchy.py

    # extract atom sites
    atoms = hierarchy.atoms()
    sites_cart = atoms.extract_xyz()
    print()
    for i in range(10):
      print(atoms[i].id_str(), atoms[i].xyz)
    print()

    # read some sequence information
    entity_poly_entity_id = cif_block.get_looped_item("_entity_poly.entity_id")
    entity_id = cif_block.get_looped_item("_entity.id")
    entity_pdbx_description = cif_block.get_looped_item("_entity.pdbx_description")
    entity_poly_one_letter_code = cif_block.get_looped_item(
      "_entity_poly.pdbx_seq_one_letter_code")

    from cctbx.array_family import flex
    for i in range(len(entity_poly_one_letter_code)):
      idx = flex.first_index(entity_id, entity_poly_entity_id[i])
      print(entity_id[idx], entity_pdbx_description[i], end=' ')
      print("".join(entity_poly_one_letter_code[i].split()))


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
