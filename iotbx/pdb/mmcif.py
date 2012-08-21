from __future__ import division
from cctbx.array_family import flex
from libtbx.containers import OrderedSet
from iotbx.pdb import hierarchy
from iotbx.pdb import hy36encode
from iotbx.cif.builders import crystal_symmetry_builder


class pdb_hierarchy_builder(crystal_symmetry_builder):

  # The recommended translation for ATOM records can be found at:
  #   http://mmcif.rcsb.org/dictionaries/pdb-correspondence/pdb2mmcif-2010.html#ATOM

  def __init__(self, cif_block):
    crystal_symmetry_builder.__init__(self, cif_block)

    self.hierarchy = hierarchy.root()

    # These items are mandatory for the _atom_site loop, all others are optional
    type_symbol = cif_block.get("_atom_site.type_symbol")
    atom_labels = cif_block.get("_atom_site.auth_atom_id")
    if atom_labels is None:
      atom_labels = cif_block.get("_atom_site.label_atom_id") # corresponds to chem comp atom name
    alt_id = cif_block.get("_atom_site.label_alt_id") # alternate conformer id
    asym_id = cif_block.get("_atom_site.label_asym_id") # chain id
    if asym_id is None:
      asym_id = cif_block.get("_atom_site.auth_asym_id")
    comp_id = cif_block.get("_atom_site.auth_comp_id")
    if comp_id is None:
      comp_id = cif_block.get("_atom_site.label_comp_id") # residue name
    entity_id = cif_block.get("_atom_site.label_entity_id")
    seq_id = cif_block.get("_atom_site.auth_seq_id")
    if seq_id is None:
      seq_id = cif_block.get("_atom_site.label_seq_id") # residue number
    assert [atom_labels, alt_id, asym_id, comp_id, entity_id, seq_id].count(None) == 0
    assert type_symbol is not None

    pdb_ins_code = cif_block.get("_atom_site.pdbx_PDB_ins_code") # insertion code
    model_ids = cif_block.get("_atom_site.pdbx_PDB_model_num")
    atom_site_id = cif_block.get("_atom_site.id")
    # only permitted values are ATOM or HETATM
    group_PDB = cif_block.get("_atom_site.group_PDB")

    # TODO: read esds
    B_iso_or_equiv = flex.double(cif_block.get("_atom_site.B_iso_or_equiv"))
    cart_x = flex.double(cif_block.get("_atom_site.Cartn_x"))
    cart_y = flex.double(cif_block.get("_atom_site.Cartn_y"))
    cart_z = flex.double(cif_block.get("_atom_site.Cartn_z"))
    occu = flex.double(cif_block.get("_atom_site.occupancy"))
    # TODO: read charge

    # anisotropic b-factors
    # TODO: read esds
    anisotrop_id = cif_block.get("_atom_site_anisotrop.id")
    adps = None
    if anisotrop_id is not None:
      u_ij = [cif_block.get("_atom_site_anisotrop.U[%s][%s]" %(ij[0], ij[1]))
              for ij in ("11", "22", "33", "12", "13", "23")]
      assert u_ij.count(None) in (0, 6)
      if u_ij.count(None) == 0:
        adps = u_ij
      else:
        assert u_ij.count(None) == 6
        b_ij = [cif_block.get("_atom_site_anisotrop.B[%s][%s]" %(ij[0], ij[1]))
                for ij in ("11", "22", "33", "12", "13", "23")]
        assert b_ij.count(None) in (0, 6)
        if b_ij.count(None) == 0:
          adps = adptbx.b_as_u(b_ij)
        assert not (u_ij.count(None) and b_ij.count(None)) # illegal for both to be present
      if adps is not None:
        try:
          adps = [flex.double(adp) for adp in adps]
        except ValueError, e:
          raise CifBuilderError("Error interpreting ADPs: " + str(e))
        adps = flex.sym_mat3_double(*adps)

    # XXX What if _atom_site.pdbx_PDB_model_num is not given?
    unique_model_ids = OrderedSet(model_ids) # XXX more efficient way to do this?
    self.hierarchy.pre_allocate_models(len(unique_model_ids))
    for i_model in unique_model_ids:
      model_sel = (model_ids == i_model)
      model = hierarchy.model(id=i_model)
      self.hierarchy.append_model(model)
      unique_chain_ids = OrderedSet(asym_id.select(model_sel))
      model.pre_allocate_chains(len(unique_chain_ids))
      for i_chain in unique_chain_ids:
        chain_sel = (asym_id == i_chain) & model_sel
        chain = hierarchy.chain(id=i_chain)
        model.append_chain(chain)
        unique_residue_ids = OrderedSet(seq_id.select(chain_sel))
        chain.pre_allocate_residue_groups(len(unique_residue_ids))
        # XXX do we need to sort the residue ids, or leave them in the order we found them?
        for i_residue in unique_residue_ids:
          residue_sel = (seq_id == i_residue) & chain_sel
          if pdb_ins_code is not None:
            ins_codes = pdb_ins_code.select(residue_sel)
            unique_pdb_ins_codes = OrderedSet(ins_codes)
          else:
            unique_pdb_ins_codes = [None]
          for ins_code in unique_pdb_ins_codes:
            if ins_code is not None:
              ins_code_sel = (ins_code == pdb_ins_code) & residue_sel
            else:
              ins_code_sel = residue_sel
            if ins_code in ("?", "."): ins_code = None
            residue_group = hierarchy.residue_group(
              resseq=hy36encode(width=4, value=int(i_residue)), icode=ins_code)
            chain.append_residue_group(residue_group)
            unique_altloc_ids = OrderedSet(alt_id.select(ins_code_sel))
            residue_group.pre_allocate_atom_groups(len(unique_altloc_ids))
            for i_altloc in unique_altloc_ids:
              atom_group_sel = (alt_id == i_altloc) & ins_code_sel
              resnames = comp_id.select(atom_group_sel)
              unique_resnames = OrderedSet(resnames)
              # by this point there should be only one resname left
              assert len(unique_resnames) == 1 # should all in the atom group have the same resname?
              for resname in unique_resnames:
                resname_sel = (comp_id == resname) & atom_group_sel
                if i_altloc == ".": i_altloc = "" # Main chain atoms
                atom_group = hierarchy.atom_group(altloc=i_altloc, resname=resname)
                residue_group.append_atom_group(atom_group)
                atom_group_isel = atom_group_sel.iselection()
                atom_group.pre_allocate_atoms(len(atom_group_isel))
                for i_atom in atom_group_isel:
                  atom = hierarchy.atom()
                  atom_group.append_atom(atom)
                  atom.set_element(type_symbol[i_atom])
                  atom.set_name(
                    format_pdb_atom_name(atom_labels[i_atom], type_symbol[i_atom]))
                  atom.set_xyz(
                    new_xyz=(cart_x[i_atom], cart_y[i_atom], cart_z[i_atom]))
                  atom.set_b(B_iso_or_equiv[i_atom])
                  atom.set_occ(occu[i_atom])
                  atom.set_serial(atom_site_id[i_atom])
                  if anisotrop_id is not None and adps is not None:
                    u_ij_index = flex.first_index(anisotrop_id, atom.serial)
                    if u_ij_index is not None:
                      u_ij = adps[u_ij_index]
                      atom.set_uij(u_ij)
                    else:
                      pass

def format_pdb_atom_name(atom_name, atom_type):
  # The PDB-format atom name is 4 characters long (columns 13 - 16):
  #   "Alignment of one-letter atom name such as C starts at column 14,
  #    while two-letter atom name such as FE starts at column 13."
  # http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
  # Here we assume that "atom name" refers to the atom/element TYPE.

  if len(atom_name) > 4:
    # XXX will bad things happen if the atom name is more than 4 characters?
    return atom_name
  elif len(atom_name) == 4:
    return atom_name
  atom_name = atom_name.strip()
  if len(atom_type) == 1:
    if atom_name.upper()[0] == atom_type.upper():
      atom_name = " %s" %atom_name
  elif len(atom_type) == 2:
    if atom_name.upper()[0] == atom_type.upper():
      atom_name = atom_name
  while len(atom_name) < 4:
    atom_name = "%s " %atom_name
  return atom_name
