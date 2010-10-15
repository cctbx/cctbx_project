import libtbx

class connectivity_table(object):

  covalent_bond_tolerance = 0.5 # Angstrom
  conformer_indices=None
  sym_excl_indices=None

  # TODO: add possibility to fine tune connectivity by addition and deletion
  #       of individual bonds, and also by user-provided radii

  def __init__(self,
               structure,
               **kwds):
    from cctbx.eltbx import covalent_radii
    from cctbx import crystal
    self.structure = structure
    libtbx.adopt_optional_init_args(self, kwds)

    radii = [
      covalent_radii.table(elt).radius() for elt in
      structure.scattering_type_registry().type_index_pairs_as_dict() ]
    buffer_thickness = 2*max(radii) + self.covalent_bond_tolerance
    asu_mappings = structure.asu_mappings(
      buffer_thickness=buffer_thickness)
    self.pair_asu_table = crystal.pair_asu_table(asu_mappings)
    self.pair_asu_table.add_covalent_pairs(
      structure.scattering_types(),
      conformer_indices=self.conformer_indices,
      sym_excl_indices=self.sym_excl_indices,
      tolerance=self.covalent_bond_tolerance)
