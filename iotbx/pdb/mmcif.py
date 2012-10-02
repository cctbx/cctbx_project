from __future__ import division
from cctbx.array_family import flex
from cctbx import crystal
from libtbx.containers import OrderedSet
from libtbx import group_args
import iotbx.pdb
from iotbx.pdb import hierarchy
from iotbx.pdb import hy36encode
from iotbx.pdb.remark_3_interpretation import \
     refmac_range_to_phenix_string_selection, tls
import iotbx.cif
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
    label_asym_id = cif_block.get("_atom_site.label_asym_id") # chain id
    auth_asym_id = cif_block.get("_atom_site.auth_asym_id")
    if label_asym_id is None: label_asym_id = auth_asym_id
    if auth_asym_id is None: auth_asym_id = label_asym_id
    comp_id = cif_block.get("_atom_site.auth_comp_id")
    if comp_id is None:
      comp_id = cif_block.get("_atom_site.label_comp_id") # residue name
    entity_id = cif_block.get("_atom_site.label_entity_id")
    seq_id = cif_block.get("_atom_site.auth_seq_id")
    if seq_id is None:
      seq_id = cif_block.get("_atom_site.label_seq_id") # residue number
    assert [atom_labels, alt_id, auth_asym_id, comp_id, entity_id, seq_id].count(None) == 0
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
    formal_charge = cif_block.get("_atom_site.pdbx_formal_charge")

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
      if len(unique_model_ids) == 1: i_model = ""
      model = hierarchy.model(id=i_model)
      self.hierarchy.append_model(model)
      unique_chain_ids = OrderedSet(label_asym_id.select(model_sel))
      model.pre_allocate_chains(len(unique_chain_ids))
      for i_chain in unique_chain_ids:
        # we use label_asym_id to identify the separate chains because this
        # separates chains properly in the absence of TER records in the mmcif,
        # however we need to use auth_asym_id for the chain id so that they match
        # e.g. the TLS selections
        chain_sel = (label_asym_id == i_chain) & model_sel
        chain_id = set(auth_asym_id.select(chain_sel))
        assert len(chain_id) == 1
        chain_id = list(chain_id)[0]
        chain = hierarchy.chain(id=chain_id)
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
            if ins_code in ("?", ".", None): ins_code = " "
            residue_group = hierarchy.residue_group(
              resseq=hy36encode(width=4, value=int(i_residue)), icode=ins_code)
            chain.append_residue_group(residue_group)
            unique_altloc_ids = OrderedSet(alt_id.select(ins_code_sel))
            if len(unique_altloc_ids) > 1 and "." in unique_altloc_ids:
              # main chain atoms should appear before altlocs in the hierarchy
              unique_altloc_ids.discard(".")
              unique_altloc_ids = ["."] + list(unique_altloc_ids)
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
                  if formal_charge is not None:
                    charge = formal_charge[i_atom]
                    if charge not in ("?", "."):
                      charge = int(charge)
                      sign = ""
                      if charge > 0: sign = "+"
                      elif charge < 0: sign = "-"
                      atom.set_charge("%i%s" %(abs(charge), sign))
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


class extract_tls_from_cif_block(object):
  def __init__(self, cif_block, pdb_hierarchy):
    self.tls_params = []
    self.tls_present = False
    self.error_string = None
    tls_ids = cif_block.get('_pdbx_refine_tls.id')

    T_ijs = [cif_block.get('_pdbx_refine_tls.T[%s][%s]' %(i, j))
            for i, j in ('11', '22', '33', '12', '13', '23')]
    L_ijs = [cif_block.get('_pdbx_refine_tls.L[%s][%s]' %(i, j))
            for i, j in ('11', '22', '33', '12', '13', '23')]
    S_ijs = [cif_block.get('_pdbx_refine_tls.S[%s][%s]' %(i, j))
            for i, j in ('11', '12', '13', '21', '22', '23', '31', '32', '33')]
    origin_xyzs = [cif_block.get('_pdbx_refine_tls.origin_%s' %x) for x in 'xyz']

    if isinstance(tls_ids, basestring):
      # in case the TLS items are not in a loop
      tls_ids = [tls_ids]
      T_ijs = [[T_ij] for T_ij in T_ijs]
      L_ijs = [[L_ij] for L_ij in L_ijs]
      S_ijs = [[S_ij] for S_ij in S_ijs]
      origin_xyzs = [[origin_xyz] for origin_xyz in origin_xyzs]

    assert T_ijs.count(None) in (0, 6)
    assert L_ijs.count(None) in (0, 6)
    assert S_ijs.count(None) in (0, 9)
    assert origin_xyzs.count(None) in (0, 3)

    if T_ijs.count(None) == 6:
      return

    tls_group_ids = cif_block.get('_pdbx_refine_tls_group.id')
    refine_tls_ids = cif_block.get('_pdbx_refine_tls_group.refine_tls_id')
    beg_chain_ids = cif_block.get('_pdbx_refine_tls_group.beg_auth_asym_id')
    beg_seq_ids = cif_block.get('_pdbx_refine_tls_group.beg_auth_seq_id')
    end_chain_ids = cif_block.get('_pdbx_refine_tls_group.end_auth_asym_id')
    end_seq_ids = cif_block.get('_pdbx_refine_tls_group.end_auth_seq_id')
    selection_details = cif_block.get('_pdbx_refine_tls_group.selection_details')

    assert tls_group_ids is not None

    if isinstance(tls_group_ids, basestring):
      # in case the TLS group items are not in a loop
      refine_tls_ids = flex.std_string([refine_tls_ids])
      beg_chain_ids = [beg_chain_ids]
      beg_seq_ids = [beg_seq_ids]
      end_chain_ids = [end_chain_ids]
      end_seq_ids = [end_seq_ids]
      if selection_details is not None:
        selection_details = flex.std_string([selection_details])

    selection_cache = pdb_hierarchy.atom_selection_cache()

    for i, tls_id in enumerate(tls_ids):
      T = [float(T_ij[i]) for T_ij in T_ijs]
      L = [float(L_ij[i]) for L_ij in L_ijs]
      S = [float(S_ij[i]) for S_ij in S_ijs]
      origin = [float(origin_x[i]) for origin_x in origin_xyzs]
      i_groups = (refine_tls_ids == tls_id).iselection()
      sel_strings = []
      for i_group in i_groups:
        if selection_details is not None and not selection_details.all_eq('?'):
          # phenix selection
          sel_strings.append(selection_details[i_group].strip())
          # check it is a valid selection string
          selection_cache.selection(sel_strings[-1])
        else:
          sel_strings.append(refmac_range_to_phenix_string_selection(
            pdb_hierarchy=pdb_hierarchy,
            chain_start=beg_chain_ids[i_group],
            resseq_start=beg_seq_ids[i_group],
            chain_end=end_chain_ids[i_group],
            resseq_end=end_seq_ids[i_group]))
      sel_str = " or ".join(sel_strings)
      self.tls_params.append(tls(
        t=T, l=L, s=S, origin=origin, selection_string=sel_str))
    self.tls_present = True

  def extract(self):
    return self.tls_params


def tls_as_cif_block(tlsos, selection_strings, cif_block=None):
  import iotbx.cif.model
  if cif_block is None:
    cif_block = iotbx.cif.model.block()

  tls_loop = iotbx.cif.model.loop(header=(
    #"_pdbx_refine_tls.pdbx_refine_id",
    "_pdbx_refine_tls.id",
    "_pdbx_refine_tls.details",
    "_pdbx_refine_tls.method",
    "_pdbx_refine_tls.origin_x",
    "_pdbx_refine_tls.origin_y",
    "_pdbx_refine_tls.origin_z",
    "_pdbx_refine_tls.T[1][1]",
    "_pdbx_refine_tls.T[2][2]",
    "_pdbx_refine_tls.T[3][3]",
    "_pdbx_refine_tls.T[1][2]",
    "_pdbx_refine_tls.T[1][3]",
    "_pdbx_refine_tls.T[2][3]",
    "_pdbx_refine_tls.L[1][1]",
    "_pdbx_refine_tls.L[2][2]",
    "_pdbx_refine_tls.L[3][3]",
    "_pdbx_refine_tls.L[1][2]",
    "_pdbx_refine_tls.L[1][3]",
    "_pdbx_refine_tls.L[2][3]",
    "_pdbx_refine_tls.S[1][1]",
    "_pdbx_refine_tls.S[1][2]",
    "_pdbx_refine_tls.S[1][3]",
    "_pdbx_refine_tls.S[2][1]",
    "_pdbx_refine_tls.S[2][2]",
    "_pdbx_refine_tls.S[2][3]",
    "_pdbx_refine_tls.S[3][1]",
    "_pdbx_refine_tls.S[3][2]",
    "_pdbx_refine_tls.S[3][3]",
  ))

  tls_group_loop = iotbx.cif.model.loop(header=(
    #_pdbx_refine_tls_group.pdbx_refine_id
    "_pdbx_refine_tls_group.id",
    "_pdbx_refine_tls_group.refine_tls_id",
    "_pdbx_refine_tls_group.selection",
    "_pdbx_refine_tls_group.selection_details",
  ))
  for i_tls, (tlso, selection_string) in enumerate(zip(tlsos, selection_strings)):
    tls_row = [i_tls+1, "?", "refined"]
    tls_row.extend(list(tlso.origin))
    tls_row.extend(list(tlso.t))
    tls_row.extend(list(tlso.l))
    tls_row.extend(list(tlso.s))
    tls_loop.add_row(tls_row)
    tls_group_loop.add_row((i_tls+1, i_tls+1, "?", selection_string))

  cif_block.add_loop(tls_loop)
  cif_block.add_loop(tls_group_loop)
  return cif_block


class cif_input(iotbx.pdb.pdb_input_mixin):
  def __init__(self,
               file_name=None,
               source_info=iotbx.pdb.Please_pass_string_or_None,
               lines=None,
               pdb_id=None,
               raise_sorry_if_format_error=False):
    if (pdb_id is not None):
      assert file_name is None
      file_name = iotbx.pdb.ent_path_local_mirror(pdb_id=pdb_id)
    if file_name is not None:
      reader = iotbx.cif.reader(file_path=file_name)
    elif lines is not None:
      reader = iotbx.cif.reader(input_string="\n".join(lines))
    self.cif_model = reader.model()
    self.cif_block = self.cif_model.values()[0]
    self.builder = pdb_hierarchy_builder(reader.model().blocks.values()[0])
    self.hierarchy = self.builder.hierarchy
    self._source_info = "file %s" %file_name

  def construct_hierarchy(self):
    return self.hierarchy

  def source_info(self):
    return self._source_info

  def atoms(self):
    return self.hierarchy.atoms()

  def atoms_with_labels(self):
    return self.hierarchy.atoms_with_labels()

  def model_indices(self):
    return flex.size_t([m.atoms_size() for m in self.hierarchy.models()])

  def crystal_symmetry(self,
        crystal_symmetry=None,
        weak_symmetry=False):
    self_symmetry = self.builder.crystal_symmetry
    if (crystal_symmetry is None):
      return self_symmetry
    if (self_symmetry is None):
      return crystal_symmetry
    return self_symmetry.join_symmetry(
      other_symmetry=crystal_symmetry,
      force=not weak_symmetry)

  def crystal_symmetry_from_cryst1(self):
    return self.builder.crystal_symmetry

  def extract_cryst1_z_columns(self):
    return self.cif_model.values()[0].get("_cell.Z_PDB")

  def connectivity_section(self):
    # XXX Should we extract something from the CIF and return a PDB-like
    # CONECT record, or rework this whole thing?
    return ""

  def connectivity_annotation_section(self):
    return ""

  def remark_section(self):
    return ""

  def heterogen_section(self):
    return ""

  def crystallographic_section(self):
    return ""

  def title_section(self):
    return ""

  def extract_header_year(self):
    yyyymmdd = self.deposition_date()
    if yyyymmdd is not None:
      return int(yyyymmdd[:4])

  def deposition_date(self):
    # date format: yyyy-mm-dd
    cif_block = self.cif_model.values()[0]
    rev_num = cif_block.get('_database_PDB_rev.num')
    if rev_num is not None:
      date_original = cif_block.get('_database_PDB_rev.date_original')
      if isinstance(rev_num, basestring):
        return date_original
      else:
        i = flex.first_index(rev_num, '1')
        if date_original is not None:
          return date_original[i]

  def get_r_rfree_sigma(self, file_name):
    return _cif_get_r_rfree_sigma_object(self.cif_block, file_name)

  def get_solvent_content(self):
    return _float_or_None(self.cif_block.get('_exptl_crystal.density_percent_sol'))

  def get_matthews_coeff(self):
    return _float_or_None(self.cif_block.get('_exptl_crystal.density_Matthews'))

  def get_program_name(self):
    software_name = self.cif_block.get('_software.name')
    software_classification = self.cif_block.get('_software.classification')
    if (isinstance(software_classification, basestring) and
        software_classification == 'refinement'):
      return software_name
    if software_classification is not None:
      i = flex.first_index(software_classification, 'refinement')
      if i > 0: return software_name[i]

  def extract_tls_params(self, hierarchy):
    return extract_tls_from_cif_block(self.cif_block, hierarchy)

  def scale_matrix(self):
    fractionalization_matrix = [
      self.cif_block.get('_atom_sites.fract_transf_matrix[%s][%s]' %(i, j))
      for i,j in ('11', '12', '13', '21', '22', '23', '31','32', '33')]
    if fractionalization_matrix.count(None) == 9:
      return None
    assert fractionalization_matrix.count(None) == 0
    fractionalization_vector = [
      self.cif_block.get('_atom_sites.fract_transf_vector[%s]' %i) for i in '123']
    assert fractionalization_matrix.count(None) == 0
    return [[float(i) for i in fractionalization_matrix],
            [float(i) for i in fractionalization_vector]]

  def model_ids(self):
    return flex.std_string([model.id for model in self.hierarchy.models()])

  def extract_f_model_core_constants(self):
    return extract_f_model_core_constants(self.cif_block)


def _float_or_None(value):
  if value is not None:
    if value == '?' or value == '.':
      return None
    return float(value)


class _cif_get_r_rfree_sigma_object(object):
  def __init__(self, cif_block, file_name):
    self.file_name = file_name
    self.r_work = _float_or_None(cif_block.get('_refine.ls_R_factor_R_work'))
    self.r_free = _float_or_None(cif_block.get('_refine.ls_R_factor_R_free'))
    self.sigma = _float_or_None(cif_block.get('_refine.pdbx_ls_sigma_F'))
    self.high = _float_or_None(cif_block.get('_refine.ls_d_res_high'))
    self.low = _float_or_None(cif_block.get('_refine.ls_d_res_low'))
    self.resolution = _float_or_None(cif_block.get('_reflns.d_resolution_high'))

  def formatted_string(self):
    result = "%s %s %s %s %s %s %s" % (
      format_value("%6s",self.file_name),
      format_value("%6s",str(self.r_work)),
      format_value("%6s",str(self.r_free)),
      format_value("%6s",str(self.sigma)),
      format_value("%6s",str(self.high)),
      format_value("%6s",str(self.low)),
      format_value("%6s",str(self.resolution)))
    return result

  def show(self, log = None):
    if(log is None): log = sys.stdout
    print >> self.log, self.formatted_string()


def extract_f_model_core_constants(cif_block):
  k_sol = _float_or_None(cif_block.get('_refine.solvent_model_param_ksol'))
  b_sol = _float_or_None(cif_block.get('_refine.solvent_model_param_bsol'))
  b_cart = [_float_or_None(cif_block.get('_refine.aniso_B[%s][%s]' %(i, j)))
            for i, j in ('11', '22', '33', '12', '13', '23')]
  assert b_cart.count(None) in (0, 6)
  r_solv = _float_or_None(cif_block.get('_refine.pdbx_solvent_vdw_probe_radii'))
  r_shrink = _float_or_None(cif_block.get('_refine.pdbx_solvent_shrinkage_radii'))
  r_work = _float_or_None(cif_block.get('_refine.ls_R_factor_R_work'))
  r_free = _float_or_None(cif_block.get('_refine.ls_R_factor_R_free'))
  # TODO: extract these from the CIF?
  twin_fraction = None
  twin_law = None
  grid_step_factor = None
  return group_args(
    k_sol            = k_sol,
    b_sol            = b_sol,
    b_cart           = b_cart,
    twin_fraction    = twin_fraction,
    twin_law         = twin_law,
    r_solv           = r_solv,
    r_shrink         = r_shrink,
    grid_step_factor = grid_step_factor,
    r_work           = r_work,
    r_free           = r_free)


class pdb_hierarchy_as_cif_block(iotbx.cif.crystal_symmetry_as_cif_block):

  def __init__(self, pdb_hierarchy, crystal_symmetry=None):
    if crystal_symmetry is None:
      crystal_symmetry = crystal.symmetry()
    iotbx.cif.crystal_symmetry_as_cif_block.__init__(self, crystal_symmetry)

    from string import ascii_uppercase

    atom_site_loop = iotbx.cif.model.loop(header=(
      '_atom_site.group_PDB',
      '_atom_site.id',
      '_atom_site.type_symbol',
      '_atom_site.label_atom_id',
      '_atom_site.label_alt_id',
      '_atom_site.label_comp_id',
      '_atom_site.label_asym_id',
      '_atom_site.label_entity_id',
      '_atom_site.label_seq_id',
      '_atom_site.pdbx_PDB_ins_code',
      '_atom_site.Cartn_x',
      '_atom_site.Cartn_y',
      '_atom_site.Cartn_z',
      '_atom_site.occupancy',
      '_atom_site.B_iso_or_equiv',
      '_atom_site.auth_seq_id',
      '_atom_site.auth_comp_id',
      '_atom_site.auth_asym_id',
      '_atom_site.auth_atom_id',
      '_atom_site.pdbx_PDB_model_num',
    ))

    model_ids = flex.std_string()
    unique_chain_ids = set()
    auth_asym_ids = flex.std_string()
    label_asym_ids = flex.std_string()
    label_seq_id = 0
    for model in pdb_hierarchy.models():
      model_id = model.id
      if model_id == '': model_id = '1'
      for chain in model.chains():
        auth_asym_id = chain.id
        label_asym_id = auth_asym_id # must be unique for a chain
        while True:
          if label_asym_id not in unique_chain_ids:
            unique_chain_ids.add(label_asym_id)
            break
          i = ascii_uppercase.find(label_asym_id)
          if i == len(ascii_uppercase):
            break # XXX what to do?
          else:
            label_asym_id = ascii_uppercase[i+1]
        for residue_group in chain.residue_groups():
          seq_id = residue_group.resseq.strip()
          label_seq_id += 1
          icode = residue_group.icode
          for atom_group in residue_group.atom_groups():
            alt_id = atom_group.altloc
            if alt_id == '': alt_id = '.'
            comp_id = atom_group.resname
            entity_id = '?' # XXX how do we determine this?
            for atom in atom_group.atoms():
              group_pdb = "ATOM"
              if atom.hetero: group_PDB = "HETATM"
              x, y, z = [str(i) for i in atom.xyz]
              atom_site_loop['_atom_site.group_PDB'].append(group_pdb)
              atom_site_loop['_atom_site.id'].append(atom.serial)
              atom_site_loop['_atom_site.type_symbol'].append(atom.element)
              atom_site_loop['_atom_site.label_atom_id'].append(atom.name.strip())
              atom_site_loop['_atom_site.label_alt_id'].append(alt_id)
              atom_site_loop['_atom_site.label_comp_id'].append(comp_id)
              atom_site_loop['_atom_site.label_asym_id'].append(label_asym_id)
              atom_site_loop['_atom_site.label_entity_id'].append(entity_id)
              atom_site_loop['_atom_site.label_seq_id'].append(str(label_seq_id))
              atom_site_loop['_atom_site.pdbx_PDB_ins_code'].append(icode)
              atom_site_loop['_atom_site.Cartn_x'].append(x)
              atom_site_loop['_atom_site.Cartn_y'].append(y)
              atom_site_loop['_atom_site.Cartn_z'].append(z)
              atom_site_loop['_atom_site.occupancy'].append(str(atom.occ))
              atom_site_loop['_atom_site.B_iso_or_equiv'].append(str(atom.b))
              atom_site_loop['_atom_site.auth_seq_id'].append(seq_id)
              atom_site_loop['_atom_site.auth_comp_id'].append(comp_id)
              atom_site_loop['_atom_site.auth_asym_id'].append(auth_asym_id)
              atom_site_loop['_atom_site.auth_atom_id'].append(atom.name.strip())
              atom_site_loop['_atom_site.pdbx_PDB_model_num'].append(model_id)
    self.cif_block.add_loop(atom_site_loop)
