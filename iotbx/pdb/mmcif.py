from __future__ import division
import sys
from cctbx.array_family import flex
from cctbx import crystal
from libtbx.containers import OrderedDict, OrderedSet
from libtbx import dict_with_default_0
from libtbx import group_args
from libtbx.table_utils import wrap_always
from libtbx.utils import null_out, Sorry
from libtbx.str_utils import format_value
import iotbx.pdb
from iotbx.pdb import hierarchy
from iotbx.pdb import hy36encode, hy36decode
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

    atom_site_fp = cif_block.get('_atom_site.phenix_scat_dispersion_real')
    atom_site_fdp = cif_block.get('_atom_site.phenix_scat_dispersion_imag')

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

    current_model_id = None
    current_label_asym_id = None
    current_auth_asym_id = None
    current_residue_id = None
    current_ins_code = None

    for i_atom in range(atom_labels.size()):
      # model(s)
      last_model_id = current_model_id
      current_model_id = model_ids[i_atom]
      assert current_model_id is not None
      if current_model_id != last_model_id:
        model = hierarchy.model(id=current_model_id)
        self.hierarchy.append_model(model)

      # chain(s)
      last_label_asym_id = current_label_asym_id
      current_label_asym_id = label_asym_id[i_atom]
      assert current_label_asym_id is not None
      last_auth_asym_id = current_auth_asym_id
      current_auth_asym_id = auth_asym_id[i_atom]
      if current_auth_asym_id == ".": current_auth_asym_id = " "
      assert current_label_asym_id is not None
      if current_label_asym_id != last_label_asym_id:
        chain = hierarchy.chain(id=current_auth_asym_id)
        model.append_chain(chain)
      else:
        assert current_auth_asym_id == last_auth_asym_id

      # residue_group(s)
      # defined by residue id and insertion code
      last_residue_id = current_residue_id
      current_residue_id = seq_id[i_atom]
      assert current_residue_id is not None
      last_ins_code = current_ins_code
      if pdb_ins_code is not None:
        current_ins_code = pdb_ins_code[i_atom]
        if current_ins_code in ("?", ".", None): current_ins_code = " "
      if (current_residue_id != last_residue_id or
          current_ins_code != last_ins_code or
          current_label_asym_id != last_label_asym_id):
        residue_group = hierarchy.residue_group(
          resseq=hy36encode(width=4, value=int(current_residue_id)),
          icode=current_ins_code)
        chain.append_residue_group(residue_group)
        atom_groups = OrderedDict() # reset atom_groups cache

      # atom_group(s)
      # defined by resname and altloc id
      current_altloc = alt_id[i_atom]
      if current_altloc == ".": current_altloc = "" # Main chain atoms
      current_resname = comp_id[i_atom]
      if (current_altloc, current_resname) not in atom_groups:
        atom_group = hierarchy.atom_group(
          altloc=current_altloc, resname=current_resname)
        atom_groups[(current_altloc, current_resname)] = atom_group
        if current_altloc == "":
          residue_group.insert_atom_group(0, atom_group)
        else:
          residue_group.append_atom_group(atom_group)
      else:
        atom_group = atom_groups[(current_altloc, current_resname)]

      # atom(s)
      atom = hierarchy.atom()
      atom_group.append_atom(atom)
      atom.set_element(type_symbol[i_atom])
      atom.set_name(
        format_pdb_atom_name(atom_labels[i_atom], type_symbol[i_atom]))
      atom.set_xyz(
        new_xyz=(cart_x[i_atom], cart_y[i_atom], cart_z[i_atom]))
      atom.set_b(B_iso_or_equiv[i_atom])
      atom.set_occ(occu[i_atom])
      # hy36encode should go once the pdb.hierarchy has been
      # modified to no longer store fixed-width strings
      atom.set_serial(
        hy36encode(width=5, value=int(atom_site_id[i_atom])))
      # some code relies on an empty segid being 4 spaces
      atom.set_segid("    ")
      if group_PDB is not None and group_PDB[i_atom] == "HETATM":
        atom.hetero = True
      if formal_charge is not None:
        charge = formal_charge[i_atom]
        if charge not in ("?", "."):
          if charge.endswith("-") or charge.startswith("-"):
            sign = "-"
          else:
            sign = "+"
          charge = charge.strip(" -+")
          charge = int(charge)
          if charge == 0: sign = ""
          atom.set_charge("%i%s" %(charge, sign))
      if atom_site_fp is not None:
        fp = atom_site_fp[i_atom]
        if fp not in ("?", "."):
          atom.set_fp(new_fp=float(fp))
      if atom_site_fdp is not None:
        fdp = atom_site_fdp[i_atom]
        if fdp not in ("?", "."):
          atom.set_fdp(new_fdp=float(fdp))
      if anisotrop_id is not None and adps is not None:
        u_ij_index = flex.first_index(anisotrop_id, atom.serial.strip())
        if u_ij_index is not None:
          u_ij = adps[u_ij_index]
          atom.set_uij(u_ij)
        else:
          pass

    if len(self.hierarchy.models()) == 1:
      # for compatibility with single-model PDB files
      self.hierarchy.models()[0].id = ""


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

  if (len(selection_strings) == 0):
    assert len(tlsos) == 1
    selection_strings = ["all"]
  else:
    assert len(tlsos) == len(selection_strings)

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
               cif_object=None,
               source_info=iotbx.pdb.Please_pass_string_or_None,
               lines=None,
               pdb_id=None,
               raise_sorry_if_format_error=False):
    if (pdb_id is not None):
      assert file_name is None
      file_name = iotbx.pdb.ent_path_local_mirror(pdb_id=pdb_id)
    if file_name is not None:
      reader = iotbx.cif.reader(file_path=file_name)
      self.cif_model = reader.model()
    elif lines is not None:
      reader = iotbx.cif.reader(input_string="\n".join(lines))
      self.cif_model = reader.model()
    elif cif_object is not None:
      self.cif_model = cif_object
    if len(self.cif_model) == 0:
      raise Sorry("mmCIF file must contain at least one data block")
    self.cif_block = self.cif_model.values()[0]
    self.builder = pdb_hierarchy_builder(self.cif_block)
    self.hierarchy = self.builder.hierarchy
    self._source_info = "file %s" %file_name

  def file_type (self) :
    return "mmcif"

  def construct_hierarchy(self, set_atom_i_seq=True):
    if (set_atom_i_seq) :
      self.hierarchy.reset_atom_i_seqs()
    return self.hierarchy

  def source_info(self):
    return self._source_info

  def atoms(self):
    return self.hierarchy.atoms()

  def atoms_with_labels(self):
    return self.hierarchy.atoms_with_labels()

  def model_indices(self):
    return flex.size_t([m.atoms_size() for m in self.hierarchy.models()])

  def ter_indices(self):
    # for compatibility with pdb_input
    return flex.size_t()

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
      if i >= 0: return software_name[i]

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

  def extract_wavelength (self, first_only=True) :
    wavelengths = self.cif_block.get('_diffrn_source.pdbx_wavelength_list')
    wavelength = _float_or_None(self.cif_block.get(
        '_diffrn_source.pdbx_wavelength'))
    if (first_only) :
      if (wavelength is not None) :
        return wavelength
      elif (wavelengths is not None) :
        wavelengths = [ float(f.strip()) for f in wavelengths.split(",") ]
        return wavelengths[0]
    elif (wavelengths is not None) :
        return [ float(f) for f in wavelengths.split(",") ]
    return None

  def get_experiment_type (self) :
    exptl_method = self.cif_block.get('_exptl.method')
    return exptl_method

  def process_mtrix_records(self,error_handle=True,eps=1e-4):
    """
    Read MTRIX records from a pdb file

    Arguments:
    ----------
    error_handle -- True: will stop execution on improper rotation matrices
                    False: will continue execution but will replace the values in the
                           rotation matrix with [0,0,0,0,0,0,0,0,0]
    eps -- Rounding accuracy for avoiding numerical issue when when testing proper rotation

    error_handle can be use if one does not want a 'Sorry' to be raised. The program will continue to execute
    but all coordinates of the bad rotation matrix will zero out

    Returns:
    --------
    result : object containing information on all NCS operations
    """
    from scitbx import matrix

    trans = []
    rots = []
    serial_number = []
    coordinates_present = []

    for i,sn in enumerate(self.cif_block.get('_struct_ncs_oper.id')):
      serial_number.append((sn,i))
      coordinates_present.append(
        self.cif_block.get('_struct_ncs_oper.code')[i] == 'given')
      r = [(self.cif_block.get('_struct_ncs_oper.matrix[%s][%s]' %(x,y))[i])
        for x,y in ('11', '12', '13', '21', '22', '23', '31','32', '33')]
      rots.append(matrix.sqr(map(float,r)))
      t = [(self.cif_block.get('_struct_ncs_oper.vector[%s]' %x)[i])
        for x in '123']
      trans.append(matrix.col(map(float,t)))
    # sort records by serial number
    serial_number.sort()
    items_order = [i for (_,i) in serial_number]
    trans = [trans[i] for i in items_order]
    rots = [rots[i] for i in items_order]
    coordinates_present = [coordinates_present[i] for i in items_order]
    serial_number = [j for (j,_) in serial_number]
    # collect results
    result = iotbx.pdb.mtrix_and_biomt_records_container()
    for sn,r,t,cp in zip(serial_number,rots,trans,coordinates_present):
      if not r.is_r3_rotation_matrix(rms_tolerance=eps):
        if error_handle:
          raise Sorry('Rotation matrices are not proper! ')
        else:
          print Sorry('Rotation matrices are not proper! ')
      ignore_transform = r.is_r3_identity_matrix() and t.is_col_zero()
      result.add(
        r=r, t=t,
        coordinates_present=(cp or ignore_transform),
        serial_number=sn)
    return result

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
    print >> log, self.formatted_string()


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

  def __init__(self, pdb_hierarchy, crystal_symmetry=None,
               coordinate_precision=5,
               occupancy_precision=3,
               b_iso_precision=5,
               u_aniso_precision=5):
    if crystal_symmetry is None:
      crystal_symmetry = crystal.symmetry()
    iotbx.cif.crystal_symmetry_as_cif_block.__init__(
      self, crystal_symmetry, format="mmcif")

    coord_fmt_str = "%%.%if" %coordinate_precision
    occ_fmt_str = "%%.%if" %occupancy_precision
    b_iso_fmt_str = "%%.%if" %b_iso_precision
    u_aniso_fmt_str = "%%.%if" %u_aniso_precision

    atom_site_loop = iotbx.cif.model.loop(header=(
      '_atom_site.group_PDB',
      '_atom_site.id',
      '_atom_site.label_atom_id',
      '_atom_site.label_alt_id',
      '_atom_site.label_comp_id',
      '_atom_site.auth_asym_id',
      '_atom_site.auth_seq_id',
      '_atom_site.pdbx_PDB_ins_code',
      '_atom_site.Cartn_x',
      '_atom_site.Cartn_y',
      '_atom_site.Cartn_z',
      '_atom_site.occupancy',
      '_atom_site.B_iso_or_equiv',
      '_atom_site.type_symbol',
      '_atom_site.pdbx_formal_charge',
      '_atom_site.phenix_scat_dispersion_real',
      '_atom_site.phenix_scat_dispersion_imag',
      '_atom_site.label_asym_id',
      '_atom_site.label_entity_id',
      '_atom_site.label_seq_id',
      #'_atom_site.auth_comp_id',
      #'_atom_site.auth_atom_id',
      '_atom_site.pdbx_PDB_model_num',
    ))

    aniso_loop = iotbx.cif.model.loop(header=(
      '_atom_site_anisotrop.id',
      '_atom_site_anisotrop.pdbx_auth_atom_id',
      '_atom_site_anisotrop.pdbx_label_alt_id',
      '_atom_site_anisotrop.pdbx_auth_comp_id',
      '_atom_site_anisotrop.pdbx_auth_asym_id',
      '_atom_site_anisotrop.pdbx_auth_seq_id',
      '_atom_site_anisotrop.pdbx_PDB_ins_code',
      '_atom_site_anisotrop.U[1][1]',
      '_atom_site_anisotrop.U[2][2]',
      '_atom_site_anisotrop.U[3][3]',
      '_atom_site_anisotrop.U[1][2]',
      '_atom_site_anisotrop.U[1][3]',
      '_atom_site_anisotrop.U[2][3]'
    ))

    # cache dictionary lookups to save time in inner loop
    atom_site_group_PDB = atom_site_loop['_atom_site.group_PDB']
    atom_site_id = atom_site_loop['_atom_site.id']
    atom_site_label_atom_id = atom_site_loop['_atom_site.label_atom_id']
    atom_site_label_alt_id = atom_site_loop['_atom_site.label_alt_id']
    atom_site_label_comp_id = atom_site_loop['_atom_site.label_comp_id']
    atom_site_auth_asym_id = atom_site_loop['_atom_site.auth_asym_id']
    atom_site_auth_seq_id = atom_site_loop['_atom_site.auth_seq_id']
    atom_site_pdbx_PDB_ins_code = atom_site_loop['_atom_site.pdbx_PDB_ins_code']
    atom_site_Cartn_x = atom_site_loop['_atom_site.Cartn_x']
    atom_site_Cartn_y = atom_site_loop['_atom_site.Cartn_y']
    atom_site_Cartn_z = atom_site_loop['_atom_site.Cartn_z']
    atom_site_occupancy = atom_site_loop['_atom_site.occupancy']
    atom_site_B_iso_or_equiv = atom_site_loop['_atom_site.B_iso_or_equiv']
    atom_site_type_symbol = atom_site_loop['_atom_site.type_symbol']
    atom_site_pdbx_formal_charge = atom_site_loop['_atom_site.pdbx_formal_charge']
    atom_site_phenix_scat_dispersion_real = \
      atom_site_loop['_atom_site.phenix_scat_dispersion_real']
    atom_site_phenix_scat_dispersion_imag = \
      atom_site_loop['_atom_site.phenix_scat_dispersion_imag']
    atom_site_label_asym_id = atom_site_loop['_atom_site.label_asym_id']
    atom_site_label_entity_id = atom_site_loop['_atom_site.label_entity_id']
    atom_site_label_seq_id = atom_site_loop['_atom_site.label_seq_id']
    #atom_site_loop['_atom_site.auth_comp_id'].append(comp_id)
    #atom_site_loop['_atom_site.auth_atom_id'].append(atom.name.strip())
    atom_site_pdbx_PDB_model_num = atom_site_loop['_atom_site.pdbx_PDB_model_num']
    atom_site_anisotrop_id = aniso_loop['_atom_site_anisotrop.id']
    atom_site_anisotrop_pdbx_auth_atom_id = \
      aniso_loop['_atom_site_anisotrop.pdbx_auth_atom_id']
    atom_site_anisotrop_pdbx_label_alt_id = \
      aniso_loop['_atom_site_anisotrop.pdbx_label_alt_id']
    atom_site_anisotrop_pdbx_auth_comp_id = \
      aniso_loop['_atom_site_anisotrop.pdbx_auth_comp_id']
    atom_site_anisotrop_pdbx_auth_asym_id = \
      aniso_loop['_atom_site_anisotrop.pdbx_auth_asym_id']
    atom_site_anisotrop_pdbx_auth_seq_id = \
      aniso_loop['_atom_site_anisotrop.pdbx_auth_seq_id']
    atom_site_anisotrop_pdbx_PDB_ins_code = \
      aniso_loop['_atom_site_anisotrop.pdbx_PDB_ins_code']
    atom_site_anisotrop_U11 = aniso_loop['_atom_site_anisotrop.U[1][1]']
    atom_site_anisotrop_U22 = aniso_loop['_atom_site_anisotrop.U[2][2]']
    atom_site_anisotrop_U33 = aniso_loop['_atom_site_anisotrop.U[3][3]']
    atom_site_anisotrop_U12 = aniso_loop['_atom_site_anisotrop.U[1][2]']
    atom_site_anisotrop_U13 = aniso_loop['_atom_site_anisotrop.U[1][3]']
    atom_site_anisotrop_U23 = aniso_loop['_atom_site_anisotrop.U[2][3]']

    model_ids = flex.std_string()
    unique_chain_ids = set()
    auth_asym_ids = flex.std_string()
    label_asym_ids = flex.std_string()
    label_seq_id = 0
    label_asym_id = ""
    for model in pdb_hierarchy.models():
      model_id = model.id
      if model_id == '': model_id = '1'
      for chain in model.chains():
        auth_asym_id = chain.id
        if auth_asym_id.strip() == '': auth_asym_id = '.'
        label_asym_id = increment_label_asym_id(label_asym_id)
        for residue_group in chain.residue_groups():
          seq_id = residue_group.resseq.strip()
          label_seq_id += 1
          icode = residue_group.icode
          if icode == ' ': icode = '?'
          for atom_group in residue_group.atom_groups():
            alt_id = atom_group.altloc
            if alt_id == '': alt_id = '.'
            comp_id = atom_group.resname.strip()
            entity_id = '?' # XXX how do we determine this?
            for atom in atom_group.atoms():
              group_pdb = "ATOM"
              if atom.hetero: group_pdb = "HETATM"
              x, y, z = [coord_fmt_str %i for i in atom.xyz]
              atom_charge = atom.charge_tidy()
              if atom_charge is None:
                atom_charge = "?"
              else:
                atom_charge = atom_charge.strip()
              if atom_charge == "": atom_charge = "?"
              fp, fdp = atom.fp, atom.fdp
              if fp == 0 and fdp == 0:
                fp = '.'
                fdp = '.'
              else:
                fp = "%.4f" %fp
                fdp = "%.4f" %fdp
              atom_site_group_PDB.append(group_pdb)
              atom_site_id.append(str(hy36decode(width=5, s=atom.serial)))
              atom_site_label_atom_id.append(atom.name.strip())
              atom_site_label_alt_id.append(alt_id)
              atom_site_label_comp_id.append(comp_id)
              atom_site_auth_asym_id.append(auth_asym_id)
              atom_site_auth_seq_id.append(seq_id)
              atom_site_pdbx_PDB_ins_code.append(icode)
              atom_site_Cartn_x.append(x)
              atom_site_Cartn_y.append(y)
              atom_site_Cartn_z.append(z)
              atom_site_occupancy.append(occ_fmt_str % atom.occ)
              atom_site_B_iso_or_equiv.append(b_iso_fmt_str % atom.b)
              atom_site_type_symbol.append(atom.element.strip())
              atom_site_pdbx_formal_charge.append(atom_charge)
              atom_site_phenix_scat_dispersion_real.append(fp)
              atom_site_phenix_scat_dispersion_imag.append(fdp)
              atom_site_label_asym_id.append(label_asym_id)
              atom_site_label_entity_id.append(entity_id)
              atom_site_label_seq_id.append(str(label_seq_id))
              #atom_site_loop['_atom_site.auth_comp_id'].append(comp_id)
              #atom_site_loop['_atom_site.auth_atom_id'].append(atom.name.strip())
              atom_site_pdbx_PDB_model_num.append(model_id.strip())

              if atom.uij_is_defined():
                u11, u22, u33, u12, u13, u23 = [
                  u_aniso_fmt_str %i for i in atom.uij]
                atom_site_anisotrop_id.append(atom.serial.strip())
                atom_site_anisotrop_pdbx_auth_atom_id.append(atom.name.strip())
                atom_site_anisotrop_pdbx_label_alt_id.append(alt_id)
                atom_site_anisotrop_pdbx_auth_comp_id.append(comp_id)
                atom_site_anisotrop_pdbx_auth_asym_id.append(auth_asym_id)
                atom_site_anisotrop_pdbx_auth_seq_id.append(seq_id)
                atom_site_anisotrop_pdbx_PDB_ins_code.append(icode)
                atom_site_anisotrop_U11.append(u11)
                atom_site_anisotrop_U22.append(u22)
                atom_site_anisotrop_U33.append(u33)
                atom_site_anisotrop_U12.append(u12)
                atom_site_anisotrop_U13.append(u13)
                atom_site_anisotrop_U23.append(u23)

    for key in ('_atom_site.phenix_scat_dispersion_real',
                '_atom_site.phenix_scat_dispersion_imag'):
      if atom_site_loop[key].all_eq('.'):
        del atom_site_loop[key]
    self.cif_block.add_loop(atom_site_loop)
    if aniso_loop.size() > 0:
      self.cif_block.add_loop(aniso_loop)


class pdb_hierarchy_as_cif_block_with_sequence(pdb_hierarchy_as_cif_block):
  def __init__(self, pdb_hierarchy,
               sequences,
               alignment_params=None,
               crystal_symmetry=None,
               coordinate_precision=5,
               occupancy_precision=3,
               b_iso_precision=5,
               u_aniso_precision=5):

    pdb_hierarchy_as_cif_block.__init__(
      self, pdb_hierarchy, crystal_symmetry=crystal_symmetry,
    coordinate_precision=coordinate_precision,
    occupancy_precision=occupancy_precision,
    b_iso_precision=b_iso_precision,
    u_aniso_precision=u_aniso_precision)

    import mmtbx.validation.sequence
    validation = mmtbx.validation.sequence.validation(
      pdb_hierarchy=pdb_hierarchy,
      sequences=sequences,
      params=alignment_params,
      extract_residue_groups=True,
      log=null_out(), # silence output
    )

    entity_loop = iotbx.cif.model.loop(header=(
      '_entity.id',
      '_entity.type',
      #'_entity.src_method',
      #'_entity.pdbx_description',
      '_entity.formula_weight',
      '_entity.pdbx_number_of_molecules',
      #'_entity.details',
      #'_entity.pdbx_mutation',
      #'_entity.pdbx_fragment',
      #'_entity.pdbx_ec'
    ))

    entity_poly_loop = iotbx.cif.model.loop(header=(
      '_entity_poly.entity_id',
      '_entity_poly.type',
      '_entity_poly.nstd_chirality',
      '_entity_poly.nstd_linkage',
      '_entity_poly.nstd_monomer',
      '_entity_poly.pdbx_seq_one_letter_code',
      '_entity_poly.pdbx_seq_one_letter_code_can',
      '_entity_poly.pdbx_strand_id',
      '_entity_poly.type_details'
    ))

    entity_poly_seq_loop = iotbx.cif.model.loop(header=(
      '_entity_poly_seq.entity_id',
      '_entity_poly_seq.num',
      '_entity_poly_seq.mon_id',
      '_entity_poly_seq.hetero',
    ))

    sequence_counts = OrderedDict()
    sequence_to_chain_ids = {}
    entity_id = 0
    sequence_to_entity_id = {}
    chain_id_to_entity_id = {}
    sequence_to_chains = {}
    residue_group_to_seq_num_mapping = {}
    aligned_pdb_chains = OrderedSet()
    non_polymer_counts = dict_with_default_0()
    non_polymer_resname_to_entity_id = OrderedDict()

    for chain in validation.chains:
      sequence = chain.alignment.b
      if sequence not in sequence_to_entity_id:
        entity_id += 1
        sequence_to_entity_id[sequence] = entity_id
      sequence_counts.setdefault(sequence, 0)
      sequence_counts[sequence] += 1
      sequence_to_chain_ids.setdefault(sequence, [])
      sequence_to_chain_ids[sequence].append(chain.chain_id)
      sequence_to_chains.setdefault(sequence, [])
      sequence_to_chains[sequence].append(chain)
      chain_id_to_entity_id[chain.chain_id] = sequence_to_entity_id[sequence]
      aligned_pdb_chains.add(chain.residue_groups[0].parent())
      unaligned_pdb_chains = OrderedSet(pdb_hierarchy.chains()) - aligned_pdb_chains

      assert len(chain.residue_groups) + chain.n_missing_start + chain.n_missing_end == len(sequence)
      residue_groups = [None] * chain.n_missing_start + chain.residue_groups + [None] * chain.n_missing_end
      i = chain.n_missing_start
      seq_num = 0
      for i, residue_group in enumerate(residue_groups):
        if residue_group is None and chain.alignment.b[i] == '-':
          # a deletion
          continue
        seq_num += 1
        if residue_group is not None:
          residue_group_to_seq_num_mapping[
            residue_group] = seq_num

    for pdb_chain in unaligned_pdb_chains:
      for residue_group in pdb_chain.residue_groups():
        for resname in residue_group.unique_resnames():
          if resname not in non_polymer_resname_to_entity_id:
            entity_id += 1
            non_polymer_resname_to_entity_id[resname] = entity_id
          non_polymer_counts[resname] += 1

    for sequence, count in sequence_counts.iteritems():
      entity_poly_seq_num = 0
      entity_id = sequence_to_entity_id[sequence]

      entity_loop.add_row((
        entity_id,
        'polymer', #polymer/non-polymer/macrolide/water
        #'?', #src_method
        #'?', # pdbx_description
        '?', # formula_weight
        len(sequence_to_chains[sequence]), # pdbx_number_of_molecules
        #'?', # details
        #'?', # pdbx_mutation
        #'?', # pdbx_fragment
        #'?' # pdbx_ec
      ))

      # The definition of the cif item _entity_poly.pdbx_seq_one_letter_code
      # says that modifications and non-standard amino acids should be encoded
      # as 'X', however in practice the PDB seem to encode them as the three-letter
      # code in parentheses.
      pdbx_seq_one_letter_code = []
      pdbx_seq_one_letter_code_can = []

      chains = sequence_to_chains[sequence]

      from iotbx.pdb import amino_acid_codes

      chain = chains[0]
      matches = chain.alignment.matches()

      for i, one_letter_code in enumerate(sequence):

        #Data items in the ENTITY_POLY_SEQ category specify the sequence
        #of monomers in a polymer. Allowance is made for the possibility
        #of microheterogeneity in a sample by allowing a given sequence
        #number to be correlated with more than one monomer ID. The
        #corresponding ATOM_SITE entries should reflect this
        #heterogeneity.

        monomer_id = None
        if i >= chain.n_missing_start and i < (len(sequence) - chain.n_missing_end):
          monomer_id = chain.resnames[i-chain.n_missing_start]

        if monomer_id is None and one_letter_code == '-': continue

        pdbx_seq_one_letter_code_can.append(one_letter_code)

        if monomer_id is None:
          if sequence_to_chains[sequence][0].chain_type == mmtbx.validation.sequence.PROTEIN:
            monomer_id = amino_acid_codes.three_letter_given_one_letter.get(
              one_letter_code, "UNK") # XXX
          else:
            monomer_id = one_letter_code
        else:
          if sequence_to_chains[sequence][0].chain_type == mmtbx.validation.sequence.PROTEIN:
            one_letter_code = amino_acid_codes.one_letter_given_three_letter.get(
              monomer_id, "(%s)" %monomer_id)

        pdbx_seq_one_letter_code.append(one_letter_code)

        entity_poly_seq_num += 1

        entity_poly_seq_loop.add_row((
          entity_id,
          entity_poly_seq_num,
          monomer_id,
          'no', #XXX
        ))

      entity_poly_type = '?'
      entity_nstd_chirality = 'n'
      # we should probably determine the chirality more correctly by examining
      # the chirality of the backbone chain rather than relying on the residue
      # names to be correct
      if chain.chain_type == mmtbx.validation.sequence.PROTEIN:
        n_d_peptides = 0
        n_l_peptides = 0
        n_achiral_peptides = 0
        n_unknown = 0
        for resname in chain.resnames:
          if resname == "GLY":
            n_achiral_peptides += 1
          elif resname in iotbx.pdb.common_residue_names_amino_acid:
            n_l_peptides += 1
          elif resname in amino_acid_codes.three_letter_l_given_three_letter_d:
            n_d_peptides += 1
          else:
            n_unknown += 1
        n_total = sum([n_d_peptides, n_l_peptides, n_achiral_peptides, n_unknown])
        if (n_l_peptides + n_achiral_peptides)/n_total > 0.5:
          entity_poly_type = 'polypeptide(L)'
          if n_d_peptides > 0:
            entity_nstd_chirality = 'y'
        elif (n_d_peptides + n_achiral_peptides)/n_total > 0.5:
          entity_poly_type = 'polypeptide(D)'
          if n_l_peptides > 0:
            entity_nstd_chirality = 'y'
      elif chain.chain_type == mmtbx.validation.sequence.NUCLEIC_ACID:
        n_dna = 0
        n_rna = 0
        n_unknown = 0
        for resname in chain.resnames:
          if resname is not None and resname.strip().upper() in (
            'AD', 'CD', 'GD', 'TD', 'DA', 'DC', 'DG', 'DT'):
            n_dna += 1
          elif resname is not None and resname.strip().upper() in (
            'A', 'C', 'G', 'T', '+A', '+C', '+G', '+T'):
            n_rna += 1
          else:
            n_unknown += 1
        n_total = sum([n_dna + n_rna + n_unknown])
        if n_dna/n_total > 0.5 and n_rna == 0:
          entity_poly_type = 'polydeoxyribonucleotide'
        elif n_rna/n_total > 0.5 and n_dna == 0:
          entity_poly_type = 'polyribonucleotide'
        elif (n_rna + n_dna)/n_total > 0.5:
          entity_poly_type = 'polydeoxyribonucleotide/polyribonucleotide hybrid'

      entity_poly_loop.add_row((
        entity_id,
        entity_poly_type,
        entity_nstd_chirality,
        'no',
        'no',
        wrap_always("".join(pdbx_seq_one_letter_code), width=80).strip(),
        wrap_always("".join(pdbx_seq_one_letter_code_can), width=80).strip(),
        ','.join(sequence_to_chain_ids[sequence]),
        '?'
      ))

    for resname, entity_id in non_polymer_resname_to_entity_id.iteritems():
      entity_type = "non-polymer"
      if resname == "HOH":
        entity_type = "water" # XXX
      entity_loop.add_row((
        entity_id,
        entity_type, #polymer/non-polymer/macrolide/water
        #'?', #src_method
        #'?', # pdbx_description
        '?', # formula_weight
        non_polymer_counts[resname], # pdbx_number_of_molecules
        #'?', # details
        #'?', # pdbx_mutation
        #'?', # pdbx_fragment
        #'?' # pdbx_ec
      ))

    self.cif_block.add_loop(entity_loop)
    self.cif_block.add_loop(entity_poly_loop)
    self.cif_block.add_loop(entity_poly_seq_loop)
    self.cif_block.update(pdb_hierarchy.as_cif_block())

    label_entity_id = self.cif_block['_atom_site.label_entity_id']
    auth_seq_id = self.cif_block['_atom_site.auth_seq_id']
    ins_code = self.cif_block['_atom_site.pdbx_PDB_ins_code']
    auth_asym_id = self.cif_block['_atom_site.auth_asym_id']
    label_seq_id = flex.std_string(auth_seq_id.size(), '.')
    ins_code = ins_code.deep_copy()
    ins_code.set_selected(ins_code == '?', '')
    for residue_group, seq_num in residue_group_to_seq_num_mapping.iteritems():
      sel = ((auth_asym_id == residue_group.parent().id) &
             (ins_code == residue_group.icode.strip()) &
             (auth_seq_id == residue_group.resseq.strip()))
      label_seq_id.set_selected(sel, str(seq_num))
      label_entity_id.set_selected(
        sel, str(chain_id_to_entity_id[residue_group.parent().id]))

    for pdb_chain in unaligned_pdb_chains:
      for residue_group in pdb_chain.residue_groups():
        sel = ((auth_asym_id == residue_group.parent().id) &
               (ins_code == residue_group.icode.strip()) &
               (auth_seq_id == residue_group.resseq.strip()))
        label_entity_id.set_selected(
          sel, str(non_polymer_resname_to_entity_id[residue_group.unique_resnames()[0]]))

    self.cif_block['_atom_site.label_seq_id'] = label_seq_id

    # reorder the loops
    atom_site_loop = self.cif_block['_atom_site']
    atom_site_aniso_loop = self.cif_block.get('_atom_site_anisotrop')
    del self.cif_block['_atom_site']
    self.cif_block.add_loop(atom_site_loop)
    if atom_site_aniso_loop is not None:
      del self.cif_block['_atom_site_anisotrop']
      self.cif_block.add_loop(atom_site_aniso_loop)


def increment_label_asym_id(asym_id):
  from string import ascii_uppercase
  if len(asym_id) == 0:
    return "A"
  asym_id = list(asym_id)
  for i in range(len(asym_id)):
    if asym_id[i] == "Z":
      asym_id[i] = "A"
      if (i+1) == len(asym_id):
        return "A" * (len(asym_id) + 1)
    else:
      while True:
        j = ascii_uppercase.find(asym_id[i])
        asym_id[i] = ascii_uppercase[j+1]
        return "".join(asym_id)
