"""Interpretation of mmCIF formatted model files"""
from __future__ import absolute_import, division, print_function
import sys
from cctbx.array_family import flex
from libtbx.containers import OrderedDict
from libtbx import group_args
from libtbx.utils import Sorry
from libtbx.str_utils import format_value
import iotbx.pdb
import cctbx.crystal
from iotbx.pdb import hierarchy
from iotbx.pdb import hy36encode
from iotbx.pdb import cryst1_interpretation
from iotbx.pdb.experiment_type import experiment_type
from iotbx.pdb.remark_3_interpretation import \
     refmac_range_to_phenix_string_selection, tls
import iotbx.cif
from iotbx.cif.builders import crystal_symmetry_builder
import iotbx.mtrix_biomt
from six import string_types
from six.moves import range, zip

class pdb_hierarchy_builder(crystal_symmetry_builder):

  # The recommended translation for ATOM records can be found at:
  #   http://mmcif.rcsb.org/dictionaries/pdb-correspondence/pdb2mmcif-2010.html#ATOM
  #   https://www.ebi.ac.uk/pdbe/docs/exchange/pdb-correspondence/pdb2mmcif.html#ATOM

  #  Note: pdb_hierarchy allows the same chain ID to be used in multiple chains.
  #    This is indicated in PDB format with a TER record indicating the end of a chain
  #    The corresponding information in mmCIF is the label_asym_id changes for each
  #    new chain (the auth_asym_id matches the chain ID).
  #    Catch cases where one chain id is split into multiple chains in the following way:
  #    If label_asym_id changes and previous residue was protein or rna/dna or modified
  #    (or auth_asym_id or model_id changes), create a chain break.

  def __init__(self, cif_block):
    crystal_symmetry_builder.__init__(self, cif_block)

    self.hierarchy = hierarchy.root()
    # These items are mandatory for the _atom_site loop, all others are optional
    type_symbol = self._wrap_loop_if_needed(cif_block, "_atom_site.type_symbol")
    atom_labels = self._wrap_loop_if_needed(cif_block, "_atom_site.label_atom_id") # corresponds to chem comp atom name
    if atom_labels is None:
      atom_labels = self._wrap_loop_if_needed(cif_block, "_atom_site.auth_atom_id")
    alt_id = self._wrap_loop_if_needed(cif_block,"_atom_site.label_alt_id") # alternate conformer id
    label_asym_id = self._wrap_loop_if_needed(cif_block, "_atom_site.label_asym_id") # chain id
    auth_asym_id = self._wrap_loop_if_needed(cif_block, "_atom_site.auth_asym_id")
    auth_segid = self._wrap_loop_if_needed(cif_block, "_atom_site.auth_segid")
    auth_break = self._wrap_loop_if_needed(cif_block, "_atom_site.auth_break")
    if label_asym_id is None: label_asym_id = auth_asym_id
    if auth_asym_id is None: auth_asym_id = label_asym_id
    comp_id = self._wrap_loop_if_needed(cif_block, "_atom_site.auth_comp_id")
    if comp_id is None:
      comp_id = self._wrap_loop_if_needed(cif_block, "_atom_site.label_comp_id") # residue name
    entity_id = self._wrap_loop_if_needed(cif_block, "_atom_site.label_entity_id")
    seq_id = self._wrap_loop_if_needed(cif_block, "_atom_site.auth_seq_id")
    if seq_id is None:
      seq_id = self._wrap_loop_if_needed(cif_block, "_atom_site.label_seq_id") # residue number
    error_message = """something is not present - not enough information
    to make a hierarcy out of this CIF file.
    It could be because this is restraints or component cif or a corrupted model cif."""
    assert [atom_labels, alt_id, auth_asym_id, comp_id, entity_id, seq_id].count(None) == 0, error_message
    assert type_symbol is not None

    atom_site_fp = cif_block.get('_atom_site.phenix_scat_dispersion_real')
    atom_site_fdp = cif_block.get('_atom_site.phenix_scat_dispersion_imag')

    pdb_ins_code = cif_block.get("_atom_site.pdbx_PDB_ins_code") # insertion code
    model_ids = cif_block.get("_atom_site.pdbx_PDB_model_num")
    atom_site_id = cif_block.get("_atom_site.id")
    # only permitted values are ATOM or HETATM
    group_PDB = cif_block.get("_atom_site.group_PDB")
    # TODO: read esds
    B_iso_or_equiv = flex.double(self._wrap_loop_if_needed(cif_block, "_atom_site.B_iso_or_equiv"))
    cart_x = flex.double(self._wrap_loop_if_needed(cif_block, "_atom_site.Cartn_x"))
    cart_y = flex.double(self._wrap_loop_if_needed(cif_block, "_atom_site.Cartn_y"))
    cart_z = flex.double(self._wrap_loop_if_needed(cif_block, "_atom_site.Cartn_z"))
    occu =   flex.double(self._wrap_loop_if_needed(cif_block, "_atom_site.occupancy"))
    formal_charge = self._wrap_loop_if_needed(cif_block, "_atom_site.pdbx_formal_charge")
    # anisotropic b-factors
    # TODO: read esds
    anisotrop_id = self._wrap_loop_if_needed(cif_block, "_atom_site_anisotrop.id")
    adps = None
    if anisotrop_id is not None:
      u_ij = [self._wrap_loop_if_needed(cif_block, "_atom_site_anisotrop.U[%s][%s]" %(ij[0], ij[1]))
              for ij in ("11", "22", "33", "12", "13", "23")]
      assert u_ij.count(None) in (0, 6)
      if u_ij.count(None) == 0:
        adps = u_ij
      else:
        assert u_ij.count(None) == 6
        b_ij = [self._wrap_loop_if_needed(cif_block, "_atom_site_anisotrop.B[%s][%s]" %(ij[0], ij[1]))
                for ij in ("11", "22", "33", "12", "13", "23")]
        assert b_ij.count(None) in (0, 6)
        if b_ij.count(None) == 0:
          adps = adptbx.b_as_u(b_ij)
        assert not (u_ij.count(None) and b_ij.count(None)) # illegal for both to be present
      if adps is not None:
        try:
          adps = [flex.double(adp) for adp in adps]
        except ValueError as e:
          raise CifBuilderError("Error interpreting ADPs: " + str(e))
        adps = flex.sym_mat3_double(*adps)
    py_adps = {}
    if anisotrop_id is not None and adps is not None:
      for an_id, adp in zip(list(anisotrop_id), list(adps)):
        py_adps[an_id] = adp
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
      assert current_auth_asym_id not in [".", "?", " "], "mmCIF file contains " + \
        "record with empty auth_asym_id, which is wrong."
      assert current_label_asym_id is not None
      if (current_auth_asym_id != last_auth_asym_id
          or current_model_id != last_model_id
          or (current_label_asym_id != last_label_asym_id and
             i_atom > 0 and is_aa_or_rna_dna(comp_id[i_atom-1]))
         ): # insert chain breaks
        chain = hierarchy.chain(id=current_auth_asym_id)
        is_first_in_chain = None
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
      if (current_residue_id != last_residue_id
          or current_ins_code != last_ins_code
          or current_auth_asym_id != last_auth_asym_id
          or current_model_id != last_model_id):
        try:
          resseq = hy36encode(width=4, value=int(current_residue_id))
        except ValueError as e:
          resseq = current_residue_id
          assert len(resseq)==4
        residue_group = hierarchy.residue_group(
          resseq=resseq,
          icode=current_ins_code)
        chain.append_residue_group(residue_group)
        if is_first_in_chain is None:
          is_first_in_chain = True
        else:
          is_first_in_chain = False
        atom_groups = OrderedDict() # reset atom_groups cache
      # atom_group(s)
      # defined by resname and altloc id
      current_altloc = alt_id[i_atom]
      if current_altloc == "." or current_altloc == "?":
        current_altloc = "" # Main chain atoms
      current_resname = comp_id[i_atom]
      if (current_altloc, current_resname) not in atom_groups:
        atom_group = hierarchy.atom_group(
          altloc=current_altloc, resname="%3s" % current_resname)
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
      if auth_segid:
        atom.set_segid(auth_segid[i_atom][:4]+(4-len(auth_segid[i_atom]))*" ")
      else:
        atom.set_segid("    ")
      if auth_break and (not is_first_in_chain) and auth_break[i_atom] == "1":
        # insert break before this residue
        residue_group.link_to_previous = False
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
        py_u_ij = py_adps.get(atom.serial.strip(), None)
        if py_u_ij is not None:
          atom.set_uij(py_u_ij)
    if len(self.hierarchy.models()) == 1:
      # for compatibility with single-model PDB files
      self.hierarchy.models()[0].id = ""

  def _wrap_loop_if_needed(self, cif_block, name):
    data = cif_block.get(name)
    if data is None:
      return data
    if isinstance(data, string_types):
      data = flex.std_string([data])
    return data


def is_aa_or_rna_dna(resname):
   from iotbx.pdb import common_residue_names_get_class
   resname = resname.strip()
   residue_class = common_residue_names_get_class(resname)
   if residue_class in ['common_amino_acid', 'modified_amino_acid',
                'common_rna_dna', 'modified_rna_dna']:
     return True
   else:
     return False

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

    if isinstance(tls_ids, string_types):
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

    if isinstance(tls_group_ids, string_types):
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

class cif_input(iotbx.pdb.pdb_input_mixin):

  def __init__(self,
               file_name=None,
               cif_object=None,
               source_info=iotbx.pdb.Please_pass_string_or_None,
               lines=None,
               raise_sorry_if_format_error=False):
    if file_name is not None:
      reader = iotbx.cif.reader(file_path=file_name)
      self.cif_model = reader.model()
    elif lines is not None:
      if not isinstance( lines, str ):
        lines = "\n".join(lines)
      reader = iotbx.cif.reader(input_string=lines)
      self.cif_model = reader.model()
    elif cif_object is not None:
      self.cif_model = cif_object
    if len(self.cif_model) == 0:
      raise Sorry("mmCIF file must contain at least one data block")
    self.cif_block = list(self.cif_model.values())[0]
    self._source_info = "file %s" %file_name
    self.hierarchy = None
    self.builder = None

  def file_type(self):
    return "mmcif"

  def construct_hierarchy(self, set_atom_i_seq=True, sort_atoms=True):
    if self.hierarchy is not None:
      return self.hierarchy
    self.builder = pdb_hierarchy_builder(self.cif_block)
    self.hierarchy = self.builder.hierarchy
    if sort_atoms:
      self.hierarchy.sort_atoms_in_place()
      if (set_atom_i_seq):
        self.hierarchy.reset_atom_i_seqs()
      self.hierarchy.atoms_reset_serial()
    return self.hierarchy

  def label_to_auth_asym_id_dictionary(self):
    auth_asym = self.cif_block.get('_atom_site.auth_asym_id')
    label_asym = self.cif_block.get('_atom_site.label_asym_id')
    assert len(label_asym) == len(auth_asym)
    return dict(zip(label_asym, auth_asym))

  def source_info(self):
    return self._source_info

  def atoms(self):
    if self.hierarchy is None:
      self.construct_hierarchy()
    return self.hierarchy.atoms()

  def atoms_with_labels(self):
    if self.hierarchy is None:
      self.construct_hierarchy()
    return self.hierarchy.atoms_with_labels()

  def model_indices(self):
    if self.hierarchy is None:
      self.construct_hierarchy()
    mi = flex.size_t([m.atoms_size() for m in self.hierarchy.models()])
    for i in range(1, len(mi)):
      mi[i] += mi[i-1]
    return mi

  def ter_indices(self):
    # for compatibility with pdb_input
    return flex.size_t()

  def crystal_symmetry(self,
                       crystal_symmetry=None,
                       weak_symmetry=False):
    self_symmetry = self.crystal_symmetry_from_cryst1()
    if (crystal_symmetry is None):
      return self_symmetry
    if (self_symmetry is None):
      return crystal_symmetry
    return self_symmetry.join_symmetry(
      other_symmetry=crystal_symmetry,
      force=not weak_symmetry)

  def crystal_symmetry_from_cryst1(self):
    if self.hierarchy is None:
      self.construct_hierarchy()
    builder_cs = self.builder.crystal_symmetry
    # check for dummy one
    if (builder_cs.unit_cell() is not None and
        cryst1_interpretation.dummy_unit_cell(
            abc = builder_cs.unit_cell().parameters()[:3],
            abg = builder_cs.unit_cell().parameters()[3:],
            sg_symbol=str(builder_cs.space_group_info()))):
      return cctbx.crystal.symmetry(
        unit_cell=None,
        space_group_info=None)
    return builder_cs


  def extract_cryst1_z_columns(self):
    return self.cif_model.values()[0].get("_cell.Z_PDB")

  def connectivity_section(self):
    # XXX Should we extract something from the CIF and return a PDB-like
    # CONECT record, or rework this whole thing?
    return ""

  def connectivity_annotation_section(self):
    return ""

  def extract_secondary_structure(self, log=None):
    from iotbx.pdb import secondary_structure
    return secondary_structure.annotation.from_cif_block(self.cif_block, log=log)

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
      try:
        return int(yyyymmdd[:4])
      except ValueError:
        pass
    return None


  def deposition_date(self, us_style=True):
    # date format: yyyy-mm-dd
    cif_block = list(self.cif_model.values())[0]
    date_orig = cif_block.get("_pdbx_database_status.recvd_initial_deposition_date")
    result = date_orig
    if not us_style:
      raise NotImplementedError
    return result
    #rev_num = cif_block.get('_database_PDB_rev.num')
    #if rev_num is not None:
    #  date_original = cif_block.get('_database_PDB_rev.date_original')
    #  if isinstance(rev_num, string_types):
    #    return date_original
    #  else:
    #    i = flex.first_index(rev_num, '1')
    #    if date_original is not None:
    #      return date_original[i]

  def get_r_rfree_sigma(self, file_name=None):
    return _cif_get_r_rfree_sigma_object(self.cif_block, file_name)

  def get_solvent_content(self):
    return _float_or_None(self.cif_block.get('_exptl_crystal.density_percent_sol'))

  def get_matthews_coeff(self):
    return _float_or_None(self.cif_block.get('_exptl_crystal.density_Matthews'))

  def get_program_name(self):
    software_name = self.cif_block.get('_software.name')
    software_classification = self.cif_block.get('_software.classification')
    if isinstance(software_classification, string_types):
      if software_classification == 'refinement':
        return software_name
    elif software_classification is not None:
      i = flex.first_index(software_classification, 'refinement')
      if (i is not None) and (i >= 0) and (software_name is not None) and (
           i < len(software_name)):
        return software_name[i]

  def resolution(self):
    result = []
    r1 = _float_or_None(self.cif_block.get('_reflns.d_resolution_high'))
    r2 = _float_or_None(self.cif_block.get('_reflns.resolution'))
    r3 = _float_or_None(self.cif_block.get('_em_3d_reconstruction.resolution'))
    r4 = _float_or_None(self.cif_block.get('_refine.ls_d_res_high'))
    for r in [r1,r2,r3,r4]:
      if(r is not None): result.append(r)
    result = list(set(result))
    if(len(result)  ==0): result = None
    elif(len(result)==1): result = result[0]
    elif(len(result)>1):
      if(r4 is None):
        result = r1
      else:
        result = min(result)
    return result

  def extract_tls_params(self, hierarchy):
    return extract_tls_from_cif_block(self.cif_block, hierarchy)

  def scale_matrix(self):
    if (not hasattr(self, "_scale_matrix")):
      fractionalization_matrix = [
        self.cif_block.get('_atom_sites.fract_transf_matrix[%s][%s]' %(i, j))
        for i,j in ('11', '12', '13', '21', '22', '23', '31','32', '33')]
      if fractionalization_matrix.count(None) == 9:
        return None
      assert fractionalization_matrix.count(None) == 0
      fractionalization_vector = [
        self.cif_block.get('_atom_sites.fract_transf_vector[%s]' %i) for i in '123']
      assert fractionalization_matrix.count(None) == 0
      self._scale_matrix = [[float(i) for i in fractionalization_matrix],
                            [float(i) for i in fractionalization_vector]]
    return self._scale_matrix

  def model_ids(self):
    return flex.std_string([model.id for model in self.hierarchy.models()])

  def extract_f_model_core_constants(self):
    return extract_f_model_core_constants(self.cif_block)

  def extract_wavelength(self, first_only=True):
    wavelengths = self.cif_block.get('_diffrn_source.pdbx_wavelength_list')
    wavelength = _float_or_None(self.cif_block.get(
        '_diffrn_source.pdbx_wavelength'))
    if (first_only):
      if (wavelength is not None):
        return wavelength
      elif (wavelengths is not None):
        wavelengths = [ float(f.strip()) for f in wavelengths.split(",") ]
        return wavelengths[0]
    elif (wavelengths is not None):
        return [ float(f) for f in wavelengths.split(",") ]
    return None

  def get_experiment_type(self):
    exptl_method = self.cif_block.get('_exptl.method')
    lines = []
    if exptl_method is not None:
      if(isinstance(exptl_method,flex.std_string)):
        lines = list(exptl_method)
      else:
        lines = [exptl_method]
    return experiment_type(lines)

  def process_BIOMT_records(self):
    import iotbx.mtrix_biomt
    return iotbx.mtrix_biomt.process_BIOMT_records_cif(
      cif_block = self.cif_block)

  def process_MTRIX_records(self):
    """
    Read MTRIX records from a pdb file

    Arguments:
    ----------

    Returns:
    --------
    result : object containing information on all NCS operations
             iotbx.mtrix_biomt.container
    """
    import iotbx.mtrix_biomt
    return iotbx.mtrix_biomt.process_MTRIX_records_cif(
      cif_block = self.cif_block)

  def get_restraints_used(self):
    return {'CDL' : self.used_cdl_restraints(),
            'omega' : self.used_omega_restraints(),
            'Amber' : self.used_amber_restraints(),
           }

  def _used_what_restraints(self, what):
    rc = False
    for cif_key, cif_block in self.cif_model.items():
      target = cif_block.get("_refine.pdbx_stereochemistry_target_values")
      if (target is not None) and (what in target):
        rc = True
        break
    return rc

  def used_cdl_restraints(self):
    return self._used_what_restraints('CDL')

  def used_omega_cdl_restraints(self):
    return self._used_what_restraints('omega-cdl')

  def used_amber_restraints(self):
    return self._used_what_restraints('Amber')

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
    print(self.formatted_string(), file=log)


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

  twin_fraction = _float_or_None(cif_block.get('_pdbx_reflns_twin.fraction'))
  twin_law = cif_block.get('_pdbx_reflns_twin.operator')
  # TODO: extract these from the CIF?
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
