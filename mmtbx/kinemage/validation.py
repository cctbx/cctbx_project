from __future__ import absolute_import, division, print_function
import os
from iotbx import pdb
from mmtbx import monomer_library
from mmtbx.validation.rotalyze import rotalyze
from mmtbx.validation.ramalyze import ramalyze
from mmtbx.validation.cbetadev import cbetadev
from mmtbx.validation import rna_validate
from mmtbx.validation import omegalyze
from mmtbx.kinemage import kin_vec
from iotbx.pdb import common_residue_names_get_class
from libtbx import Auto
from scitbx import matrix
from cctbx import geometry_restraints
from mmtbx.monomer_library import pdb_interpretation
import iotbx.phil
from libtbx.utils import Sorry, Usage
from libtbx.str_utils import show_string
import libtbx.load_env

def get_master_phil():
  return iotbx.phil.parse(
    input_string="""
      kinemage {
        pdb = None
         .type = path
         .optional = True
         .help = '''Enter a PDB file name'''
        cif = None
         .type = path
         .optional = True
         .multiple = True
         .help = '''Enter a CIF file for ligand parameters'''
        out_file = None
         .type = path
         .optional = True
         .help = '''Enter a .kin output name'''
        keep_hydrogens = False
        .type = bool
        .help = '''Keep hydrogens in input file'''
        pdb_interpretation
          .short_caption = Model interpretation
        {
          include scope mmtbx.monomer_library.pdb_interpretation.master_params_str
        }
  }
    """,
    process_includes=True,
  )

def build_name_hash(pdb_hierarchy):
  i_seq_name_hash = dict()
  for atom in pdb_hierarchy.atoms():
    labels = atom.fetch_labels()
    i_seq_name_hash[atom.i_seq] = {
      'name':     labels.name,
      'altloc':   labels.altloc,
      'resname':  labels.resname,
      'chain_id': labels.chain_id,
      'resseq':   labels.resseq,
      'icode':    labels.icode,
    }
  return i_seq_name_hash


def angle_outlier_as_kinemage(self):
  """
  Represent the angle outlier as a kinemage 'fan'.
  """
  from mmtbx.kinemage import kin_vec
  from scitbx import matrix
  kin_text = ""
  atom0 = self.atoms_info[0]
  angle_key = self.kinemage_key()
  num_sigmas = - self.delta / self.sigma
  angle_key_full = "%s %.3f sigma" % (angle_key, num_sigmas)
  if (num_sigmas < 0):
    color = "blue"
  else:
    color = "red"
  sites = self.sites_cart()
  delta = self.delta
  a = matrix.col(sites[0])
  b = matrix.col(sites[1])
  c = matrix.col(sites[2])
  normal = (a-b).cross(c-b).normalize()
  r = normal.axis_and_angle_as_r3_rotation_matrix(angle=delta, deg=True)
  new_c = tuple( (r*(c-b)) +b)
  kin_text += "@vectorlist {%s} color= %s width= 4 master= {angle dev}\n" \
               % (angle_key_full, color)
  kin_text += kin_vec(angle_key_full,sites[0],angle_key_full,sites[1])
  kin_text += kin_vec(angle_key_full,sites[1],angle_key_full,new_c)
  r = normal.axis_and_angle_as_r3_rotation_matrix(angle=(delta*.75), deg=True)
  new_c = tuple( (r*(c-b)) +b)
  kin_text += "@vectorlist {%s} color= %s width= 3 master= {angle dev}\n" \
               % (angle_key_full, color)
  kin_text += kin_vec(angle_key_full,sites[1],angle_key,new_c)
  r = normal.axis_and_angle_as_r3_rotation_matrix(angle=(delta*.5), deg=True)
  new_c = tuple( (r*(c-b)) +b)
  kin_text += "@vectorlist {%s} color= %s width= 2 master= {angle dev}\n" \
               % (angle_key_full, color)
  kin_text += kin_vec(angle_key_full,sites[1],angle_key,new_c)
  r = normal.axis_and_angle_as_r3_rotation_matrix(angle=(delta*.25), deg=True)
  new_c = tuple( (r*(c-b)) +b)
  kin_text += "@vectorlist {%s} color= %s width= 1 master= {angle dev}\n" \
              % (angle_key_full, color)
  kin_text += kin_vec(angle_key_full,sites[1],angle_key,new_c)
  return kin_text

def chiral_outlier_as_kinemage(self):
  """
  Represent a chiral volume outlier in kinemage
  """
  from mmtbx.kinemage import kin_vec
  outlier_type = self.outlier_type()
  atoms = self.atoms_info
  chiral_center = atoms[0]
  chiral_coords = matrix.rec((chiral_center.xyz[0], chiral_center.xyz[1], chiral_center.xyz[2]),(3,1))
  legs = [None] #None fills the 0 index, this indexes the same as atoms_info
  for atom in atoms[1:]:
    legs.append(matrix.rec((atom.xyz[0], atom.xyz[1], atom.xyz[2]),(3,1)))
  #legs need to be shortened for visual/selection clarity
  leg_vs = [None]
  for leg in legs[1:]:
    leg_vs.append(leg - chiral_coords)
  leg_ends = [None,0,0,0]
  shorten_by = 0.3
  leg_ends[1] = legs[1] - leg_vs[1].normalize()*shorten_by
  leg_ends[2] = legs[2] - leg_vs[2].normalize()*shorten_by
  leg_ends[3] = legs[3] - leg_vs[3].normalize()*shorten_by
  kin_text = ""
  #tetrahedron is drawn for all markups
  kin_text += "@vectorlist {%s %s} color= yellow width= 4 master={chiral dev}\n" %(chiral_center.id_str(), "lines")
  #draw legs from chiral center
  kin_text += "{%s %s: %.3f sigma} P %.3f %.3f %.3f\n" % (chiral_center.id_str(),outlier_type,self.score,chiral_coords[0],chiral_coords[1],chiral_coords[2])
  kin_text += "{%s %s: %.3f sigma} %.3f %.3f %.3f\n" % (atoms[1].id_str(),outlier_type,self.score,leg_ends[1][0],leg_ends[1][1],leg_ends[1][2])
  kin_text += "{%s %s: %.3f sigma} %.3f %.3f %.3f\n" % (atoms[2].id_str(),outlier_type,self.score,leg_ends[2][0],leg_ends[2][1],leg_ends[2][2])
  kin_text += "{%s %s: %.3f sigma} %.3f %.3f %.3f\n" % (atoms[3].id_str(),outlier_type,self.score,leg_ends[3][0],leg_ends[3][1],leg_ends[3][2])
  kin_text += "{%s %s: %.3f sigma} %.3f %.3f %.3f\n" % (atoms[1].id_str(),outlier_type,self.score,leg_ends[1][0],leg_ends[1][1],leg_ends[1][2])
  kin_text += "{%s %s: %.3f sigma} P %.3f %.3f %.3f\n" % (atoms[2].id_str(),outlier_type,self.score,leg_ends[2][0],leg_ends[2][1],leg_ends[2][2])
  kin_text += "{%s %s: %.3f sigma} %.3f %.3f %.3f\n" % (chiral_center.id_str(),outlier_type,self.score,chiral_coords[0],chiral_coords[1],chiral_coords[2])
  kin_text += "{%s %s: %.3f sigma} %.3f %.3f %.3f\n" % (atoms[3].id_str(),outlier_type,self.score,leg_ends[3][0],leg_ends[3][1],leg_ends[3][2])

  if outlier_type == "Chiral handedness swap":
    #Also add an arrow to suggest mirror operation
    #create normal to plane of 3 legs, scale to suitable length
    #m = rec((1, 2, 3, 4), (2, 2))
    v1 = matrix.rec((atoms[1].xyz[0]-atoms[2].xyz[0], atoms[1].xyz[1]-atoms[2].xyz[1], atoms[1].xyz[2]-atoms[2].xyz[2]), (3,1))
    v2 = matrix.rec((atoms[1].xyz[0]-atoms[3].xyz[0], atoms[1].xyz[1]-atoms[3].xyz[1], atoms[1].xyz[2]-atoms[3].xyz[2]), (3,1))
    norm = v1.cross(v2)
    norm = norm.normalize()
    p = matrix.rec((chiral_center.xyz[0]+norm[0], chiral_center.xyz[1]+norm[1], chiral_center.xyz[2]+norm[2]),(3,1))
    arrow_text = "%s %s: %.3f sigma" % (chiral_center.id_str(), outlier_type, self.score)
    kin_text += kin_vec(chiral_center.id_str(), chiral_center.xyz, arrow_text, p)
    #Add arrowhead
    #Move back lone that some line and out along the atoms[1]-atoms[2] line
    arrow_width = 0.125 * v1.normalize()
    arrow_end_1 = p - 0.125*norm + arrow_width
    arrow_end_2 = p - 0.125*norm - arrow_width
    #draw second arrow rotated 90 degrees, so arrowheaad is visible from more orientations
    arrow_end_3 = matrix.rotate_point_around_axis(
        axis_point_1 = chiral_coords,
        axis_point_2 = p,
        point        = arrow_end_1,
        angle        = 90,
        deg          = True)
    arrow_end_4 = matrix.rotate_point_around_axis(
        axis_point_1 = chiral_coords,
        axis_point_2 = p,
        point        = arrow_end_2,
        angle        = 90,
        deg          = True)
    kin_text += kin_vec(arrow_text, p, outlier_type, arrow_end_1)
    kin_text += kin_vec(arrow_text, p, outlier_type, arrow_end_2)
    kin_text += kin_vec(arrow_text, p, outlier_type, arrow_end_3)
    kin_text += kin_vec(arrow_text, p, outlier_type, arrow_end_4)

  #draw balls onto atoms of greatest interest
  kin_text += "@balllist {%s %s} color= yellow radius= 0.15 master={chiral dev}\n" %(chiral_center.id_str(), "balls")
  if outlier_type != 'Pseudochiral naming error':
    #"Chiral handedness swap" or "Tetrahedral geometry outlier"
    #error is located at chiral center, draw ball at chiral center of interest
    kin_text += "{%s %s: %.3f sigma} %.3f %.3f %.3f\n" % (chiral_center.id_str(),outlier_type,self.score,atoms[0].xyz[0],atoms[0].xyz[1],atoms[0].xyz[2])
  else:
    #Pseudochiral naming error
    #error is probably in naming of the non-center atoms
    #draw balls on legs, add atomnames as labels to draw attention to naming
    i = 1
    while i <= 3:
      kin_text += "{%s %s: %.3f sigma} %.3f %.3f %.3f\n" % (atoms[i].id_str(),outlier_type,self.score,leg_ends[i][0],leg_ends[i][1],leg_ends[i][2])
      i+=1
    kin_text += "@labellist {%s %s} color= white master={chiral dev}\n" % (chiral_center.id_str(), "labels")
    #needs to be a different color than balls for readability during overlap
    #atomnames go on atoms rather than balls to reduce overlap and increase clarity
    kin_text += "{  %s} %.3f %.3f %.3f\n" % (chiral_center.name.strip(),chiral_center.xyz[0],chiral_center.xyz[1],chiral_center.xyz[2])
    i = 1
    while i <= 3:
      kin_text += "{  %s} %.3f %.3f %.3f\n" % (atoms[i].name.strip(),atoms[i].xyz[0],atoms[i].xyz[1],atoms[i].xyz[2])
      #leading spaces reduce overlap with the ball already drawn on labeled atom
      i+=1
  return kin_text


def bond_outlier_as_kinemage(self):
  """
  Represent the bond outlier as a kinemage spring.
  """
  from mmtbx.kinemage import kin_vec
  from scitbx import matrix
  kin_text = ""
  bond_key = self.kinemage_key()
  num_sigmas = - self.delta / self.sigma
  if (num_sigmas < 0):
    color = "blue"
  else:
    color = "red"
  sites = self.sites_cart()
  a = matrix.col(sites[0])
  b = matrix.col(sites[1])
  c = matrix.col( (1,0,0) )
  normal = ((a-b).cross(c-b).normalize())*0.2
  current = a+normal
  new = tuple(current)
  kin_text += \
    "@vectorlist {%s %.3f sigma} color= %s width= 3 master= {length dev}\n"\
              % (bond_key, num_sigmas, color)
  bond_key_long = "%s %.3f sigma" % (bond_key, num_sigmas)
  kin_text += kin_vec(bond_key_long,sites[0],bond_key_long,new)
  angle = 36
  dev = num_sigmas
  if dev > 10.0:
    dev = 10.0
  if dev < -10.0:
    dev = -10.0
  if dev <= 0.0:
    angle += 1.5*abs(dev)
  elif dev > 0.0:
    angle -= 1.5*dev
  i = 0
  n = 60
  axis = b-a
  step = axis*(1.0/n)
  r = axis.axis_and_angle_as_r3_rotation_matrix(angle=angle, deg=True)
  while i < n:
    next = (r*(current-b) +b)
    next = next + step
    kin_text += kin_vec(bond_key_long,tuple(current),bond_key_long,
      tuple(next))
    current = next
    i += 1
  kin_text += kin_vec(bond_key_long,tuple(current),bond_key_long,sites[1])
  return kin_text


def make_probe_dots(hierarchy, keep_hydrogens=False):
  """Generate probe dot kinemage output using probe2 Python API.

  Uses mmtbx.reduce (reduce2) for hydrogen placement and mmtbx.programs.probe2
  for contact analysis, producing kinemage-format dot output.
  """
  try:
    from mmtbx.hydrogens import reduce_hydrogen
    from mmtbx.reduce import Optimizers
    from mmtbx.programs import probe2
    import mmtbx.model
    from libtbx.utils import null_out
    import tempfile
  except ImportError:
    return ""

  probe_return = ""
  for i_mod, m in enumerate(hierarchy.models()):
    r = pdb.hierarchy.root()
    mdc = m.detached_copy()
    r.append_model(mdc)

    # Build a model manager for this sub-model
    model_manager = mmtbx.model.manager(
      model_input=None,
      pdb_hierarchy=r,
      stop_for_unknowns=False,
      log=null_out())
    model_manager.add_crystal_symmetry_if_necessary()

    # Add hydrogens if needed
    if not keep_hydrogens:
      try:
        reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
          model=model_manager,
          use_neutron_distances=False,
          n_terminal_charge="residue_one",
          exclude_water=True,
          stop_for_unknowns=False,
          keep_existing_H=False)
        reduce_add_h_obj.run()
        model_manager = reduce_add_h_obj.get_model()

        # Build probe parameters for optimizer using probe2 defaults
        # Optimizer expects the probe sub-params, not the top-level params
        import iotbx.phil
        probe_phil = iotbx.phil.parse(probe2.master_phil_str, process_includes=True)
        probe_params = probe_phil.extract()
        opt = Optimizers.Optimizer(probe_params.probe, False, model_manager,
          modelIndex=None, fillAtomDump=False)
      except Exception:
        # If hydrogen addition fails, continue with existing atoms
        pass

    # Rebuild model manager after hydrogen changes
    model_manager = mmtbx.model.manager(
      model_input=None,
      pdb_hierarchy=model_manager.get_hierarchy(),
      stop_for_unknowns=False,
      crystal_symmetry=model_manager.crystal_symmetry(),
      log=null_out())

    # Run probe2 in kinemage output mode
    try:
      import iotbx.cli_parser

      tempName = tempfile.mktemp()
      parser = iotbx.cli_parser.CCTBXParser(
        program_class=probe2.Program, logger=null_out())
      args = [
        "approach=self",
        "output.format=kinemage",
        "output.filename='%s'" % tempName,
        "output.separate_worse_clashes=True",
        "output.report_vdws=False",
        "output.write_files=False",
        "count_dots=False",
        "ignore_lack_of_explicit_hydrogens=True",
      ]
      parser.parse_args(args)
      dm = parser.data_manager
      p2 = probe2.Program(dm, parser.working_phil.extract(),
                          master_phil=parser.master_phil, logger=null_out())
      p2.overrideModel(model_manager)
      dots, output = p2.run()
      probe_return += output
      if os.path.exists(tempName):
        os.unlink(tempName)
    except Exception:
      # If probe2 fails, return what we have so far
      pass
  return probe_return

def cbeta_dev(outliers, chain_id=None):
  cbeta_out = "@subgroup {CB dev} dominant\n"
  cbeta_out += "@balllist {CB dev Ball} color= magenta master= {Cbeta dev}\n"
  for outlier in outliers :
    if (chain_id is None) or (chain_id == outlier.chain_id):
      cbeta_out += outlier.as_kinemage() + "\n"
  if len(cbeta_out.splitlines()) == 2:
    cbeta_out = ""
  return cbeta_out

def midpoint(p1, p2):
  mid = [0.0, 0.0, 0.0]
  mid[0] = (p1[0]+p2[0])/2
  mid[1] = (p1[1]+p2[1])/2
  mid[2] = (p1[2]+p2[2])/2
  return mid

def get_chain_color(index):
  chain_colors = ['white',
                  'yellowtint',
                  'peachtint',
                  'pinktint',
                  'lilactint',
                  'bluetint',
                  'greentint']
  return chain_colors[int(index) % len(chain_colors)]

def get_ions(ion_list):
  ion_txt = "@spherelist {het M} color= gray  radius= 0.5 nobutton master= {hets}\n"
  ion_txt += ion_list
  return ion_txt

def _get_prev_connection(prev_key_hash, prev_xyz_hash, altloc):
  """Look up a previous residue's key/xyz for backbone connection, falling
  back to the blank altloc if the current altloc is not found."""
  prev_key = prev_key_hash.get(altloc)
  prev_xyz = prev_xyz_hash.get(altloc)
  if prev_key is None:
    prev_key = prev_key_hash.get(' ')
    prev_xyz = prev_xyz_hash.get(' ')
  return prev_key, prev_xyz

def _track_amino_acid_atom(atom, key, altloc, residue_group, prev_resid,
                           cur_C_xyz, cur_C_key, cur_CA_xyz, cur_CA_key,
                           prev_C_key, prev_C_xyz, prev_CA_key, prev_CA_xyz,
                           mc_parts, ca_parts):
  """Track backbone atoms (C, CA, N) for amino acids and add inter-residue
  connections to mc_parts and ca_parts lists."""
  if atom.name == ' C  ':
    cur_C_xyz[altloc] = atom.xyz
    cur_C_key[altloc] = key
  if atom.name == ' CA ':
    cur_CA_xyz[altloc] = atom.xyz
    cur_CA_key[altloc] = key
    if len(prev_CA_key) > 0 and len(prev_CA_xyz) > 0:
      if prev_resid is not None and \
         int(residue_group.resseq_as_int()) - int(prev_resid[0:4]) == 1:
        prev_key, prev_xyz = _get_prev_connection(prev_CA_key, prev_CA_xyz, altloc)
        if prev_key is not None:
          ca_parts.append(kin_vec(prev_key, prev_xyz, key, atom.xyz))
  if atom.name == ' N  ':
    if len(prev_C_key) > 0 and len(prev_C_xyz) > 0:
      if prev_resid is not None and \
         int(residue_group.resseq_as_int()) - int(prev_resid[0:4]) == 1:
        prev_key, prev_xyz = _get_prev_connection(prev_C_key, prev_C_xyz, altloc)
        if prev_key is not None:
          mc_parts.append(kin_vec(prev_key, prev_xyz, key, atom.xyz))

def _track_rna_dna_atom(atom, key, altloc, residue_group, prev_resid,
                        cur_O3_xyz, cur_O3_key,
                        prev_O3_key, prev_O3_xyz,
                        p_hash_key, p_hash_xyz,
                        c1_hash_key, c1_hash_xyz,
                        c4_hash_key, c4_hash_xyz,
                        mc_parts):
  """Track backbone atoms for RNA/DNA and add O3'-P connections."""
  if atom.name == " O3'":
    cur_O3_xyz[altloc] = atom.xyz
    cur_O3_key[altloc] = key
  elif atom.name == ' P  ':
    if len(prev_O3_key) > 0 and len(prev_O3_xyz) > 0:
      if prev_resid is not None and \
         int(residue_group.resseq_as_int()) - int(prev_resid[0:4]) == 1:
        prev_key, prev_xyz = _get_prev_connection(prev_O3_key, prev_O3_xyz, altloc)
        if prev_key is not None:
          mc_parts.append(kin_vec(prev_key, prev_xyz, key, atom.xyz))
    resseq = residue_group.resseq_as_int()
    p_hash_key[resseq] = key
    p_hash_xyz[resseq] = atom.xyz
  elif atom.name == " C1'":
    c1_hash_key[residue_group.resseq_as_int()] = key
    c1_hash_xyz[residue_group.resseq_as_int()] = atom.xyz
  elif atom.name == " C4'":
    c4_hash_key[residue_group.resseq_as_int()] = key
    c4_hash_xyz[residue_group.resseq_as_int()] = atom.xyz

def _draw_rna_virtual_backbone(residue_group, p_hash_key, p_hash_xyz,
                               c1_hash_key, c1_hash_xyz,
                               c4_hash_key, c4_hash_xyz):
  """Generate virtual backbone vectors for RNA/DNA residues (C4'->P->C4'->C1')."""
  vbb = ""
  resseq = residue_group.resseq_as_int()
  # C4'(prev) -> P(cur)
  if (resseq - 1) in c4_hash_key and resseq in p_hash_key:
    vbb += kin_vec(c4_hash_key[resseq-1], c4_hash_xyz[resseq-1],
                   p_hash_key[resseq], p_hash_xyz[resseq])
  # P(cur) -> C4'(cur)
  if resseq in p_hash_key and resseq in c4_hash_key:
    vbb += kin_vec(p_hash_key[resseq], p_hash_xyz[resseq],
                   c4_hash_key[resseq], c4_hash_xyz[resseq])
  # C4'(cur) -> C1'(cur)
  if resseq in c4_hash_key and resseq in c1_hash_key:
    vbb += kin_vec(c4_hash_key[resseq], c4_hash_xyz[resseq],
                   c1_hash_key[resseq], c1_hash_xyz[resseq])
  return vbb

def _draw_residue_bonds(residue, bond_hash, i_seq_name_hash, key_hash,
                        xyz_hash, het_hash, iseq_altloc, altloc,
                        mc_atoms, show_hydrogen, drawn_bonds):
  """Draw bonds for a residue, sorting them into the appropriate kin lists.

  Returns a dict with keys: mc, sc, mc_h, sc_h, het, het_h for the
  kinemage vectors generated."""
  result = {'mc': '', 'sc': '', 'mc_h': '', 'sc_h': '', 'het': '', 'het_h': ''}
  res_class = common_residue_names_get_class(residue.resname)
  for atom in residue.atoms():
    cur_bonds = bond_hash.get(atom.i_seq)
    if cur_bonds is None:
      continue
    for bond in cur_bonds:
      info_1 = i_seq_name_hash.get(atom.i_seq)
      atom_1 = info_1['name'].strip() if info_1 is not None else None
      info_2 = i_seq_name_hash.get(bond)
      atom_2 = info_2['name'].strip() if info_2 is not None else None
      if atom_1 is None or atom_2 is None:
        continue
      # handle altlocs
      if key_hash.get(atom_1) is None or key_hash.get(atom_2) is None:
        continue
      drawn_key = key_hash[atom_1]+key_hash[atom_2]
      if drawn_key in drawn_bonds:
        continue
      altloc_2 = iseq_altloc.get(bond)
      if altloc_2 != altloc and altloc_2 != '':
        continue
      is_hydrogen = (atom_1.startswith('H') or atom_2.startswith('H') or
                     atom_1.startswith('D') or atom_2.startswith('D'))
      if res_class in ("common_amino_acid", "common_rna_dna"):
        if atom_1 in mc_atoms and atom_2 in mc_atoms:
          # Skip inter-residue bonds handled separately (C-N, O3'-P)
          if not ((atom_1 == "C" and atom_2 == "N") or
                  (atom_1 == "O3'" and atom_2 == "P")):
            result['mc'] += kin_vec(key_hash[atom_1], xyz_hash[atom_1],
                                    key_hash[atom_2], xyz_hash[atom_2])
        elif is_hydrogen:
          if show_hydrogen:
            if atom_1 in mc_atoms or atom_2 in mc_atoms:
              result['mc_h'] += kin_vec(key_hash[atom_1], xyz_hash[atom_1],
                                        key_hash[atom_2], xyz_hash[atom_2])
            else:
              result['sc_h'] += kin_vec(key_hash[atom_1], xyz_hash[atom_1],
                                        key_hash[atom_2], xyz_hash[atom_2])
        else:
          result['sc'] += kin_vec(key_hash[atom_1], xyz_hash[atom_1],
                                  key_hash[atom_2], xyz_hash[atom_2])
        drawn_bonds.append(drawn_key)
      else:
        # Het/ligand bonds (covers 'other', 'common_small_molecule',
        # 'common_saccharide', and any other non-polymer classes)
        if is_hydrogen:
          if show_hydrogen and atom_1 in het_hash and atom_2 in het_hash:
            result['het_h'] += kin_vec(het_hash[atom_1][0], het_hash[atom_1][1],
                                       het_hash[atom_2][0], het_hash[atom_2][1])
        else:
          if atom_1 in het_hash and atom_2 in het_hash:
            result['het'] += kin_vec(het_hash[atom_1][0], het_hash[atom_1][1],
                                     het_hash[atom_2][0], het_hash[atom_2][1])
  return result

def get_kin_lots(chain, bond_hash, i_seq_name_hash, pdbID=None, index=0,
                 show_hydrogen=True, ss_bonds=None, sites_cart=None):
  mc_atoms = ["N", "CA", "C", "O", "OXT",
              "P", "OP1", "OP2", "OP3", "O5'", "C5'", "C4'", "O4'", "C1'",
              "C3'", "O3'", "C2'", "O2'"]
  kin_out = ""
  color = get_chain_color(index)
  mc_veclist = "@vectorlist {mc} color= %s  master= {mainchain}\n" % color
  sc_veclist = "@vectorlist {sc} color= cyan  master= {sidechain}\n"
  ca_trace = "@vectorlist {Calphas} color= %s master= {Calphas}\n" % color
  virtual_bb = "@vectorlist {Virtual BB} color= %s  off   master= {Virtual BB}\n" % color
  water_list = "@balllist {water O} color= peachtint  radius= 0.15  master= {water}\n"
  hets = "@vectorlist {het} color= pink  master= {hets}\n"
  het_h = "@vectorlist {ht H} color= gray  nobutton master= {hets} master= {H's}\n"
  mc_h_veclist = ""
  sc_h_veclist = ""
  if show_hydrogen:
    mc_h_veclist = \
      "@vectorlist {mc H} color= gray nobutton master= {mainchain} master= {H's}\n"
    sc_h_veclist = \
      "@vectorlist {sc H} color= gray nobutton master= {sidechain} master= {H's}\n"
  ion_list = ""
  prev_resid = None
  prev_C_xyz = {}
  prev_C_key = {}
  prev_CA_xyz = {}
  prev_CA_key = {}
  prev_O3_xyz = {}
  prev_O3_key = {}
  p_hash_key = {}
  p_hash_xyz = {}
  c1_hash_key = {}
  c1_hash_xyz = {}
  c4_hash_key = {}
  c4_hash_xyz = {}
  drawn_bonds = []

  for residue_group in chain.residue_groups():
    altloc_hash = {}
    iseq_altloc = {}
    cur_C_xyz = {}
    cur_C_key = {}
    cur_CA_xyz = {}
    cur_CA_key = {}
    cur_O3_xyz = {}
    cur_O3_key = {}
    for atom_group in residue_group.atom_groups():
      ag_altloc = atom_group.altloc
      for atom in atom_group.atoms():
        if altloc_hash.get(atom.name.strip()) is None:
          altloc_hash[atom.name.strip()] = []
        altloc_hash[atom.name.strip()].append(ag_altloc)
        iseq_altloc[atom.i_seq] = ag_altloc
    cur_resid = residue_group.resid()
    for conformer in residue_group.conformers():
      for residue in conformer.residues():
        cur_resid = residue.resid()
        key_hash = {}
        xyz_hash = {}
        het_hash = {}
        altloc = conformer.altloc
        if altloc == '':
          altloc = ' '
        for atom in residue.atoms():
          cur_altlocs = altloc_hash.get(atom.name.strip())
          if cur_altlocs == ['']:
            cur_altloc = ' '
          elif altloc in cur_altlocs:
            cur_altloc = altloc
          else:
            cur_altloc = ' '
          key = "%s%s%s %s%s  B%.2f %s" % (
                atom.name.lower(),
                cur_altloc.lower(),
                residue.resname.lower(),
                chain.id,
                residue_group.resid(),
                atom.b,
                pdbID)
          key_hash[atom.name.strip()] = key
          xyz_hash[atom.name.strip()] = atom.xyz
          res_class = common_residue_names_get_class(residue.resname)
          if res_class == "common_amino_acid":
            mc_parts = []
            ca_parts = []
            _track_amino_acid_atom(
              atom, key, altloc, residue_group, prev_resid,
              cur_C_xyz, cur_C_key, cur_CA_xyz, cur_CA_key,
              prev_C_key, prev_C_xyz, prev_CA_key, prev_CA_xyz,
              mc_parts, ca_parts)
            for part in mc_parts:
              mc_veclist += part
            for part in ca_parts:
              ca_trace += part
          elif res_class == "common_rna_dna":
            mc_parts = []
            _track_rna_dna_atom(
              atom, key, altloc, residue_group, prev_resid,
              cur_O3_xyz, cur_O3_key,
              prev_O3_key, prev_O3_xyz,
              p_hash_key, p_hash_xyz,
              c1_hash_key, c1_hash_xyz,
              c4_hash_key, c4_hash_xyz,
              mc_parts)
            for part in mc_parts:
              mc_veclist += part
          elif res_class == "common_element":
            ion_list += "{%s} %.3f %.3f %.3f\n" % (
              key, atom.xyz[0], atom.xyz[1], atom.xyz[2])
          elif res_class == "common_water":
            if atom.name == ' O  ':
              water_list += "{%s} P %.3f %.3f %.3f\n" % (
                key, atom.xyz[0], atom.xyz[1], atom.xyz[2])
          elif len(residue.atoms()) == 1:
            # Single-atom non-polymer residue (e.g., lone S from SO4,
            # unknown ligand with one atom) — draw as a sphere
            ion_list += "{%s} %.3f %.3f %.3f\n" % (
              key, atom.xyz[0], atom.xyz[1], atom.xyz[2])
          else:
            het_hash[atom.name.strip()] = [key, atom.xyz]

        # Virtual backbone for RNA/DNA
        if common_residue_names_get_class(residue.resname) == "common_rna_dna":
          virtual_bb += _draw_rna_virtual_backbone(
            residue_group, p_hash_key, p_hash_xyz,
            c1_hash_key, c1_hash_xyz, c4_hash_key, c4_hash_xyz)

        # Draw bonds
        bond_result = _draw_residue_bonds(
          residue, bond_hash, i_seq_name_hash, key_hash,
          xyz_hash, het_hash, iseq_altloc, altloc,
          mc_atoms, show_hydrogen, drawn_bonds)
        mc_veclist += bond_result['mc']
        sc_veclist += bond_result['sc']
        mc_h_veclist += bond_result['mc_h']
        sc_h_veclist += bond_result['sc_h']
        hets += bond_result['het']
        het_h += bond_result['het_h']

    prev_CA_xyz = cur_CA_xyz
    prev_CA_key = cur_CA_key
    prev_C_xyz = cur_C_xyz
    prev_C_key = cur_C_key
    prev_resid = cur_resid
    prev_O3_key = cur_O3_key
    prev_O3_xyz = cur_O3_xyz

  ion_kin = None
  if len(ion_list) > 1:
    ion_kin = get_ions(ion_list)

  # Only include non-empty lists
  if len(mc_veclist.splitlines()) > 1:
    kin_out += mc_veclist
  if len(mc_h_veclist.splitlines()) > 1:
    kin_out += mc_h_veclist
  if len(ca_trace.splitlines()) > 1:
    kin_out += ca_trace
  if len(sc_veclist.splitlines()) > 1:
    kin_out += sc_veclist
  if len(sc_h_veclist.splitlines()) > 1:
    kin_out += sc_h_veclist
  # Draw disulfide bonds for this chain
  if ss_bonds and sites_cart is not None:
    ss_veclist = "@vectorlist {SS} color= yellow  master= {sidechain}\n"
    # Collect i_seqs involved in SS bonds for this chain
    ss_i_seqs = set()
    for i_seq_0, i_seq_1 in ss_bonds:
      info_0 = i_seq_name_hash.get(i_seq_0)
      if info_0 is not None and info_0['chain_id'] == chain.id:
        ss_i_seqs.add(i_seq_0)
        ss_i_seqs.add(i_seq_1)
    # Build B-factor lookup for just the SS atoms we need
    b_hash = {}
    if ss_i_seqs:
      for atom in chain.parent().parent().atoms():
        if atom.i_seq in ss_i_seqs:
          b_hash[atom.i_seq] = atom.b
    def _ss_key(info, i_seq):
      return "%s%s%s %s%s%s  B%.2f %s" % (
        info['name'].lower(), info['altloc'].lower(), info['resname'].lower(),
        info['chain_id'], info['resseq'], info['icode'],
        b_hash.get(i_seq, 0.0), pdbID)
    for i_seq_0, i_seq_1 in ss_bonds:
      info_0 = i_seq_name_hash.get(i_seq_0)
      info_1 = i_seq_name_hash.get(i_seq_1)
      if info_0 is None or info_1 is None:
        continue
      # Only draw when the first atom belongs to this chain (avoid duplicates)
      if info_0['chain_id'] != chain.id:
        continue
      key_0 = _ss_key(info_0, i_seq_0)
      key_1 = _ss_key(info_1, i_seq_1)
      ss_veclist += kin_vec(key_0, sites_cart[i_seq_0],
                            key_1, sites_cart[i_seq_1])
    if len(ss_veclist.splitlines()) > 1:
      kin_out += ss_veclist
  if len(water_list.splitlines()) > 1:
    kin_out += water_list
  if len(virtual_bb.splitlines()) > 1:
    kin_out += virtual_bb
  if len(hets.splitlines()) > 1:
    kin_out += hets
  if ion_kin is not None:
    kin_out += ion_kin
  if len(het_h.splitlines()) > 1:
    kin_out += het_h
  return kin_out

def get_default_header():
  header = """
@kinemage 1
@onewidth
@1viewid {Overview}
@master {mainchain} indent
@master {Calphas} indent
@master {Virtual BB} indent
@master {sidechain} indent
@master {H's} indent
@pointmaster 'H' {H's}
@master {water} indent
@master {protein ribbon} indent
@master {NA ribbon} indent
"""
  return header

def get_footer():
  footer = """
@master {mainchain} off
@master {sidechain} off
@master {H's} off
@master {water} off
@master {Rota outliers} on
@master {Rama outliers} on
@master {Calphas} on
@master {Virtual BB} on
@master {vdw contact} off
@master {small overlap} off
@master {H-bonds} off
@master {length dev} on
@master {angle dev} on
@master {Cbeta dev} on
@master {base-P perp} on
@master {hets} on
@master {protein ribbon} off
@master {NA ribbon} off
"""
  return footer

def get_altid_controls(hierarchy):
  #print hierarchy.get_conformer_indices()
  #STOP()
  altids = []
  altid_controls = ""
  txt = "@pointmaster 'a' {alta} on"
  for model in hierarchy.models():
    for chain in model.chains():
      for conformer in chain.conformers():
        altloc = conformer.altloc
        if altloc not in altids and altloc != '':
          altids.append(altloc)
  for i, altid in enumerate(altids):
    altid_controls += "@pointmaster '%s' {alt%s} " % (altid.lower(),
                                                      altid.lower())
    if i == 0:
      altid_controls += "on\n"
    else:
      altid_controls += "off\n"
  return altid_controls

def _build_ss_bond_list(bond_proxies, i_seq_name_hash):
  """Find disulfide bonds (CYS SG-SG pairs) from geometry restraints.

  Returns a list of (i_seq_0, i_seq_1) tuples for each SS bond."""
  ss_bonds = []
  for bp in bond_proxies.simple:
    info_0 = i_seq_name_hash.get(bp.i_seqs[0])
    info_1 = i_seq_name_hash.get(bp.i_seqs[1])
    if info_0 is None or info_1 is None:
      continue
    if (info_0['name'].strip() == "SG" and info_1['name'].strip() == "SG" and
        info_0['resname'].strip() == "CYS" and info_1['resname'].strip() == "CYS"):
      ss_bonds.append((bp.i_seqs[0], bp.i_seqs[1]))
  return ss_bonds

def _same_residue(info_a, info_b):
  """Check whether two atom info dicts belong to the same residue."""
  return (info_a['chain_id'] == info_b['chain_id'] and
          info_a['resseq'] == info_b['resseq'] and
          info_a['icode'] == info_b['icode'])

def _build_bond_hash(bond_proxies, i_seq_name_hash):
  """Build a hash mapping atom i_seq to bonded atom i_seqs within the same
  residue. Shared by make_multikin and export_molprobity_result_as_kinemage."""
  quick_bond_hash = {}
  for bp in bond_proxies.simple:
    if _same_residue(i_seq_name_hash[bp.i_seqs[0]],
                     i_seq_name_hash[bp.i_seqs[1]]):
      if quick_bond_hash.get(bp.i_seqs[0]) is None:
        quick_bond_hash[bp.i_seqs[0]] = []
      quick_bond_hash[bp.i_seqs[0]].append(bp.i_seqs[1])
  return quick_bond_hash

def _build_kinemage(hierarchy, bond_hash, i_seq_name_hash, pdbID,
                    rot_outliers, rama_result, cb_result,
                    restraints_result, keep_hydrogens,
                    ss_bonds=None, sites_cart=None,
                    ss_annotation=None,
                    omega_result=None, rna_puckers_result=None,
                    cablam_result=None):
  """Shared logic for building the kinemage string.

  Args:
    hierarchy: PDB hierarchy
    bond_hash: mapping of atom i_seq to bonded i_seqs (from _build_bond_hash)
    i_seq_name_hash: mapping of i_seq to pdb_label_columns (from build_name_hash)
    pdbID: structure identifier string
    rot_outliers: rotalyze result object
    rama_result: ramalyze result object
    cb_result: cbetadev result object
    restraints_result: mmtbx.validation.restraints.combined result object
    keep_hydrogens: whether to keep hydrogens for probe dots
    ss_bonds: list of (i_seq_0, i_seq_1) tuples for disulfide bonds
    sites_cart: Cartesian coordinates for all atoms
    ss_annotation: iotbx.pdb.secondary_structure.annotation object, or None
    omega_result: omegalyze result object, or None (computed on-the-fly if None)
    rna_puckers_result: rna_puckers result object, or None
    cablam_result: cablamalyze result object, or None
  """
  kin_out = get_default_header()
  altid_controls = get_altid_controls(hierarchy=hierarchy)
  if altid_controls != "":
    kin_out += altid_controls
  kin_out += "@group {%s} dominant animate\n" % pdbID
  initiated_chains = []
  validated_chains = []
  counter = 0
  for model in hierarchy.models():
    for chain in model.chains():
      if chain.id not in initiated_chains:
        kin_out += "@subgroup {%s} dominant master= {chain %s}\n" % (
                  pdbID,
                  chain.id)
        initiated_chains.append(chain.id)
      kin_out += get_kin_lots(chain=chain,
                              bond_hash=bond_hash,
                              i_seq_name_hash=i_seq_name_hash,
                              pdbID=pdbID,
                              index=counter,
                              ss_bonds=ss_bonds,
                              sites_cart=sites_cart)
      # Validation overlays filter by chain_id, so they only need to be
      # emitted once per unique chain ID (not once per chain segment).
      if chain.id not in validated_chains:
        if (chain.is_protein()):
          kin_out += rot_outliers.as_kinemage(chain_id=chain.id,
            pdb_hierarchy=hierarchy)
          kin_out += rama_result.as_kinemage(chain_id=chain.id)
        kin_out += restraints_result.as_kinemage(chain_id=chain.id)
        if (chain.is_protein()):
          kin_out += cb_result.as_kinemage(chain_id=chain.id)
        if rna_puckers_result is not None:
          kin_out += rna_puckers_result.as_kinemage(chain_id=chain.id,
            pdb_hierarchy=hierarchy)
        validated_chains.append(chain.id)
      counter += 1
  if omega_result is not None:
    kin_out += omega_result.as_kinemage()
  else:
    kin_out += omegalyze.omegalyze(pdb_hierarchy=hierarchy,nontrans_only=True,
      out=None,quiet=False).as_kinemage()
  if cablam_result is not None:
    # Suppress self.out write for assembled kinemage; use returned string only
    saved_out = cablam_result.out
    cablam_result.out = None
    kin_out += cablam_result.as_kinemage()
    cablam_result.out = saved_out
  # Ribbon rendering
  try:
    from mmtbx.kinemage.ribbon_rendering import (
        build_secondary_structure_map, consolidate_sheets,
        generate_chain_ribbons)
    from mmtbx.kinemage.ribbons import chain_has_DNA, chain_has_RNA

    ss_map = build_secondary_structure_map(hierarchy, annotation=ss_annotation)
    consolidate_sheets(ss_map)

    ribbon_kin = ""
    ribbon_counter = 0
    for model in hierarchy.models():
      has_dna = any(chain_has_DNA(c) for c in model.chains())
      has_rna = any(chain_has_RNA(c) for c in model.chains())
      for chain in model.chains():
        chain_color = get_chain_color(ribbon_counter)
        ribbon_kin += generate_chain_ribbons(
            chain=chain, secondary_structure=ss_map,
            chain_id=chain.id, chain_color=chain_color,
            has_dna=has_dna, has_rna=has_rna)
        ribbon_counter += 1

    if ribbon_kin:
      kin_out += ribbon_kin
  except Exception as e:
    import sys
    print("Warning: ribbon rendering failed: %s" % str(e), file=sys.stderr)
  kin_out += make_probe_dots(hierarchy=hierarchy, keep_hydrogens=keep_hydrogens)
  kin_out += get_footer()
  return kin_out

def make_multikin(f, processed_pdb_file, pdbID=None, keep_hydrogens=False):
  if pdbID is None:
    pdbID = "PDB"
  hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  i_seq_name_hash = build_name_hash(pdb_hierarchy=hierarchy)
  sites_cart = processed_pdb_file.all_chain_proxies.sites_cart
  geometry = processed_pdb_file.geometry_restraints_manager()
  flags = geometry_restraints.flags.flags(default=True)
  pair_proxies = geometry.pair_proxies(flags=flags,
                                       sites_cart=sites_cart)
  bond_proxies = pair_proxies.bond_proxies
  quick_bond_hash = _build_bond_hash(bond_proxies, i_seq_name_hash)
  ss_bonds = _build_ss_bond_list(bond_proxies, i_seq_name_hash)

  # Run validators
  rot_outliers = rotalyze(pdb_hierarchy=hierarchy, outliers_only=True)
  cb = cbetadev(pdb_hierarchy=hierarchy, outliers_only=True)
  rama = ramalyze(pdb_hierarchy=hierarchy, outliers_only=True)
  omega = omegalyze.omegalyze(pdb_hierarchy=hierarchy, nontrans_only=True,
    out=None, quiet=True)

  # RNA puckers (only if RNA is present)
  rna_puckers_result = None
  has_rna = any(chain.is_na() for model in hierarchy.models()
                for chain in model.chains())
  if has_rna:
    rna_puckers_result = rna_validate.rna_puckers(pdb_hierarchy=hierarchy)

  # CaBLAM
  from mmtbx.validation.cablam import cablamalyze
  cablam_result = None
  has_protein = any(chain.is_protein() for model in hierarchy.models()
                    for chain in model.chains())
  if has_protein:
    from libtbx.utils import null_out
    cablam_result = cablamalyze(
      pdb_hierarchy=hierarchy,
      outliers_only=True,
      out=null_out(),
      quiet=True)

  # Build restraints validation using mmtbx.validation.restraints
  from mmtbx.validation.restraints import combined as restraints_combined
  xray_structure = processed_pdb_file.xray_structure()
  restraints_result = restraints_combined(
    pdb_hierarchy=hierarchy,
    xray_structure=xray_structure,
    geometry_restraints_manager=geometry,
    ignore_hd=True,
    outliers_only=True)

  # Compute secondary structure annotation for ribbon rendering
  try:
    import mmtbx.secondary_structure
    sec_str_from_pdb_file = None
    if hasattr(processed_pdb_file, 'all_chain_proxies'):
      pdb_inp = processed_pdb_file.all_chain_proxies.pdb_inp
      if hasattr(pdb_inp, 'extract_secondary_structure'):
        sec_str_from_pdb_file = pdb_inp.extract_secondary_structure()
    ss_params = mmtbx.secondary_structure.manager.get_default_ss_params()
    ss_params.secondary_structure.protein.search_method = "from_ca"
    ss_params = ss_params.secondary_structure
    from libtbx.utils import null_out
    ss_manager = mmtbx.secondary_structure.manager(
        hierarchy,
        params=ss_params,
        sec_str_from_pdb_file=sec_str_from_pdb_file,
        log=null_out())
    ss_annotation = ss_manager.actual_sec_str
  except Exception:
    ss_annotation = None

  kin_out = _build_kinemage(
    hierarchy=hierarchy,
    bond_hash=quick_bond_hash,
    i_seq_name_hash=i_seq_name_hash,
    pdbID=pdbID,
    rot_outliers=rot_outliers,
    rama_result=rama,
    cb_result=cb,
    restraints_result=restraints_result,
    keep_hydrogens=keep_hydrogens,
    ss_bonds=ss_bonds,
    sites_cart=sites_cart,
    ss_annotation=ss_annotation,
    omega_result=omega,
    rna_puckers_result=rna_puckers_result,
    cablam_result=cablam_result)

  outfile = open(f, 'w')
  outfile.write(kin_out)
  outfile.close()
  return f

def usage():
  return """
phenix.kinemage file.pdb [params.eff] [options ...]

Options:

  pdb=input_file        input PDB file
  cif=cif_file          input custom definitions (ligands, etc.)
  keep_hydrogens=False  keep input hydrogen files (otherwise regenerate)

Example:

  phenix.kinemage pdb=1ubq.pdb cif=ligands.cif
"""

def run(args, pdb_interpretation_params=None):
  if (len(args) == 0 or "--help" in args or "--h" in args or "-h" in args):
      raise Usage(usage())
  pdbID = None
  master_phil = get_master_phil()
  import iotbx.phil
  input_objects = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    pdb_file_def="kinemage.pdb",
    cif_file_def="kinemage.cif")
  work_params = input_objects.work.extract()
  # Stange workaround?
  # auto_cdl=True
  # for po in input_objects["phil"]:
  #   if po.as_str().find(".cdl")>-1:
  #     auto_cdl=False
  #     break
  # if auto_cdl:
  work_params.kinemage.pdb_interpretation.restraints_library.cdl = Auto
  if work_params.kinemage.pdb == None:
    assert len(input_objects["pdb"]) == 1
    file_obj = input_objects["pdb"][0]
    file_name = file_obj.file_name
  else:
    file_name = work_params.kinemage.pdb
  if file_name and os.path.exists(file_name):
    pdb_io = pdb.input(file_name)
    pdbID = os.path.basename(pdb_io.source_info().split(' ')[1]).split('.')[0]
  else :
    raise Sorry("PDB file does not exist")
  assert pdb_io is not None
  cif_file = None
  cif_object = None
  cif_file = work_params.kinemage.cif
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  if cif_file != None:
    for cif in cif_file:
      try:
        cif_object = monomer_library.server.read_cif(file_name=cif)
      except Exception:
        raise Sorry("Unknown file format: %s" % show_string(cif))
    if cif_object != None:
      for srv in [mon_lib_srv, ener_lib]:
        srv.process_cif_object(cif_object=cif_object)
  if pdb_interpretation_params is None:
    #pdb_int_work_params = pdb_interpretation.master_params.extract()
    pdb_int_work_params = work_params.kinemage.pdb_interpretation
  else:
    pdb_int_work_params = pdb_interpretation_params
  pdb_int_work_params.clash_guard.nonbonded_distance_threshold = None
  processed_pdb_file = pdb_interpretation.process(
        mon_lib_srv=mon_lib_srv,
        ener_lib=ener_lib,
        pdb_inp=pdb_io,
        params=pdb_int_work_params,
        substitute_non_crystallographic_unit_cell_if_necessary=True)
  if work_params.kinemage.out_file is not None:
    outfile = work_params.kinemage.out_file
  else :
    outfile = pdbID+'.kin'
  outfile = make_multikin(f=outfile,
                          processed_pdb_file=processed_pdb_file,
                          pdbID=pdbID,
                          keep_hydrogens=work_params.kinemage.keep_hydrogens)
  return outfile

def export_molprobity_result_as_kinemage(
    result,
    pdb_hierarchy,
    geometry,
    probe_file,
    keep_hydrogens=False,
    pdbID="PDB",
    ss_annotation=None,
    cablam_result=None):
  assert (result.restraints is not None)
  i_seq_name_hash = build_name_hash(pdb_hierarchy=pdb_hierarchy)
  sites_cart = pdb_hierarchy.atoms().extract_xyz()
  flags = geometry_restraints.flags.flags(default=True)
  pair_proxies = geometry.pair_proxies(flags=flags,
                                       sites_cart=sites_cart)
  bond_proxies = pair_proxies.bond_proxies
  quick_bond_hash = _build_bond_hash(bond_proxies, i_seq_name_hash)
  ss_bonds = _build_ss_bond_list(bond_proxies, i_seq_name_hash)
  rna_puckers_result = None
  if hasattr(result, 'rna') and result.rna is not None:
    rna_puckers_result = result.rna.puckers
  omega_result = None
  if hasattr(result, 'omegalyze'):
    omega_result = result.omegalyze
  return _build_kinemage(
    hierarchy=pdb_hierarchy,
    bond_hash=quick_bond_hash,
    i_seq_name_hash=i_seq_name_hash,
    pdbID=pdbID,
    rot_outliers=result.rotalyze,
    rama_result=result.ramalyze,
    cb_result=result.cbetadev,
    restraints_result=result.restraints,
    keep_hydrogens=keep_hydrogens,
    ss_bonds=ss_bonds,
    sites_cart=sites_cart,
    ss_annotation=ss_annotation,
    omega_result=omega_result,
    rna_puckers_result=rna_puckers_result,
    cablam_result=cablam_result)
