import os, sys, string
from iotbx import pdb
from mmtbx import monomer_library
from mmtbx.chemical_components import get_bond_pairs
from mmtbx.refinement import fit_rotamers
from mmtbx.validation.rotalyze import rotalyze
from mmtbx.validation.ramalyze import ramalyze
from mmtbx.validation.cbetadev import cbetadev
from iotbx.pdb import common_residue_names_get_class
from libtbx import easy_run

def get_residue_bonds(residue):
  if residue is None: return []
  if residue is False: return []
  if residue in never_do_residues: return []
  monomer_lib_entry = mon_lib_query(residue)
  ml = mon_lib_query(residue)
  if ml is None: return []
  bonds = []
  for bond in ml.bond_list:
    bonds.append([bond.atom_id_1, bond.atom_id_2])
  return bonds

def kin_vec(start_key, start_xyz, end_key, end_xyz):
  return "{%s} P %.3f %.3f %.3f {%s} L %.3f %.3f %.3f\n" % (
           start_key,
           start_xyz[0],
           start_xyz[1],
           start_xyz[2],
           end_key,
           end_xyz[0],
           end_xyz[1],
           end_xyz[2])

def make_probe_dots(hierarchy, keep_hydrogens=False):
  probe = 'phenix.probe -4H -quiet -noticks -nogroup -dotmaster -mc -self "alta" -'
  trim = "phenix.reduce -quiet -trim -"
  build = "phenix.reduce -oh -his -flip -pen9999 -keep -allalt -"
  probe_return = ""
  for i,m in enumerate(hierarchy.models()):
    r = pdb.hierarchy.root()
    mdc = m.detached_copy()
    r.append_model(mdc)
    if keep_hydrogens is False:
      clean_out = easy_run.fully_buffered(trim,
                                  stdin_lines=r.as_pdb_string())
      build_out = easy_run.fully_buffered(build,
                                  stdin_lines=clean_out.stdout_lines)
      input_str = string.join(build_out.stdout_lines, '\n')
    else:
      input_str = r.as_pdb_string()
    probe_out = easy_run.fully_buffered(probe,
                         stdin_lines=input_str).stdout_lines
    for line in probe_out:
      probe_return += line+'\n'
  return probe_return

def cbeta_dev(chain, pdbID, deviations, ideal):
  cbeta_out = "@subgroup {CB dev} dominant\n"
  cbeta_out += "@balllist {CB dev Ball} color= gold radius= 0.0020   master= {Cbeta dev}\n"
  outlier_list = []
  angle_dict = {}
  deviation_dict = {}
  for outlier in deviations.splitlines():
    if outlier.startswith('pdb:alt:res'):
      continue
    PDBfileStr,altchar,res,sub,resnum,dev,dihedralNABB,occ,segid,last = outlier.split(':')
    key = '%s%4s %s' % \
             (sub.strip(),
             resnum.rstrip(),
             altchar+res.upper())
    outlier_list.append(key)
    angle_dict[key] = float(dihedralNABB)
    deviation_dict[key] = float(dev)
  for residue_group in chain.residue_groups():
    for atom_group in residue_group.atom_groups():
      altloc = atom_group.altloc
      if len(altloc) < 1:
        altloc = " "
      check_key = '%s%4s %s' % \
                    (chain.id,
                     residue_group.resseq,
                     altloc+atom_group.resname.strip())
      if check_key not in outlier_list:
        continue
      for atom in atom_group.atoms():
        if atom.name == ' CB ':
          altloc = atom_group.altloc
          if len(altloc) < 1:
            altloc = " "
          chainid = chain.id
          if len(chainid) == 1:
            chainid = " "+chainid
          ideal_key = altloc+atom_group.resname.lower()+ \
                      chainid+residue_group.resseq+ "  " 
          ideal_xyz = ideal[ideal_key]
          key = "%s %s %s%s  %.3f %.2f" % (
              atom.name.lower(),
              atom_group.resname.lower(),
              chain.id,
              residue_group.resseq,
              deviation_dict[check_key],
              angle_dict[check_key])
          cbeta_out += '{%s} r=%.3f magenta  %.3f, %.3f, %.3f\n' % (
              key,
              deviation_dict[check_key],
              ideal_xyz[0],
              ideal_xyz[1],
              ideal_xyz[2])
  if len(cbeta_out.splitlines()) == 2:
    cbeta_out = ""
  return cbeta_out

def rotamer_outliers(chain, pdbID, rot_outliers):
  mc_atoms = ["N", "C", "O", "OXT"]
  rot_out = "@subgroup {rotamer outliers} dominant\n"
  rot_out += "@vectorlist {chain %s} color= gold  master= {rotamer outlie}\n" % chain.id
  outlier_list = []
  for outlier in rot_outliers.splitlines():
    outlier_list.append(outlier.split(':')[0])
  for residue_group in chain.residue_groups():
    for atom_group in residue_group.atom_groups():
      check_key = '%s%4s %s' % \
                    (chain.id,
                     residue_group.resseq,
                     atom_group.altloc+atom_group.resname.strip())
      if check_key not in outlier_list:
        continue
      key_hash = {}
      xyz_hash = {}
      for atom in atom_group.atoms():      
        key = "%s %s %s%s  B%.2f %s" % (
              atom.name.lower(),
              atom_group.resname.lower(),
              chain.id,
              residue_group.resseq,
              atom.b,
              pdbID)
        key_hash[atom.name.strip()] = key
        xyz_hash[atom.name.strip()] = atom.xyz
      bonds = get_bond_pairs(code=atom_group.resname)
      for bond in bonds:
        if bond[0] in mc_atoms or bond[1] in mc_atoms:
          continue
        elif bond[0].startswith('H') or bond[1].startswith('H'):
          continue
        rot_out += kin_vec(key_hash[bond[0]],
                           xyz_hash[bond[0]],
                           key_hash[bond[1]],
                           xyz_hash[bond[1]])
  if len(rot_out.splitlines()) == 2:
    rot_out = ""
  return rot_out
          
def get_chain_color(index):
  chain_colors = ['white',
                  'yellowtint',
                  'peachtint',
                  'pinktint',
                  'lilactint',
                  'bluetint',
                  'greentint']
  match = False
  while not match:
    if index > 6:
      index = index - 6
    else:
      match = True
  return chain_colors[index]
                  

def get_kin_lots(chain, pdbID=None, index=0, show_hydrogen=True):
  mc_atoms = ["N", "CA", "C", "O", "OXT"]
  mc_veclist = ""
  sc_veclist = ""
  mc_h_veclist = ""
  sc_h_veclist = ""
  ca_trace = ""
  water_list = ""
  kin_out = ""
  color = get_chain_color(index)
  mc_veclist = "@vectorlist {mc} color= %s  master= {mainchain}\n" % color
  sc_veclist = "@vectorlist {sc} color= cyan  master= {sidechain}\n"
  ca_trace = "@vectorlist {Calphas} color= %s master= {Calphas}\n" % color
  water_list = "@balllist {water O} color= peachtint  radius= 0.15  master= {water}\n"
  if show_hydrogen:
    mc_h_veclist = \
      "@vectorlist {mc H} color= gray nobutton master= {mainchain} master= {H's}\n"
    sc_h_veclist = \
      "@vectorlist {sc H} color= gray nobutton master= {sidechain} master= {H's}\n"
  prev_resid = None
  prev_C_xyz = None
  prev_C_key = None
  prev_CA_xyz = None
  prev_CA_key = None
  cur_resid = None
  cur_C_xyz = None
  cur_C_key = None
  cur_CA_xyz = None
  cur_CA_key = None
  for residue_group in chain.residue_groups():
    cur_resid = residue_group.resseq
    for atom_group in residue_group.atom_groups():
      key_hash = {}
      xyz_hash = {}
      for atom in atom_group.atoms():      
        key = "%s %s %s%s  B%.2f %s" % (
              atom.name.lower(),
              atom_group.resname.lower(),
              chain.id,
              residue_group.resseq,
              atom.b,
              pdbID)
        key_hash[atom.name.strip()] = key
        xyz_hash[atom.name.strip()] = atom.xyz
        if atom.name == ' C  ':
          cur_C_xyz = atom.xyz
          cur_C_key = key
        if atom.name == ' CA ':
          cur_CA_key = key
          cur_CA_xyz = atom.xyz
          if prev_CA_key != None and prev_CA_xyz != None:
            if int(residue_group.resid()) - int(prev_resid) == 1:
              try:
                ca_trace += kin_vec(prev_CA_key, prev_CA_xyz, key, atom.xyz)
              except:
                continue
        if atom.name == ' N  ':
          if prev_C_key != None and prev_C_xyz != None:
            if int(residue_group.resid()) - int(prev_resid) == 1:
              try:
                mc_veclist += kin_vec(prev_C_key, prev_C_xyz, key, atom.xyz)
              except:
                continue
        if atom_group.resname.lower() == 'hoh':
          if atom.name == ' O  ':
            water_list += "{%s} P %.3f %.3f %.3f\n" % (
                       key,
                       atom.xyz[0],
                       atom.xyz[1],
                       atom.xyz[2])
      bonds = get_bond_pairs(code=atom_group.resname)

      prev_CA_xyz = cur_CA_xyz
      prev_CA_key = cur_CA_key
      prev_C_xyz = cur_C_xyz
      prev_C_key = cur_C_key
      prev_resid = cur_resid
      
      for bond in bonds:
        if bond[0] in mc_atoms and bond[1] in mc_atoms:
          try:
            mc_veclist += kin_vec(key_hash[bond[0]],
                                  xyz_hash[bond[0]],
                                  key_hash[bond[1]],
                                  xyz_hash[bond[1]]) 
          except:
            continue

        elif (bond[0].startswith('H') or bond[1].startswith('H')):
          if show_hydrogen:
            if (bond[0] in mc_atoms or bond[1] in mc_atoms): 
              try:
                mc_h_veclist += kin_vec(key_hash[bond[0]],
                                        xyz_hash[bond[0]],
                                        key_hash[bond[1]],
                                        xyz_hash[bond[1]])
              except:
                continue
            else:
              try:
                sc_h_veclist += kin_vec(key_hash[bond[0]],
                                        xyz_hash[bond[0]],
                                        key_hash[bond[1]],
                                        xyz_hash[bond[1]])
              except:
                continue
        else:
          try:
            sc_veclist += kin_vec(key_hash[bond[0]],
                                  xyz_hash[bond[0]],
                                  key_hash[bond[1]],
                                  xyz_hash[bond[1]])
          except:
            continue
  
  #clean up empty lists:
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
  if len(water_list.splitlines()) > 1:
    kin_out += water_list
  return kin_out
              
def get_default_header():
  header = """
@kinemage 1
@onewidth
@1viewid {Overview}
@master {mainchain} indent
@master {Calphas} indent
@master {sidechain} indent
@master {H's} indent
@pointmaster 'H' {H's}
@master {water} indent
"""
  return header

def get_footer():
  footer = """
@master {mainchain} off
@master {sidechain} off
@master {H's} off
@master {water} off
@master {Calphas} on
@master {vdw contact} off
@master {small overlap} off
@master {H-bonds} off
@master {Cbeta dev} on
"""
  return footer

def kingen(f, pdb_io):
  pdbID = None
  pdbID = os.path.basename(pdb_io.source_info().split(' ')[1]).split('.')[0]
  assert pdb_io is not None
  if(pdb_io is not None):
    hierarchy = pdb_io.construct_hierarchy()
  kin_out = get_default_header()
  kin_out += "@group {%s} dominant animate\n" % pdbID
  initiated_chains = []
  rt = rotalyze()
  rot_outliers, output_list = rt.analyze_pdb(hierarchy=hierarchy, outliers_only=True)
  cb =cbetadev()
  deviations, summary, output_list = cb.analyze_pdb(hierarchy=hierarchy,
                                                           outliers_only=True)
  counter = 0
  for model in hierarchy.models():
    for chain in model.chains():
      if chain.id not in initiated_chains:
        kin_out += "@subgroup {%s} dominant master= {chain %s}\n" % (
                  pdbID,
                  chain.id)
        initiated_chains.append(chain.id)
      kin_out += get_kin_lots(chain=chain, pdbID=pdbID, index=counter)
      kin_out += rotamer_outliers(chain=chain, pdbID=pdbID, rot_outliers=rot_outliers)
      kin_out += cbeta_dev(chain=chain, 
                           pdbID=pdbID,
                           deviations=deviations,
                           ideal=cb.get_beta_ideal())
      counter += 1
  kin_out += make_probe_dots(hierarchy=hierarchy)
  kin_out += get_footer()
  outfile = file(f, 'w')
  for line in kin_out:
    outfile.write(line)
  outfile.close()

def run(args):
  assert len(args)==1
  file_name = args[0]
  if file_name and os.path.exists(file_name):
      pdb_io = pdb.input(file_name)
  pdbID = os.path.basename(pdb_io.source_info().split(' ')[1]).split('.')[0]
  outfile = pdbID+'.kin'
  kingen(f=outfile, pdb_io=pdb_io)

if __name__ == "__main__":
  run(args=sys.argv[1:])