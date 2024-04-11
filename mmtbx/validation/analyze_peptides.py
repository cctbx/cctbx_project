from __future__ import absolute_import, division, print_function
import sys
import mmtbx.rotamer
from libtbx.utils import Sorry

def get_master_phil():
  import libtbx.phil
  return libtbx.phil.parse(
    input_string="""
    analyze_peptides {
      pdb = None
        .type = path
        .multiple = True
        .help = '''Enter a PDB file name'''
      cis_only = True
        .type = bool
        .help = '''Only print cis peptides'''
      verbose = True
        .type = bool
        .help = '''Verbose'''
}
""")

def usage():
  return """
mmtbx.analyze_peptides file.pdb [params.eff] [options ...]

Options:

  pdb=input_file        input PDB file
  cis_only=True         only list cis peptides
  verbose=True          verbose text output

Example:

  mmtbx.analyze_peptides pdb=1xwl.pdb

"""

def analyze(pdb_hierarchy):
  cis_peptides = []
  trans_peptides = []
  outliers = []
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      prev_rg = None
      chain_id = chain.id
      for rg in chain.residue_groups():
        if prev_rg is None:
          prev_rg = rg
          continue
        for conformer in rg.conformers():
          if conformer.is_protein():
            altloc = conformer.altloc
            for residue in conformer.residues():
              atoms = residue.atoms()
              resname = residue.resname
              resid = residue.resid()
              #get appropriate prev_atoms
              for prev_con in prev_rg.conformers():
                if prev_con.is_protein():
                  if( (prev_con.altloc == altloc) or
                      (prev_con.altloc == ' ' or
                       altloc == ' ') ):
                    prev_altloc = prev_con.altloc
                    for prev_res in prev_con.residues():
                      prev_atoms = prev_res.atoms()
                      prev_resname = prev_res.resname
                      prev_resid = prev_res.resid()
                      is_cis = is_cis_peptide(prev_atoms=prev_atoms,
                                              atoms=atoms)
                      is_trans = is_trans_peptide(prev_atoms=prev_atoms,
                                                  atoms=atoms)

                      res1 = "%s%5s %s%s" % (chain.id,
                                             prev_resid,
                                             prev_altloc,
                                             prev_resname)
                      res2 = "%s%5s %s%s" % (chain.id,
                                             resid,
                                             altloc,
                                             resname)
                      if is_cis:
                        cis_peptides.append( (res1, res2) )
                      elif is_trans:
                        trans_peptides.append( (res1, res2) )
                      elif ( (not is_cis) and
                             (not is_trans) ):
                        outliers.append( (res1, res2) )
                      else:
                        raise Sorry("%s:%s " % (res1, res2) +
                                    "have uncharacteristic peptide geometry")
        prev_rg = rg
  return cis_peptides, trans_peptides, outliers

def get_omega(prev_atoms, atoms):
  prevCA, prevC, curN, curCA = None, None, None, None
  if (prev_atoms is not None):
    for atom in prev_atoms:
      if (atom.name == " CA "): prevCA = atom
      if (atom.name == " C  "): prevC = atom
  if (atoms is not None):
    for atom in atoms:
      if (atom.name == " N  "): curN = atom
      if (atom.name == " CA "): curCA = atom
  if (prevCA is not None and
      prevC is not None and
      curN is not None and
      curCA is not None):
    return mmtbx.rotamer.omega_from_atoms(prevCA, prevC, curN, curCA)

def is_cis_peptide(prev_atoms, atoms):
  omega = get_omega(prev_atoms, atoms)
  if(omega > -30 and omega < 30):
    return True
  else:
    return False

def is_trans_peptide(prev_atoms, atoms):
  omega = get_omega(prev_atoms, atoms)
  if( (omega > 150 and omega <= 180) or
      (omega >= -180 and omega < -150) ):
    return True
  else:
    return False

def show(
      cis_peptides,
      trans_peptides,
      outliers,
      log=None,
      cis_only=True):
  print_list = []
  for pair in cis_peptides:
    txt = "%s:%s:cis" % (pair[0], pair[1])
    print_list.append(txt)
  if not cis_only:
    for pair in trans_peptides:
      txt = "%s:%s:trans" % (pair[0], pair[1])
      print_list.append(txt)
    for pair in outliers:
      txt = "%s:%s:OUTLIER" % (pair[0], pair[1])
      print_list.append(txt)
  def get_sort_key(pair):
    return pair[0:6]
  print_list.sort(key=get_sort_key)
  if log is None:
    log = sys.stdout
  print("#\n#peptide geometry analysis\n#", file=log)
  print("residue1:residue2:classification", file=log)
  for line in print_list:
    print(line, file=log)
