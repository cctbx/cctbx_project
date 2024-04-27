from __future__ import absolute_import, division, print_function
#!/usr/bin/python
# remediator.py - version 1.61.110622 6/22/11
# Copyright 2007-2011, Jeffrey J. Headd and Robert Immormino

# revision 1.56 - JJH 070808 - added support for DU DNA base
#               - JJH 070808 - added compiled RE object for HN2 RES special case
#       - JJH 070815 - updated name of hash dictionary file
#       - JJH 070823 - added support for CNS Xplor and Coot RNA names
#       - JJH 070908 - added REMARK   4 comment addition
#               - JJH 070913 - added support for left-justified RNA/DNA old names
#       - JJH 070913 - added support for all left-justified residue names
#       - JJH 080328 - fixed REMARK   4 comment addition to work w/ PHENIX header info
# revision 1.57 - vbc 080620 - lines get output with at least 80 columns
#                              master_hash.txt not found bug fixed, must be in same directory
#                              as this script
# revision 1.58 - JJH 080630 - fixed --old functionality for DNA, this did not work previously
#                              and would sometimes incorrectly rename protein atoms
#               - JJH 080827 - fixed a bug in this version only by restricting the number of
#                              substitutions to 1 per line
# revision 1.59 - JJH 081015 - added a --dict flag, which allows user to input a custom
#                              dictionary for detecting non-standard atom names and converting
#                              them to version 3.0
# revision 1.60 - JJH 081120 - reorganized code into functions and cleaned up flag usage
# revision 1.61 - JJH 110622 - moved to iotbx
# revision 1.70 - VBC 201116 - switch over to using model object and chemical components library
# revision 1.71 - VBC 230923 - switch to using atom.parent() instead of relying on id_str()

import sys
import os
import re
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate

import libtbx.load_env
import iotbx.pdb
import iotbx.pdb.mmcif
import mmtbx.model
from libtbx.utils import null_out

import mmtbx.chemical_components as chemical_components

class Program(ProgramTemplate):

  description = """
Program for converting a PDB format file to version 3.2 or 2.3.

iotbx.pdb_remediator model.pdb [params.eff] [options ...]

Options:

  file_name     specify input pdb file
  output_file   specify output file name (defaults to stdout)
  version       3.2 (default) or 2.3
  dict          optionally supply custom definition file

Examples:

Convert version 2.3 file to version 3.2 naming:

  iotbx.pdb_remediator model.pdb > model_3.2.pdb
  iotbx.pdb_remediator file_name=model.pdb output_file=model_3.2.pdb_file

Convert version 3.2 file to version 2.3 naming:

  iotbx.pdb_remediator model.pdb version=2.3 > model_2.3.pdb_file
"""
  datatypes = ['model']
  results = ""

  def run(self):
    model = self.data_manager.get_model()
    convert_to_new = True
    if self.params.version == 2.3:
      convert_to_new = False
    remed = Remediator()
    self.results = remed.remediate_model_object(model, convert_to_new)

  def get_results(self):
    return self.results

  def validate(self):
    if not self.params.file_name or (not os.path.isfile(self.params.file_name)):
      raise Sorry("The file %s is missing" %(self.params.file_name))


#{{{ pre_screen_file
def pre_screen_file(filename, atom_exch, alt_atom_exch):
  count = 0
  res_count = 0
  if not os.path.isfile(filename):
    from libtbx.utils import Sorry
    raise Sorry("The file %s is missing" %(filename))
  pdb_file = open(filename)
  for line in pdb_file:
    line=line.rstrip()
    line=line.ljust(80)
    type_test = line[0:6]
    if type_test in ("ATOM  ", "HETATM", "TER   ", "ANISOU", "SIGATM", "SIGUIJ", "LINK  "):
      adjust_res = False
      #--make any left-justified residue names right-justified------------------
      if re.match(r'.{17}([a-zA-Z0-9])  ',line):
        line = re.sub(r'\A(.{17})(.)\s\s',r'\g<1>  \g<2>',line)
        adjust_res = True
      elif re.match(r'.{17}([a-zA-Z0-9][a-zA-Z0-9]) ',line):
        line = re.sub(r'\A(.{17})(..)\s',r'\g<1> \g<2>',line)
        adjust_res = True
      #-------------------------------------------------------------------------

      #--pre-screen for CNS Xplor RNA base names and Coot RNA base names--------
      if re.match(r'.{17}(GUA|ADE|CYT|THY|URI)',line):
        line = re.sub(r'\A(.{17})(.)..',r'\g<1>  \g<2>',line)
        adjust_res = True
      elif re.match(r'.{17}(OIP| Ar| Gr| Cr| Ur)',line):
        line = re.sub(r'\A(.{17}).(.).',r'\g<1>  \g<2>',line)
        adjust_res = True
      #-------------------------------------------------------------------------

      if adjust_res:
        res_count += 1
      entry = line[12:20]
      clean_entry = entry[0:4] + " " + entry[5:8]
      if clean_entry in atom_exch:
        if clean_entry in alt_atom_exch:
          pass
        else:
          count += 1
  return count, res_count
#}}}

#{{{ build_hash
#--Build Hash Table------------------------------------------------
def build_hash(remediated_out, custom_dict, user_dict):
  atom_exch = {}
  file_name = os.path.join(libtbx.env.dist_path("iotbx"), "pdb",
    "remediation", "remediation.dict")
  f = open(file_name, "r")
  if remediated_out == True: #converting to remediated
    for line in f:
      line=line.rstrip()
      new, old = line.split(':')
      atom_exch[old] = new
  else: #converting to old
    for line in f:
      new, old = line.split(':')
      atom_exch[new] = old
  if custom_dict == True:
    user_f = open(user_dict)
    if remediated_out == True: #converting to remediated
      for line in user_f:
        line=line.rstrip()
        new, old = line.split(':')
        atom_exch[old] = new
    else:
      for line in user_f:
        line=line.rstrip()
        new, old = line.split(':')
        atom_exch[new] = old
    user_f.close()
  f.close()
  return atom_exch
#}}}
#------------------------------------------------------------------

def justify_atom_names(atom_name):
  if (len(atom_name) == 4):
    return atom_name
  if (len(atom_name) == 3):
    if (atom_name == "H3T"):
      return atom_name.ljust(4)
    if (atom_name[0:2] == "NH"):
      return atom_name.rjust(4)
    if (atom_name[1] == "H"):
      return atom_name.ljust(4)
    return atom_name.rjust(4)
  else:
    return atom_name.center(4)

def get_model_from_file(file_path):
  pdb_inp = iotbx.pdb.input(file_name = file_path)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    stop_for_unknowns = False,
    log = null_out())
  model.process(make_restraints=False)
  return model

class Remediator():
  amino_acids = [ "ALA","ARG","ASN","ASP","ASX","CSE","CYS","GLN","GLU","GLX","GLY","HIS","ILE",
    "LEU","LYS","MET","MSE","PHE","PRO","SER","THR","TRP","TYR","VAL" ]
  na_bases = ["  A", "  C", "  T", "  G", "  I", "  U"]
  dna_bases = [" DA", " DC", " DT", " DG", " DI", " DU"]
  residues_dict = {}

  def build_hash_from_chem_components(self, residue_name, convert_to_new=True, build_all_atoms=False):
    atom_exch = {}
    old_atom_name_exceptions = ["1H  ", "2H  ", "3H  "]
    new_atom_name_exceptions = [" H1 ", " H2 ", " H3 "]
    if residue_name in self.amino_acids: # covers the case where the N terminal Hs aren't in chem comps
      for new_atom, old_atom in zip(new_atom_name_exceptions, old_atom_name_exceptions):
        if convert_to_new:
          atom_exch[old_atom+" "+residue_name] = new_atom+" "+residue_name
        else:
          atom_exch[new_atom+" "+residue_name] = old_atom+" "+residue_name
    na_old_atom_name_exceptions = ["H5T ","5HO*", "H3T ", "3HO*"] # it's difficult in 2020 to figure out the proper spacing for this old style atom name
    na_new_atom_name_exceptions = ["HO5'","HO5'", "HO3'", "HO3'"] # old remediator uses 5HO* but other examples on the web use H5T for HO5'
    residues_to_test = [ residue_name ]
    if (residue_name in self.na_bases) or (residue_name in self.dna_bases):
      if (residue_name in self.na_bases):
        residues_to_test.append(" D"+residue_name[2])
      else:
        residues_to_test.append("  "+residue_name[2])
      for nucleic_acid in residues_to_test:
        for new_atom, old_atom in zip(na_new_atom_name_exceptions, na_old_atom_name_exceptions):
          if convert_to_new:
            atom_exch[old_atom+" "+nucleic_acid] = new_atom+" "+nucleic_acid
          else:
            atom_exch[new_atom+" "+nucleic_acid] = old_atom+" "+nucleic_acid
    for residue in residues_to_test:
      if (chemical_components.is_code(residue)):
        #sys.stderr.write(residue_name+" is in chem_components\n")
        new_atom_names = chemical_components.get_atom_names(residue, alternate=False)
        old_atom_names = chemical_components.get_atom_names(residue, alternate=True)
        atom_type_symbols = chemical_components.get_atom_type_symbol(residue)
        if (len(new_atom_names) == len(old_atom_names)):
          for new_atom, old_atom, atom_type in zip(new_atom_names, old_atom_names, atom_type_symbols):
            if build_all_atoms or not (new_atom == old_atom):
              for res_name in residues_to_test:
                justified_old_atom = iotbx.pdb.mmcif.format_pdb_atom_name(old_atom, atom_type)
                justified_new_atom = iotbx.pdb.mmcif.format_pdb_atom_name(new_atom, atom_type)
                new_entry = justified_new_atom+" "+res_name
                old_entry = justified_old_atom+" "+res_name
                #print(old_entry + "->" + new_entry)
                if convert_to_new:
                  atom_exch[old_entry] = new_entry
                  #check for 1HA, 2HA, etc, which don't seem to always be in chem components as possible old names
                  if not build_all_atoms or re.match(r' H[A-Z]\d', justified_old_atom):
                    digit_first_hydrogen = justified_old_atom[3]+justified_old_atom[1:3]+" "
                    atom_exch[digit_first_hydrogen+" "+res_name] = new_entry
                else:
                  atom_exch[new_entry] = old_entry
    #print("finished making hash:"+str(atom_exch))
    return atom_exch

  def is_model_v3(self, model):
    self.residues_dict = {}
    pdb_hierarchy = model.get_hierarchy()
    atoms = pdb_hierarchy.atoms()
    non_v3_atoms_count = 0
    for atom in atoms:
      #print(atom.name)
      res_name = atom.parent().resname
      if not res_name in self.residues_dict:
        self.residues_dict[res_name] = self.build_hash_from_chem_components(res_name, convert_to_new=False, build_all_atoms=True)
        #print(res_name+"\n"+str(residues_dict[res_name]))
      atom_exch_dict = self.residues_dict[res_name]
      if not atom.name+" "+res_name in atom_exch_dict:
        #print(atom_exch_dict)
        print("|"+atom.name +"| does not seem to be v3 in "+res_name)
        non_v3_atoms_count=non_v3_atoms_count+1
      #print(atom.name)
    return non_v3_atoms_count == 0

  def remediate_model_object(self, model, convert_to_new):
    self.residues_dict = {}
    non_v3_atoms_count = 0
    pdb_hierarchy = model.get_hierarchy()
    #pdb_hierarchy.atoms().reset_i_seq()
    residue_groups = pdb_hierarchy.residue_groups()
    for residue_group in residue_groups:
      atom_groups = residue_group.atom_groups()
      for atom_group in atom_groups:
        res_name = atom_group.resname
        if res_name in self.dna_bases or res_name in self.na_bases:
          self.remediate_na_atom_group(atom_group, convert_to_new)
          res_name = atom_group.resname
        for atom in atom_group.atoms():
          if not res_name in self.residues_dict:
            self.residues_dict[res_name] = self.build_hash_from_chem_components(res_name, convert_to_new=convert_to_new, build_all_atoms=False)
          atom_exch_dict = self.residues_dict[res_name]
          #print(atom.name + str(atom.name+" "+res_name in atom_exch_dict))
          if atom.name+" "+res_name in atom_exch_dict:
            new_entry = atom_exch_dict.get(atom.name+" "+res_name)
            atom.set_name(new_entry[0:4])
            non_v3_atoms_count=non_v3_atoms_count+1
    pdb_hierarchy.sort_atoms_in_place()
    pdb_hierarchy.atoms_reset_serial()
    return model

  def remediate_na_atom_group(self, atom_group, convert_to_new):
    rna_group = False
    if convert_to_new:
      if atom_group.resname in self.na_bases:
        for atom in atom_group.atoms():
          if atom.name == " O2'" or atom.name == " O2*":
            rna_group = True
        if not rna_group:
          atom_group.resname = " D" + atom_group.resname[2]
    else:
      if atom_group.resname in self.dna_bases:
        atom_group.resname = "  " + atom_group.resname[2]


def remediate_atomic_line(line, atom_exch):
  #--make any left-justified residue names right-justified------------------
  if re.match(r'.{17}([a-zA-Z0-9])  ',line):
    line = re.sub(r'\A(.{17})(.)\s\s',r'\g<1>  \g<2>',line)
  elif re.match(r'.{17}([a-zA-Z0-9][a-zA-Z0-9]) ',line):
    line = re.sub(r'\A(.{17})(..)\s',r'\g<1> \g<2>',line)
  #-------------------------------------------------------------------------

  #--pre-screen for CNS Xplor RNA base names and Coot RNA base names--------
  if re.match(r'.{17}(GUA|ADE|CYT|THY|URI)',line):
    line = re.sub(r'\A(.{17})(.)..',r'\g<1>  \g<2>',line)
  elif re.match(r'.{17}(OIP| Ar| Gr| Cr| Ur)',line):
    line = re.sub(r'\A(.{17}).(.).',r'\g<1>  \g<2>',line)
  #-------------------------------------------------------------------------

  entry = line[12:20]
  clean_entry = entry[0:4] + " " + entry[5:8]
  if clean_entry in atom_exch:
    line = line.replace(clean_entry[0:4],atom_exch[clean_entry][0:4],1)
  return line

#{{{ remediate
#----PDB routine---------------------------------------------------
def remediate(filename, remediated_out, f=None):
  if remediated_out:
    remark4 = "REMARK   4 REMEDIATOR VALIDATED PDB VERSION 3.2 COMPLIANT".ljust(80)
  else:
    remark4 = "REMARK   4 REMEDIATOR VALIDATED PDB VERSION 2.3 COMPLIANT".ljust(80)

  if f == None:
    f = sys.stdout
  previous = None
  current = ""
  print_line = ""
  remark_flag = False
  remark_block = False

  pdb_file = open(filename)

  aa_re = re.compile(
    ' HN2 (ALA|ARG|ASN|ASP|ASX|CSE|CYS|GLN|GLU|GLX|GLY|HIS|ILE|'+
    'LEU|LYS|MET|MSE|PHE|PRO|SER|THR|TRP|UNK|TYR|VAL)')

  for line in pdb_file:
    line=line.rstrip()
    line=line.ljust(80)
    type_test = line[0:6]
    if remark_flag == False:
      if type_test == "REMARK":
        if remark_block == False:
          remark_block = True
        if re.search(remark4,line):
          remark_flag = True
        elif re.match(r'REMARK...\D',line):
          print_line += line + "\n"
          continue
        elif re.match('REMARK   4 REMEDIATOR',line):
          continue
        elif int(line[6:10]) > 4:
          print_line += remark4 + "\n"
          remark_flag = True
        else:
          print_line += line + "\n"
          continue

    if type_test in ("ATOM  ", "HETATM", "TER   ", "ANISOU", "SIGATM", "SIGUIJ", "LINK  "):
      if remark_flag == False:
        print_line += remark4 + "\n"
        remark_flag = True
      previous = current
      current = line[18:26]
      residue_name = line[17:20]
      remed = Remediator()
      chem_comp_atom_exch = remed.build_hash_from_chem_components(residue_name, remediated_out)
      #sys.stderr.write(str(chem_comp_atom_exch)+"\n")
      line = remediate_atomic_line(line, chem_comp_atom_exch)
    elif (remark_flag == False) and (remark_block == True): #deal with non-remark lines stuck in the top before main record types
      print_line += remark4 + "\n"
      remark_flag = True
    if previous == current:
      print_line += line + "\n"
    elif previous != current: # appears to check an entire residue for dna residue/atom names
      if re.search(r'^.{12}.\S..  .[ACTGIU]',print_line):
        if re.search(r'O2[\'|\*]   .',print_line) == None and previous != None:
          DNA_base = previous[1]
          if remediated_out == True:
            print_line = re.sub(r'(?m)(^.{12}.\S..)   '+DNA_base+' ',r'\g<1>  D'+DNA_base+' ',print_line)
            print_line = re.sub(r'(?m)(^TER.{15}) '+DNA_base+' ',r'\g<1>D'+DNA_base+' ',print_line)
          else:
            print_line = re.sub(r'(?m)(^.{12}.\S..)  D'+DNA_base+' ',r'\g<1>   '+DNA_base+' ',print_line)
            print_line = re.sub(r'(?m)(^TER.{15})D'+DNA_base+' ',r'\g<1> '+DNA_base+' ',print_line)

      if remediated_out == False:
        m = aa_re.search(print_line)
        if m:
          res = m.group(1)
          if re.search('1H   '+res,print_line) or re.search('2H   '+res,print_line):
            print_line = re.sub(' HN2 '+res,'2H   '+res,print_line)
      print_line=print_line.rstrip("\n")

      if not (print_line == ""):
        print(print_line, file=f)
      print_line = line + "\n"
  print_line=print_line.rstrip("\n")
  print(print_line, file=f)
  pdb_file.close()

#}}}

def remediator(params, log=None):
  if log == None:
    log = sys.stderr
  custom_dict = False
  remediated_out = True
  user_dict = ""
  file_name = params.file_name
  if params.version == "3.2":
    remediated_out = True
  elif params.version == "2.3":
    remediated_out = False
  if params.dict != None:
    custom_dict = True
    user_dict = params.dict
  if params.output_file != None:
    f = open(params.output_file, "w")
  else:
    f = sys.stdout
  atom_exch = build_hash(remediated_out,
                                  custom_dict,
                                  user_dict)
  if remediated_out:
    remediated_alt = False
  else:
    remediated_alt = True
  alt_atom_exch = build_hash(remediated_alt,
                                           custom_dict,
                                           user_dict)
  count, res_count = \
    pre_screen_file(file_name, atom_exch, alt_atom_exch)
  if count > 0 or res_count > 0:
    remediate(params.file_name, remediated_out, f)
    if params.output_file != None:
      f.close()
  else:
    if remediated_out:
      print("All atoms conform to PDB v3.x standard  **skipping**", file=log)
    else:
      print("All atoms conform to PDB v2.3 standard  **skipping**", file=log)
