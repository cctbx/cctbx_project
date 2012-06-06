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

import sys
import os
import re
import libtbx.load_env

#{{{ get_summary
def get_summary():
  summary = """
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
  return summary
#}}}

#{{{ pre_screen_file
def pre_screen_file(filename, atom_exch, alt_atom_exch):
  count = 0
  pdb_file = open(filename)
  for line in pdb_file:
    line=line.rstrip()
    line=line.ljust(80)
    type_test = line[0:6]
    if type_test in ("ATOM  ", "HETATM", "TER   ", "ANISOU", "SIGATM", "SIGUIJ", "LINK  "):
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
      if atom_exch.has_key(clean_entry):
        if alt_atom_exch.has_key(clean_entry):
          pass
        else:
          count += 1
  return count
#}}}

#{{{ build_hash
#--Build Hash Table------------------------------------------------
def build_hash(remediated_out, custom_dict, user_dict):
  atom_exch = {}
  f = open(libtbx.env.find_in_repositories(
        relative_path="cctbx_project/iotbx/pdb/remediation/remediation.dict",
        test=os.path.isfile), "rb")
  if remediated_out == True: #converting to remediated
    for line in f:
      line=line.rstrip()
      new, old = line.split(':')
      atom_exch[old] = new
    remark4 = "REMARK   4 REMEDIATOR VALIDATED PDB VERSION 3.2 COMPLIANT".ljust(80)
  else: #converting to old
    for line in f:
      new, old = line.split(':')
      atom_exch[new] = old
    remark4 = "REMARK   4 REMEDIATOR VALIDATED PDB VERSION 2.3 COMPLIANT".ljust(80)
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
  return atom_exch, remark4
#}}}
#------------------------------------------------------------------

#{{{ remediate
#----PDB routine---------------------------------------------------
def remediate(filename, atom_exch, remediated_out, remark4, f=None):
  if f == None:
    f = sys.stdout
  previous = None
  current = None
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
      previous = current
      current = line[18:26]
      clean_entry = entry[0:4] + " " + entry[5:8]
      if atom_exch.has_key(clean_entry):
        line = line.replace(clean_entry[0:4],atom_exch[clean_entry][0:4],1)
    elif (remark_flag == False) and (remark_block == True): #deal with non-remark lines stuck in the top before main record types
      print_line += remark4 + "\n"
      remark_flag = True
    if previous == None:
      previous = current
    if previous == current:
      print_line += line + "\n"
    elif previous != current:
      if re.search(r'^.{12}.\S..  .[ACTGIU]',print_line):
        if re.search(r'O2[\'|\*]   .',print_line) == None:
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
      print >> f, print_line
      print_line = line + "\n"
  pdb_file.close()

  if re.search(r'^.{12}.\S..  .[ACTGIU]',print_line):
    if re.search(r'O2[\'|\*]   .',print_line) == None:
      DNA_base = previous[1]
      if remediated_out == True:
        print_line = re.sub(r'(?m)(^.{12}.\S..)   '+DNA_base,r'\g<1>  D'+DNA_base,print_line)
        print_line = re.sub(r'(?m)(^TER.{15}) '+DNA_base+' ',r'\g<1>D'+DNA_base+' ',print_line)
      else:
        print_line = re.sub(r'(?m)(^.{12}.\S..)  D'+DNA_base,r'\g<1>   '+DNA_base,print_line)
        print_line = re.sub(r'(?m)(^TER.{15})D'+DNA_base+' ',r'\g<1> '+DNA_base+' ',print_line)

    if remediated_out == False:
      m = aa_re.search(print_line)
      if m:
        res = m.group(1)
        if re.search('1H   '+res,print_line) or re.search('2H   '+res,print_line):
          print_line = re.sub(' HN2 '+res,'2H   '+res,print_line)

  print_line=print_line.rstrip("\n")
  print >> f, print_line
#}}}

def remediator(params, log=None):
  if log == None:
    log = sys.stderr
  custom_dict = False
  remedidated_out = True
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
  atom_exch, remark4 = build_hash(remediated_out,
                                  custom_dict,
                                  user_dict)
  if remediated_out:
    remediated_alt = False
  else:
    remediated_alt = True
  alt_atom_exch, remark4_dump = build_hash(remediated_alt,
                                           custom_dict,
                                           user_dict)
  count = pre_screen_file(file_name, atom_exch, alt_atom_exch)
  if count > 0:
    remediate(params.file_name, atom_exch, remediated_out, remark4, f)
    if params.output_file != None:
      f.close()
  else:
    if remediated_out:
      print >> log, "All atoms conform to PDB v3.x standard  **skipping**"
    else:
      print >> log, "All atoms conform to PDB v2.3 standard  **skipping**"
