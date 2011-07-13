# (jEdit options) :folding=explicit:collapseFolds=1:
#This module contains tools for annotating protein structures with probabilities
#  of secondry structure classification. This is the central feature of the
#  cablam system.
#In a nutshell, cablam_training is used to create data for each dssp code. These
#  data are binned, and the probability of point membership in each
#  conformational bin is calculated. Bayesian statistics are then used to
#  calculate the likelihood that a point would be given a particular dssp code
#  (see cablam_datamanagement.py). The probabilities are stored as NDimTable
#  contours accessible though mmtbx.rotamer.n_dim_table . There may be an
#  external help file or a forthcoming publication which will explain in greater
#  detail.
#This module is the most likely to be under heavy development. Function names
#  may be subject to change.

from mmtbx.rotamer.n_dim_table import NDimTable
from libtbx import easy_pickle
#from libtbx.utils import Sorry
import libtbx.load_env
#import weakref

import sys, os
#from mmtbx.cablam import cablam_res, cablam_math

def annote_dssp_3d(resdata, dsspcode):
  picklefile = libtbx.env.find_in_repositories(
    relative_path=("chem_data/cablam_data/general."+dsspcode+".calc3d.pickle"),
    test=os.path.isfile)
  if (picklefile is None):
    sys.stderr.write("""
      Could not find a needed pickle file for dsspcode "+dsspcode+" in chem_data
      Exiting.
      """)
    sys.exit()
  ndt = easy_pickle.load(file_name=picklefile)
  for resid in resdata:
    residue = resdata[resid]
    try:
      CA_d_in, CA_d_out = residue.results['CA_d_in'],residue.results['CA_d_out']
      CA_a = residue.results['CA_a']
    except KeyError:
      continue
    residue.dssprob[dsspcode] = ndt.valueAt([CA_d_in,CA_d_out,CA_a])

def print_annoted_text_human(resdata, writeto=sys.stdout):
#Note: "human" vs "machine" is not very descriptive, and should be changed
  reskeys = resdata.keys()
  reskeys.sort()
  writeto.write("pdb:chainID:resnum:1st:2nd:3rd:4th:5th:6th:7th:8th\n")
  for resid in reskeys:
    residue = resdata[resid]
    if not residue.dssprob:
      continue
    residue.dssprob_sort()
#Do I want to cull out anything without predictions to make text parsing easier?
    outlist = [residue.pdbid,"%2s"%residue.chain,"%4s"%residue.resnum]
    for dsspcode in residue.dssprob_decending:
      outlist.append(dsspcode + ' ' + "%.4f"%residue.dssprob[dsspcode])
    outline = ':'.join(outlist)
    writeto.write(outline)
    writeto.write("\n")

def print_annoted_text_machine(resdata, writeto=sys.stdout):
  reskeys = resdata.keys()
  reskeys.sort()
  writeto.write("pdb:chainID:resnum:H:G:E:S:T:X:B:I\n")
  for resid in reskeys:
    residue = resdata[resid]
    if not residue.dssprob:
      continue
#Do I want to cull out anything without predictions to make text parsing easier?
    outlist = [residue.pdbid,"%2s"%residue.chain,"%4s"%residue.resnum]
    for dsspcode in ['H','G','E','S','T','X','B','I']:
      outlist.append("%.4f"%residue.dssprob[dsspcode])
    outline = ':'.join(outlist)
    writeto.write(outline)
    writeto.write("\n")

#{{{ print_annoted_kin function
#This function outputs a kinemage-format annotation of secondary structure
#  probabilities as calculated by cablam.
#-------------------------------------------------------------------------------
def print_annoted_kin(resdata, writeto=sys.stdout):
  reskeys = resdata.keys()
  reskeys.sort()
  writeto.write("""
@pointmaster 'h' {H (red)}
@pointmaster 'g' {G (purple)}
@pointmaster 'i' {I (hotpink)}
@pointmaster 'e' {E (green)}
@pointmaster 'b' {B (blue)}
@pointmaster 's' {S (yellow)}
@pointmaster 't' {T (orange)}
@pointmaster 'x' {X (white)}
@group {dssprob stacks}
@ringlist {100% rings} color=gray radius=1.0
""")
  for resid in reskeys:
    residue = resdata[resid]
    firstalt = residue.firstalt('CA')
    try:
      CAxyz = residue.atomxyz[firstalt]['CA']
    except KeyError:
      continue
    writeto.write(
      "{ 100% ring } "+str(CAxyz[0])+" "+str(CAxyz[1])+" "+str(CAxyz[2])+"\n")
  writeto.write("""
@balllist {balls} nohighlight
""")
  for resid in reskeys:
    residue = resdata[resid]
    firstalt = residue.firstalt('CA')
    try:
      CAxyz = residue.atomxyz[firstalt]['CA']
    except KeyError:
      continue
    residue.dssprob_sort()
    printorder = residue.dssprob_decending
    #printorder = residue.dssprob.keys()  #list of dsspkeys
    #printorder.sort(reverse=True, key=lambda x: residue.dssprob[x])
    for order in printorder:
      if residue.dssprob[order] < 0.0001:
        continue #skip display for very low dssprobabilities
      if   order == "H": color="red"
      elif order == "G": color="purple"
      elif order == "I": color="hotpink"
      elif order == "E": color="green"
      elif order == "B": color="blue"
      elif order == "S": color="yellow"
      elif order == "T": color="orange"
      else: color="white" # case X
      writeto.write(
        "{"+order+" %.4f"%residue.dssprob[order]+"} "+str(CAxyz[0])+" "
        +str(CAxyz[1])+" "+str(CAxyz[2])+ " r=%.4f"%residue.dssprob[order]+" "
        +color+" \'"+order.lower()+"\'\n")
#-------------------------------------------------------------------------------
#}}}

#{{{ output_to_kin_all (kinemage output for)
def output_to_kin_all(resdata, writeto=sys.stdout):
  reskeys = resdata.keys()
  reskeys.sort()
  writeto.write("\n\n@group {H} animate\n")
  output_to_kin_dssp(resdata,reskeys,"H", probrange=(0.10,0.25), writeto=writeto, color="pinktint")
  output_to_kin_dssp(resdata,reskeys,"H", probrange=(0.25,0.50), writeto=writeto, color="pink")
  output_to_kin_dssp(resdata,reskeys,"H", probrange=(0.50,1.00), writeto=writeto, color="red")
  writeto.write("\n@group {G} animate\n")
  output_to_kin_dssp(resdata,reskeys,"G", probrange=(0.10,0.25), writeto=writeto, color="lilactint")
  output_to_kin_dssp(resdata,reskeys,"G", probrange=(0.25,0.50), writeto=writeto, color="lilac")
  output_to_kin_dssp(resdata,reskeys,"G", probrange=(0.50,1.00), writeto=writeto, color="purple")
  writeto.write("\n@group {I} animate\n")
  output_to_kin_dssp(resdata,reskeys,"I", probrange=(0.10,0.25), writeto=writeto, color="hotpink")
  output_to_kin_dssp(resdata,reskeys,"I", probrange=(0.25,0.50), writeto=writeto, color="hotpink")
  output_to_kin_dssp(resdata,reskeys,"I", probrange=(0.50,1.00), writeto=writeto, color="hotpink")
  writeto.write("\n@group {E} animate\n")
  output_to_kin_dssp(resdata,reskeys,"E", probrange=(0.10,0.25), writeto=writeto, color="greentint")
  output_to_kin_dssp(resdata,reskeys,"E", probrange=(0.25,0.50), writeto=writeto, color="sea")
  output_to_kin_dssp(resdata,reskeys,"E", probrange=(0.50,1.00), writeto=writeto, color="green")
  writeto.write("\n@group {B} animate\n")
  output_to_kin_dssp(resdata,reskeys,"B", probrange=(0.10,0.25), writeto=writeto, color="bluetint")
  output_to_kin_dssp(resdata,reskeys,"B", probrange=(0.25,0.50), writeto=writeto, color="sky")
  output_to_kin_dssp(resdata,reskeys,"B", probrange=(0.50,1.00), writeto=writeto, color="blue")
  writeto.write("\n@group {S} animate\n")
  output_to_kin_dssp(resdata,reskeys,"S", probrange=(0.10,0.25), writeto=writeto, color="yellowtint")
  output_to_kin_dssp(resdata,reskeys,"S", probrange=(0.25,0.50), writeto=writeto, color="yellow")
  output_to_kin_dssp(resdata,reskeys,"S", probrange=(0.50,1.00), writeto=writeto, color="gold")
  writeto.write("\n@group {T} animate\n")
  output_to_kin_dssp(resdata,reskeys,"T", probrange=(0.10,0.25), writeto=writeto, color="peachtint")
  output_to_kin_dssp(resdata,reskeys,"T", probrange=(0.25,0.50), writeto=writeto, color="peach")
  output_to_kin_dssp(resdata,reskeys,"T", probrange=(0.50,1.00), writeto=writeto, color="orange")
  writeto.write("\n@group {X} animate\n")
  output_to_kin_dssp(resdata,reskeys,"X", probrange=(0.10,0.25), writeto=writeto, color="brown")
  output_to_kin_dssp(resdata,reskeys,"X", probrange=(0.25,0.50), writeto=writeto, color="gray")
  output_to_kin_dssp(resdata,reskeys,"X", probrange=(0.50,1.00), writeto=writeto, color="white")
#}}}

#{{{ output_to_kin_dssp  (kinemage output for a single dsspcode, single range)
def output_to_kin_dssp(resdata, reskeys, dsspcode, probrange=(0.0,0.0),
  color="white", writeto=sys.stdout):
  #This function by itself lacks a .write("@group") and may make untidy .kins
  balllistline= "\n@balllist {"+ str(probrange[0]) +" to "+ str(probrange[1])+"} radius= 0.45 color="+ color+"\n"
  writeto.write(balllistline)
  for resid in reskeys:
    residue = resdata[resid]
    try:
      dssprob = residue.dssprob[dsspcode]
    except KeyError:
      continue
    if dssprob < probrange[0] or dssprob >= probrange[1]:
      continue
    else:
      firstalt = residue.firstalt('CA')
      CAxyz = residue.atomxyz[firstalt]['CA']
      #print "{",residue.alts(firstalt)["resname"],residue.resnum,"}",CAxyz[0],CAxyz[1],CAxyz[2]
      kinline = "{"+dsspcode+str(dssprob)+"} "+ " ".join(CAxyz[0],CAxyz[1],CAxyz[2]) +"\n"
      writeto.write(kinline)
      #print "{",dsspcode,dssprob,"}",CAxyz[0],CAxyz[1],CAxyz[2]
#}}}
