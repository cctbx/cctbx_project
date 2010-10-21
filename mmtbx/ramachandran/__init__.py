
from __future__ import division
from scitbx import matrix
import libtbx.load_env
from libtbx.math_utils import ifloor, iceil
import libtbx.phil
from libtbx import adopt_init_args
import math
import sys
import os

import boost.python
ext = boost.python.import_ext("mmtbx_ramachandran_ext")
from mmtbx_ramachandran_ext import *

master_phil = libtbx.phil.parse("""
  rama_weight = 1.0
    .type = float
    .short_caption = Ramachandran gradients weight
    .expert_level = 1
  use_finite_differences = False
    .type = bool
    .help = Used for testing - not suitable for real structures.
    .short_caption = Use finite differences (DEVELOPERS ONLY)
    .expert_level = 3
""")

refine_opt_params = libtbx.phil.parse("""
#  min_allowed_d_min = 3.0
#    .type = float
#    .short_caption = Resolution cutoff for Ramachandran restraints
#    .expert_level = 2
  rama_selection = None
    .type = str
    .short_caption = Atom selection for Ramachandran restraints
    .style = selection
    .expert_level = 1
  exclude_secondary_structure = False
    .type = str
    .expert_level = 1
""")

def load_tables (params=None) :
  if (params is None) :
    params = master_phil.fetch().extract()
  from scitbx.array_family import flex
  tables = {}
  for residue_type in ["ala", "gly", "prepro", "pro"] :
    file_name = libtbx.env.find_in_repositories(
      relative_path="chem_data/rotarama_data/%s.rama.combined.data" %
        residue_type,
      test=os.path.isfile)
    f = open(file_name, "r")
    data = flex.double()
    for line in f.readlines() :
      val, phi, psi = line.split()
      assert ((int(phi) % 2 == 1) and (int(psi) % 2 == 1))
      data.append(float(val))
    t = lookup_table(data, 180)
    tables[residue_type] = t
  return tables

class generic_proxy (object) :
  restraint_type = None

class proxy (generic_proxy) :
  restraint_type = "ramachandran"
  def __init__ (self, i_seqs, residue_type) :
    assert (len(i_seqs) == 5)
    self.i_seqs = i_seqs
    self.residue_type = residue_type

class generic_restraints_helper (object) :
  def __init__ (self, params) :
    adopt_init_args(self, locals())
    self.tables = load_tables(params)

  def restraints_residual_sum (self,
                               sites_cart,
                               proxies,
                               gradient_array=None,
                               unit_cell=None) :
    ramachandran_proxies = []
    for proxy in proxies :
      if (proxy.restraint_type == "ramachandran") :
        ramachandran_proxies.append(proxy)
    return self._phi_psi_restraints_residual_sum(
      sites_cart=sites_cart,
      proxies=ramachandran_proxies,
      gradient_array=gradient_array)

  def _phi_psi_restraints_residual_sum (self,
                                        proxies,
                                        sites_cart,
                                        gradient_array=None) :
    if (gradient_array is None) :
      from scitbx.array_family import flex
      gradient_array = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
    sum = 0
    assert (self.params.rama_weight >= 0.0);
    for proxy in proxies :
      rama_table = self.tables[proxy.residue_type]
      if self.params.use_finite_differences :
        sum += rama_table.compute_gradients_finite_differences(
          gradient_array=gradient_array,
          sites_cart=sites_cart,
          i_seqs=proxy.i_seqs,
          weight=self.params.rama_weight,
          epsilon=0.001)
      else :
        sum += rama_table.compute_gradients(
          gradient_array=gradient_array,
          sites_cart=sites_cart,
          i_seqs=proxy.i_seqs,
          weight=self.params.rama_weight,
          epsilon=0.001)
    return sum

def extract_proxies (pdb_hierarchy,
                     atom_selection=None,
                     log=sys.stdout) :
  from iotbx.pdb.amino_acid_codes import one_letter_given_three_letter
  from cctbx import geometry_restraints
  from scitbx.array_family import flex
  if (atom_selection is None) :
    atom_selection = flex.bool(pdb_hierarchy.atoms().size(), True)
  proxies = []
  for model in pdb_hierarchy.models() :
    for chain in model.chains() :
      for conformer in chain.conformers() :
        residues = conformer.residues()
        for i, residue in enumerate(residues) :
          if (not residue.resname in one_letter_given_three_letter) :
            continue
          next_res, prev_res = None, None
          resseq2 = residue.resseq_as_int()
          resseq1, resseq3 = None, None
          if (i > 0):
            resseq1 = residues[i-1].resseq_as_int()
            if (resseq2 != (resseq1 + 1)) :
              continue
            prev_res = residues[i-1]
          if (i < (len(residues) - 1)) :
            resseq3 = residues[i+1].resseq_as_int()
            if (resseq2 != (resseq3 - 1)) :
              continue
            next_res = residues[i+1]
          if (next_res is not None) and (prev_res is not None) :
            c1, n2, ca2, c2, n3 = None, None, None, None, None
            for atom in prev_res.atoms() :
              if (atom.name == " C  ") :
                c1 = atom
                break
            for atom in residue.atoms() :
              if (atom.name == " N  ") :
                n2 = atom
              elif (atom.name == " CA ") :
                ca2 = atom
              elif (atom.name == " C  ") :
                c2 = atom
            for atom in next_res.atoms() :
              if (atom.name == " N  ") :
                n3 = atom
            if (None in [c1, n2, ca2, c2, n3]) :
              #print >> log, "  incomplete backbone for %s %d-%d, skipping." % \
              #  (chain.id, resseq1, resseq3)
              continue
            i_seqs = [c1.i_seq,n2.i_seq,ca2.i_seq,c2.i_seq,n3.i_seq]
            for i_seq in i_seqs :
              if (not atom_selection[i_seq]) :
                continue
            pep1 = geometry_restraints.bond(
              sites=[c1.xyz,n2.xyz],
              distance_ideal=1,
              weight=1)
            pep2 = geometry_restraints.bond(
              sites=[c2.xyz,n3.xyz],
              distance_ideal=1,
              weight=1)
            if (pep1.distance_model > 4) or (pep2.distance_model > 4) :
              continue
            if (residue.resname == "PRO") :
              residue_type = "pro"
            elif (residue.resname == "GLY") :
              residue_type = "gly"
            elif (residues[i+1].resname == "PRO") :
              residue_type = "prepro"
            else :
              residue_type = "ala"
            phi_psi = proxy(i_seqs, residue_type)
            proxies.append(phi_psi)
  print >> log, "%d Ramachandran restraints generated." % len(proxies)
  return proxies

def process_refinement_settings (
    params,
    pdb_hierarchy,
    secondary_structure_manager,
    d_min=None,
    log=sys.stdout,
    scope_name="refinement.ramachandran_restraints") :
  if (params.rama_selection is None) :
    atom_selection = flex.bool(pdb_hierarchy.atoms().size(), True)
  else :
    cache = pdb_hierarchy.atom_selection_cache()
    try :
      sele = cache.selection(params.rama_selection)
    except Exception, e :
      raise Sorry("""Atom selection error:
  %s
Selection string resulting in error:
  %s""" % (str(e), params.rama_selection))
    else :
      if (sele.count(True) == 0) :
        raise Sorry("""Empty atom selection for %s.rama_selection.
Current selection string:
  %s""" % (scope_name, params.rama_selection))
      atom_selection = sele
  if params.exclude_secondary_structure :
    alpha_sele = secondary_structure_manager.alpha_selection()
    beta_sele = secondary_structure_manager.beta_selection()
    ss_sele = alpha_sele | beta_sele
    atom_selection &= ~ss_sele
  return atom_selection
