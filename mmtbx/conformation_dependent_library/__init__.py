from __future__ import absolute_import, division, print_function
import sys
import copy
from string import digits
from string import ascii_letters as letters

import iotbx.pdb

from mmtbx.conformation_dependent_library.cdl_database import cdl_database
import mmtbx.conformation_dependent_library.cdl_utils

from mmtbx.conformation_dependent_library.multi_residue_class import \
  ThreeProteinResidues, RestraintsRegistry
from mmtbx.conformation_dependent_library.multi_residue_class import \
  TwoProteinResidues, FourProteinResidues, FiveProteinResidues
from mmtbx.conformation_dependent_library.multi_residue_cdl_class import \
  ThreeProteinResiduesWithCDL

from mmtbx.conformation_dependent_library.multi_base_class import \
  TwoNucleicResidues
from six.moves import range

chararcters_36 = letters[:26]+digits

registry = RestraintsRegistry()

def restraints_show(restraints_values):
  from mmtbx.conformation_dependent_library.cdl_setup import headers
  outl = ""
  for i, item in enumerate(restraints_values):
    if i%2==0:
      if i==0:
        s = "  %s, %s : %s %s\n"
      elif i<15:
        s = "  %-25s%s: %9.2f %9.2f\n"
      else:
        s = "  %-25s%s:   %9.4f %9.4f\n"
      outl += s % (headers[i],
                   headers[i+1],
                   restraints_values[i],
                   restraints_values[i+1],
        )
  return outl

def get_restraint_values(threes, interpolate=False):
  from mmtbx.conformation_dependent_library import utils
  res_type_group = cdl_utils.get_res_type_group(
    threes[1].resname,
    threes[2].resname,
  )
  if res_type_group is None: return None
  if interpolate:
    restraint_values = ["2", -1]
    key = threes.get_cdl_key(exact=interpolate)
    for i in range(2,26):
      grid = utils.get_grid_values(res_type_group, key[0], key[1], column=i)
      index = utils.get_index(*key)
      r = utils.interpolate_2d(grid, index)
      restraint_values.append(r)
  else:
    key = threes.get_cdl_key()
    if key is None: return None
    restraint_values = cdl_database[res_type_group][key]
  return restraint_values

def generate_residue_tuples(hierarchy,
                            geometry,
                            length,
                            dna_rna_residues=False,
                            include_non_linked=False,
                            backbone_only=True,
                            include_non_standard_residues=False,
                            # CDL specific
                            cdl_class=False,
                            omega_cdl=False,
                            #
                            retain_selection="name ca or name c or name n or name o or name cb or name h",
                            verbose=False,
                            ):
  assert length
  assert length>1
  if dna_rna_residues:
    assert length<=2
    LinkedResidues = TwoNucleicResidues
    residue_lookup = ['common_rna_dna']
    assert not include_non_standard_residues
  else:
    assert length<=5
    if length==3:
      if cdl_class:
        LinkedResidues = ThreeProteinResiduesWithCDL
      else:
        LinkedResidues = ThreeProteinResidues
    elif length==2: LinkedResidues = TwoProteinResidues
    elif length==4: LinkedResidues = FourProteinResidues
    elif length==5: LinkedResidues = FiveProteinResidues
    residue_lookup = ['common_amino_acid']
    if include_non_standard_residues:
      residue_lookup.append('modified_amino_acid')
  if backbone_only:
    backbone_asc = hierarchy.atom_selection_cache()
    backbone_sel = backbone_asc.selection(retain_selection)
    backbone_hierarchy = hierarchy.select(backbone_sel)
  get_class = iotbx.pdb.common_residue_names_get_class
  threes = LinkedResidues(geometry, registry=registry, length=length)
  loop_hierarchy=hierarchy
  if backbone_only: loop_hierarchy=backbone_hierarchy
  for model in loop_hierarchy.models():
    if verbose: print('model: "%s"' % model.id)
    for chain in model.chains():
      if verbose: print('chain: "%s"' % chain.id)
      for conformer in chain.conformers():
        if verbose: print('  conformer: altloc="%s"' % (
          conformer.altloc))
        while threes: del threes[0]
        threes.start=None
        threes.end=None
        list_of_threes = []
        for residue in conformer.residues():
          if verbose:
            if residue.resname not in ["HOH"]:
              print('    residue: resname="%s" resid="%s"' % (
                residue.resname, residue.resid()))
          if verbose: print('      residue class : %s' % get_class(residue.resname))
          if get_class(residue.resname) not in residue_lookup:
            continue
          if include_non_linked:
            list.append(threes, residue)
            if len(threes)>length: del threes[0]
          else:
            threes.append(residue)
          if len(threes)!=length:
            if omega_cdl:
              if len(threes)==length-1:
                threes.insert(0,None)
              else: continue
            else: continue
          assert len(threes)<=length
          #
          if 0: list_of_threes.append(copy.copy(threes))
          else:
            # transfer residues to new class because copy.copy is too deep
            tmp = LinkedResidues(geometry,
                                 registry=registry,
                                 length=length,
                                 include_non_linked=include_non_linked,
                                 )
            for pr in threes: tmp.append(pr)
            list_of_threes.append(tmp)
            #
        # per conformer
        for i, threes in enumerate(list_of_threes):
          if i==0:
            threes.start = True
          if i==len(list_of_threes)-1:
            threes.end = True
          else:
            if len(threes)!=length:
              pass
            elif threes[1] != list_of_threes[i+1][0]:
              threes.end = True
              list_of_threes[i+1].start = True
          yield threes
      threes = LinkedResidues(geometry,
                              registry=registry,
                              length=length,
                              include_non_linked=include_non_linked,
                              )

def generate_protein_tuples(hierarchy,
                            geometry,
                            length,
                            include_non_linked=False,
                            backbone_only=True,
                            include_non_standard_peptides=False,
                            # CDL specific
                            cdl_class=False,
                            omega_cdl=False,
                            #
                            retain_selection="name ca or name c or name n or name o or name cb or name h",
                            verbose=False,
                            ):
  for item in generate_residue_tuples(hierarchy,
                                      geometry,
                                      length,
                                      include_non_linked=include_non_linked,
                                      backbone_only=backbone_only,
                                      include_non_standard_residues=include_non_standard_peptides,
                                      # CDL specific
                                      cdl_class=cdl_class,
                                      omega_cdl=omega_cdl,
                                      #
                                      retain_selection=retain_selection,
                                      verbose=verbose,
                                      ):
    yield item

# retained for backwards compatibility
def generate_protein_threes(hierarchy,
                            geometry,
                            include_non_linked=False,
                            backbone_only=True,
                            include_non_standard_peptides=False,
                            # CDL specific
                            cdl_class=False,
                            omega_cdl=False,
                            verbose=False,
                            ):
  for threes in generate_protein_tuples(
    hierarchy,
    geometry,
    include_non_linked=include_non_linked,
    backbone_only=backbone_only,
    include_non_standard_peptides=include_non_standard_peptides,
    cdl_class=cdl_class,
    omega_cdl=omega_cdl,
    length=3,
    ):
    yield threes

def generate_protein_fragments(hierarchy,
                               geometry,
                               length,
                               include_non_linked=False,
                               backbone_only=True,
                               include_non_standard_peptides=False,
                               # include_non_protein_linked=False, # NH2 1KYC
                               verbose=False,
                               ):
  for fragment in generate_residue_tuples(
    hierarchy,
    geometry,
    length,
    include_non_linked=include_non_linked,
    backbone_only=backbone_only,
    include_non_standard_residues=include_non_standard_peptides,
    # include_non_protein_linked=include_non_protein_linked,
    verbose=verbose,
    ):
    yield fragment

def generate_dna_rna_fragments(hierarchy,
                               geometry,
                               length,
                               include_non_linked=False,
                               backbone_only=False,
                               include_non_standard_bases=False,
                               retain_selection='all',
                               verbose=False,
                               ):
  assert not backbone_only, 'backbone_only not available with DNA/RNA'
  for item in generate_residue_tuples(hierarchy,
                                      geometry,
                                      length,
                                      dna_rna_residues=True,
                                      include_non_linked=include_non_linked,
                                      backbone_only=False,
                                      include_non_standard_residues=include_non_standard_bases,
                                      #retain_selection=retain_selection,
                                      verbose=verbose,
                                      ):
    yield item

def update_restraints(hierarchy,
                      geometry, # restraints_manager,
                      current_geometry=None, # xray_structure!!
                      sites_cart=None,
                      cdl_proxies=None,
                      ideal=True,
                      esd=True,
                      esd_factor=1.,
                      interpolate=False,
                      log=None,
                      verbose=False,
                      ):
  global registry
  registry = RestraintsRegistry()
  if current_geometry:
    assert not sites_cart
    sites_cart = current_geometry.sites_cart()
  if sites_cart:
    pdb_atoms = hierarchy.atoms()
    # XXX PDB_TRANSITION VERY SLOW
    for j_seq, atom in enumerate(pdb_atoms):
      atom.xyz = sites_cart[j_seq]

  threes = None
  average_updates = 0
  total_updates = 0
  for threes in generate_protein_threes(hierarchy,
                                        geometry, #restraints_manager,
                                        cdl_class=True,
                                        #verbose=verbose,
                                        ):
    if threes.cis_group():
      continue

    restraint_values = get_restraint_values(threes, interpolate=interpolate)

    if restraint_values is None: continue

    if restraint_values[0]=="I":
      average_updates += 1
    else:
      total_updates += 1
    threes.apply_updates(restraint_values,
                         cdl_proxies,
                         ideal=ideal,
                         esd=esd,
                         esd_factor=esd_factor,
                         )
  if registry.n: threes.apply_average_updates(registry)
#  restraints_manager.geometry.reset_internals()
  geometry.reset_internals()
  if verbose and threes and threes.errors:
    if log:
      log.write("  Residues not completely updated with CDL restraints\n\n")
    for line in threes.errors:
      if log:
        log.write("%s\n" % line)
      else:
        print(line)
#  print 'average updates',average_updates,total_updates
#  assert average_updates==0
  return geometry #restraints_manager

def run(filename):
  if False:
    for i in range(-188,188):
      print(i,round_to_ten(i),abs(i-round_to_ten(i)))
    assert 0

  pdb_inp = iotbx.pdb.input(filename)
  hierarchy = pdb_inp.construct_hierarchy()
  hierarchy.atoms().reset_serial()
  assert 0, "broken run method"
  update_restraints(hierarchy,
                    #verbose=True,
                    )

if __name__=="__main__":
  if 0:
    psi = -180
    lookup = "Gly_nonxpro"
    print(lookup)
    for phi in range(170,181):
      key = (round_to_ten(psi),round_to_ten(phi))
      print('key',psi,phi,round_to_ten(psi),round_to_ten(phi),key, end=' ')
      print(cdl_database[lookup][key][:4])

  run(sys.argv[1])
