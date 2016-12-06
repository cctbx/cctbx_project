def get_non_plugin_selection(plugin_object, #restraints_manager,
                             sites_cart,
                             hd_selection,
                             ignore_hd,
                             verbose=False,
                         ):
  if ignore_hd:
    general_selection = ~hd_selection
  else:
    general_selection = hd_selection|~hd_selection
  ligand_i_seqs = []
  print plugin_object
  print dir(plugin_object)
  print plugin_object.sites_cart_ptrs
  for ligand in plugin_object.sites_cart_ptrs:
    print ligand
    for group in ligand:
      print 'group',group
      ligand_i_seqs += group
  for i_seq in ligand_i_seqs:
    general_selection[i_seq] = False
  if verbose:
    print plugin_object
    print "\nNumber of atoms in selection : %d" % len(filter(None,
                                                             general_selection))
    print list(general_selection)
  return general_selection

def get_plugin_selection(plugin_object, #restraints_manager,
                         sites_cart,
                         hd_selection,
                         ignore_hd,
                         verbose=False,
                     ):
  if ignore_hd:
    hd_selection = ~hd_selection
  else:
    hd_selection = hd_selection|~hd_selection
  general_selection = hd_selection&~hd_selection
  ligand_i_seqs = []
  for ligand in plugin_object.sites_cart_ptrs:
    for group in ligand:
      ligand_i_seqs += group
  for i_seq in ligand_i_seqs:
    general_selection[i_seq] = True
  rc = general_selection&hd_selection
  if verbose:
    print plugin_object
    print "\nNumber of atoms in selection : %d" % len(filter(None,
                                                             general_selection))
  return rc

def _show_gradient(g):
  return "(%9.3f %9.3f %9.3f)" % (g)
