"""BLAST search and fetch models"""

from __future__ import absolute_import, division, print_function

import mmtbx.model
import sys
import iotbx.pdb

import iotbx.bioinformatics.pdb_info
from iotbx.bioinformatics import local_blast
from iotbx.bioinformatics.structure import summarize_blast_output
from iotbx.pdb.fetch import fetch
import six

def get_best_homologues(model, chain_ids=None):
  """
  Do local BLAST search for homologues, pick the best resolution, fetch models
  chain_ids - list of chains to search. If none - do for all of them.
  As usual, multi-model is not supported.
  Returns dictionary of chain_id: model
  """

  # Load data
  pdb_info = iotbx.bioinformatics.pdb_info.pdb_info_local()

  h = model.get_hierarchy()
  result = {}
  for chain in h.only_model().chains():
    print("Working with chain '%s'" % chain.id)
    if not chain.is_protein():
      print("  skipping, not protein. Maybe water or NA")
      continue
    if chain_ids is not None and chain.id not in chain_ids():
      continue
    sequence = chain.as_padded_sequence()
    l_blast = local_blast.pdbaa(seq=sequence) # why sequence in init???
    blast_xml_result = l_blast.run()
    # Probably alignment info is lost in this call,
    # but we have tools to do alignment, see mmtbx/alignment/__init__.py
    blast_summary = summarize_blast_output("\n".join(blast_xml_result))
    pdb_ids_to_study = {}
    for hit in blast_summary:
      # print dir(hit)
      hit.show(out=sys.stdout)
      # Don't have clear idea what these values mean right now,
      # but it would be reasonable to filter somehow.
      if hit.identity < 70:
        continue
      pdb_ids_to_study[hit.pdb_id] = hit.chain_id # can add more info
    #
    info_list = pdb_info.get_info_list(list(pdb_ids_to_study.keys()))
    # It would be good to merge info_list with hits in blast_summary in
    # one data structure here and do better filtering.
    # Let's go without it for now.

    # Sort by resolution in place
    info_list.sort(key=lambda tup: tup[1])
    best_pdb_id = info_list[0][0]
    best_pdb_chain = pdb_ids_to_study[info_list[0][0]]
    print("Best pdb:", info_list[0], "chain:", pdb_ids_to_study[info_list[0][0]])
    # print info_list

    # Get actual selected model from PDB.
    data = fetch(id=best_pdb_id)
    m = mmtbx.model.manager(
        model_input=iotbx.pdb.input(source_info=None, lines=data.readlines()))
    sel = m.selection("chain '%s'" % best_pdb_chain)
    result[chain.id] = m.select(sel)
  return result

def run(args):
  # prepare dict. DO IT ONCE, THEN COMMENT OUT!!! This will fetch data from RCSB.
  # iotbx.bioinformatics.pdb_info.get_all_experimental_pdb_info_to_pkl()
  #------------------------


  model = mmtbx.model.manager(model_input=iotbx.pdb.input("1ucs.pdb"))
  r = get_best_homologues(model)
  print(r)
  for chain, model in six.iteritems(r):
    model.pdb_or_mmcif_string_info(
        target_format='pdb',
        target_filename="chain_%s.pdb" % chain,
        write_file=True)


if __name__ == '__main__':
  run(sys.argv[1:])
