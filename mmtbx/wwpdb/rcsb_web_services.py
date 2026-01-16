
"""
Module for querying the RCSB web server using the REST API, as described here:
https://search.rcsb.org/index.html#search-api
"""

from __future__ import absolute_import, division, print_function
import libtbx.utils
import json
import requests

search_base_url = "https://search.rcsb.org/rcsbsearch/v2/query?json="
report_base_url = "https://data.rcsb.org/graphql"

def value_attribute_filter(attribute_name, operator, value):
  assert operator in ["greater", "less", "less_or_equal", "greater_or_equal", "exact_match"]
  filt = {
    "type": "terminal",
    "service": "text",
    "parameters": {
      "attribute": attribute_name,
      "operator": operator,
      "value" : value
    }
  }
  return filt

def data_only_filter():
  return value_attribute_filter(
      "rcsb_accession_info.has_released_experimental_data", "exact_match", "Y")

def xray_only_filter():
  return value_attribute_filter(
      "exptl.method", "exact_match", "X-RAY DIFFRACTION")

def polymeric_type_filter(value="Protein (only)"):
  assert value in ["Protein (only)", "Protein/NA", "Nucleic acid (only)", "Other"]
  return value_attribute_filter(
      "rcsb_entry_info.selected_polymer_entity_types", "exact_match", value)

def clashscore_filter(operator, value):
  return value_attribute_filter(
      "pdbx_vrpt_summary_geometry.clashscore", operator, value)

def rama_outliers_filter(operator, value):
  return value_attribute_filter(
      "pdbx_vrpt_summary_geometry.percent_ramachandran_outliers", operator, value)

def rota_outliers_filter(operator, value):
  return value_attribute_filter(
      "pdbx_vrpt_summary_geometry.percent_rotamer_outliers", operator, value)

def resolution_filter(operator, value):
  return value_attribute_filter(
      "rcsb_entry_info.diffrn_resolution_high.value", operator, value)

sort_by_res = \
    {
      "sort_by": "rcsb_entry_info.resolution_combined",
      "direction": "asc"
    }

def add_nodes_to_query_if_needed_in_place(query_json):
  if "nodes" not in query_json["query"].keys():
    query_json["query"]["type"] = "group"
    query_json["query"]["logical_operator"] = "and"
    del query_json["query"]["service"]
    query_json["query"]["nodes"] = []

def post_query(query_json=None, xray_only=True, d_max=None, d_min=None,
    protein_only=False, data_only=False, log=None,
    sort_by_resolution=False, clashscore_range=None,
    rama_outliers_range=None, rota_outliers_range=None):
  """  Make request to RCSB search API and return list of PDB ids, optionally with
  chain IDs. If query_json is not supplied, generic one will be used which
  searches for everything in PDB. It will be enhanced according to other parameters.


  Args:
      query_json (dict, optional): _description_. Defaults to None.
      xray_only (bool, optional): Return only xray structures. Defaults to True.
      d_max (_type_, optional): Max resolution. Defaults to None.
      d_min (_type_, optional): Min resolution. Defaults to None.
      protein_only (bool, optional): Return only protein entries. Defaults to False.
      data_only (bool, optional): Return only entries with experimental data. Defaults to False.
      log (_type_, optional): Handler for log. Defaults to None.
      sort_by_resolution (bool, optional): Sort by entry resolution. Defaults to False.
      clashscore_range (tuple): tuple of min and max clashscore, e.g. (0,10) or None

  Returns:
      list: PDB ids
  """
  if query_json is None:
    query_json = {
        "query": {
          "type": "terminal",
          "service": "text"
        },
        "return_type": "entry",
        "request_options": {
          "return_all_hits": True,
        }
      }
  if d_max is not None and d_min is not None:
    assert d_max > d_min

  if (log is None):
    log = libtbx.utils.null_out()
  print("Setting up RCSB server query:", file=log)
  if (xray_only):
    print("  limiting to X-ray structures", file=log)
    add_nodes_to_query_if_needed_in_place(query_json)
    query_json["query"]["nodes"].append(xray_only_filter())
    if (data_only):
      add_nodes_to_query_if_needed_in_place(query_json)
      query_json["query"]["nodes"].append(data_only_filter())
  if d_max is not None:
    add_nodes_to_query_if_needed_in_place(query_json)
    query_json["query"]["nodes"].append(resolution_filter("less", d_max))
  if d_min is not None:
    add_nodes_to_query_if_needed_in_place(query_json)
    query_json["query"]["nodes"].append(resolution_filter("greater", d_min))
  if clashscore_range is not None:
    add_nodes_to_query_if_needed_in_place(query_json)
    query_json["query"]["nodes"].append(clashscore_filter("greater", clashscore_range[0]))
    query_json["query"]["nodes"].append(clashscore_filter("less", clashscore_range[1]))
  if rama_outliers_range is not None:
    add_nodes_to_query_if_needed_in_place(query_json)
    query_json["query"]["nodes"].append(rama_outliers_filter("greater", rama_outliers_range[0]))
    query_json["query"]["nodes"].append(rama_outliers_filter("less", rama_outliers_range[1]))
  if rota_outliers_range is not None:
    add_nodes_to_query_if_needed_in_place(query_json)
    query_json["query"]["nodes"].append(rota_outliers_filter("greater", rota_outliers_range[0]))
    query_json["query"]["nodes"].append(rota_outliers_filter("less", rota_outliers_range[1]))
  if (protein_only):
    add_nodes_to_query_if_needed_in_place(query_json)
    query_json["query"]["nodes"].append(polymeric_type_filter("Protein (only)"))
  if (sort_by_resolution):
    if "sort" not in query_json["request_options"].keys():
      query_json["request_options"]["sort"] = []
    query_json["request_options"]["sort"].append(sort_by_res)
    print("  will sort by resolution", file=log)
  if "results_verbosity" not in query_json["request_options"].keys():
    query_json["request_options"]["results_verbosity"] = "compact"
  print("  executing HTTP request...", file=log)
  # print(json.dumps(query_json, indent=4))
  r = requests.post(search_base_url, json=query_json)
  res_ids = []
  # print('r.status_code', r.status_code)
  if r.status_code == 200:
    r_json = r.json()
    # print(json.dumps(r_json, indent=4))
    res_ids = r_json["result_set"]
  return res_ids

def sequence_search(
    sequence,
    identity_cutoff=90,
    target="pdb_protein_sequence",
    e_value_cutoff=1000000,
    **kwds):
  """
  Homology search for an amino acid sequence.  The advantage of using this
  service over the NCBI/EBI BLAST servers (in iotbx.pdb.fetch) is the ability
  to exclude non-Xray structures.

  identity_cutoff : int, optional
      Apply sequence filtering to remove structures with greater than this
      percentage sequence identity.
  e_value_cutoff : float, optional
      Hits with an E-Value above the cutoff value are filtered out
  """
  sequence_query = """
{
  "query": {
    "type": "group",
    "logical_operator": "and",
    "nodes": [
      {
        "type": "terminal",
        "service": "sequence",
        "parameters": {
          "evalue_cutoff": %s,
          "identity_cutoff": %s,
          "target": "%s",
          "value": "%s"
        }
      }
    ]
  },
  "return_type": "polymer_entity",
  "request_options": {
    "return_all_hits": true,
    "scoring_strategy": "sequence",
    "sort": [
      {
        "sort_by": "score",
        "direction": "desc"
      }
    ],
    "results_verbosity":"minimal"
  }
}"""
  assert target in ["pdb_protein_sequence", "pdb_dna_sequence", "pdb_rna_sequence"]
  assert 0 < identity_cutoff < 100
  sqr = sequence_query % (e_value_cutoff, identity_cutoff/100, target, sequence)
  jsq = json.loads(sqr)
  return post_query(query_json=jsq, **kwds)


def reference_chain_search(sequence, identity_cutoff=0.9, include_xray=True, include_csm=False, **kwds):
  """ Searches sequence optionally include computed models,
  returns pdb_id with chain id that matches.

  Args:
      sequence (str): _description_
      identity_cutoff (float, optional): _description_. Defaults to 0.9.
  """
  model_choice = ""
  if include_xray:
    model_choice = '"experimental"'
  if include_csm and include_xray:
    model_choice += ', "computational"'
  if include_csm and not include_xray:
    model_choice = '"computational"'

  query= """
{
  "query": {
    "type": "terminal",
    "service": "sequence",
    "parameters": {
      "evalue_cutoff": 0.1,
      "identity_cutoff": %s,
      "sequence_type": "protein",
      "value": "%s"
    }
  },
  "return_type": "polymer_instance",
  "request_options": {
    "return_all_hits": true,
    "results_content_type": [ %s ],
    "scoring_strategy": "combined",
    "sort": [
      {
        "sort_by": "reflns.d_resolution_high",
        "direction": "asc"
      }
    ]
  }
}
"""
  sqr = query % (identity_cutoff, sequence, model_choice)
  # print(sqr)
  jsq = json.loads(sqr)
  return post_query(query_json=jsq, xray_only=False, **kwds)


def chemical_id_search(resname, **kwds):
  """
  Find all entry IDs with the specified chemical ID.

  Parameters
  ----------
  resname : str
  kwds : dict
      Keyword arguments passed to post_query.

  Examples
  --------
  >>> from mmtbx.wwpdb.rcsb_web_services import chemical_id_search
  >>> len(chemical_id_search("ZN", data_only=True))
  2874
  """
  chem_comp_query = """
{
  "query": {
    "type": "group",
    "nodes": [
      {
        "type": "group",
        "logical_operator": "or",
        "nodes": [
          {
            "type": "terminal",
            "service": "text",
            "parameters": {
              "attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
              "operator": "exact_match",
              "value": "%s"
            }
          },
          {
            "type": "terminal",
            "service": "text",
            "parameters": {
              "attribute": "rcsb_polymer_entity_container_identifiers.chem_comp_monomers",
              "operator": "exact_match",
              "value": "%s"
            }
          }
        ]
      }
    ],
    "logical_operator": "and"
  },
  "return_type": "entry",
  "request_options": {
    "return_all_hits": true,
    "scoring_strategy": "combined",
    "sort": [
      {
        "sort_by": "score",
        "direction": "desc"
      }
    ]
  }
}
"""
  assert (1 <= len(resname) <= 3) or (len(resname)==5)
  sqr = chem_comp_query % (resname, resname)
  jsq = json.loads(sqr)
  return post_query(query_json=jsq, **kwds)

def get_high_resolution_for_structures(pdb_ids):
  with_res_count = get_high_resolution_and_residue_count_for_structures(pdb_ids)
  result = [r[:2] for r in with_res_count]
  return result

def post_report_query_with_pdb_list(query, pdb_ids):
  pdb_list = "%s" % pdb_ids
  pdb_list = pdb_list.replace("'", '"')
  request = query.format(pdb_list=pdb_list)
  r = requests.post(report_base_url, json={"query":request})
  return r.json()

def get_r_work_rfree_for_structures(pdb_ids):
  """ Get Rwork and Rfree for list of pdb ids
  Args:
      pdb_ids (list): pdb ids
  Returns:
      list: [[pdb_id, rwork, rfree], [pdb_id, rwork, rfree], ...]
  """
  query = """
  {{
    entries(entry_ids: {pdb_list} )
    {{
      rcsb_id
      refine
      {{
        ls_R_factor_R_free,
        ls_R_factor_R_work,
      }}
    }}
  }}"""
  r_json = post_report_query_with_pdb_list(query, pdb_ids)
  result = []
  for res in r_json["data"]["entries"]:
    pdb_id = str(res["rcsb_id"])
    rwork = None
    rfree = None
    if res["refine"] is not None:
      rwork = res["refine"][0]["ls_R_factor_R_work"]
      if rwork is not None:
        rwork = float(rwork)
      rfree = res["refine"][0]["ls_R_factor_R_free"]
      if rfree is not None:
        rfree = float(rfree)
    result.append([pdb_id, rwork, rfree])
  return result

def get_high_resolution_and_residue_count_for_structures(pdb_ids):
  query = """
  {{
    entries(entry_ids: {pdb_list} )
    {{
      rcsb_id
      rcsb_entry_info {{
        deposited_polymer_monomer_count
      }}
      refine {{
        ls_d_res_high
      }}
    }}
  }}"""

  r_json = post_report_query_with_pdb_list(query, pdb_ids)
  result = []
  for res in r_json["data"]["entries"]:
    pdb_id = str(res["rcsb_id"])
    resol = None
    n_res = None
    if res["refine"] is not None:
      res_resol = res["refine"][0]["ls_d_res_high"]
      resol = None if res_resol is None else float(res_resol)
    if res["rcsb_entry_info"] is not None:
      res_n_res = res["rcsb_entry_info"]["deposited_polymer_monomer_count"]
      n_res = None if res_n_res is None else int(res_n_res)
    result.append([pdb_id, resol, n_res])
  return result

def get_ligand_info_for_structures(pdb_ids):
  """Return a list of ligands in the specified structure(s), including the
  SMILES strings.

  Returns list of lists with
  [PDB ID, chain id, lig id, lig MW, lig Formula, lig name, lig SMILES]
  If the same ligand is present in different chains, it will produce several
  entries which will be different in chain ids.

  [[u'1MRU', u'A', u'MG', u'24.31', u'Mg 2', u'MAGNESIUM ION', u'[Mg+2]'],
   [u'1MRU', u'B', u'MG', u'24.31', u'Mg 2', u'MAGNESIUM ION', u'[Mg+2]']]

  """
  query = """
  {{
    entries(entry_ids: {pdb_list} )
    {{
      rcsb_id
      nonpolymer_entities {{
        nonpolymer_comp {{
          chem_comp {{
            formula
            formula_weight
            id
            name
          }}
          rcsb_chem_comp_descriptor {{
            SMILES
          }}
        }}
        rcsb_nonpolymer_entity_container_identifiers {{
          auth_asym_ids
        }}
      }}
    }}
  }}"""

  r_json = post_report_query_with_pdb_list(query, pdb_ids)
  result = []
  for res in r_json["data"]["entries"]:
    pdb_id = str(res["rcsb_id"])
    for npe in res["nonpolymer_entities"]:
      smiles = str(npe["nonpolymer_comp"]["rcsb_chem_comp_descriptor"]["SMILES"])
      lig_id = str(npe["nonpolymer_comp"]["chem_comp"]["id"])
      lig_mw = float(npe["nonpolymer_comp"]["chem_comp"]["formula_weight"])
      lig_formula = str(npe["nonpolymer_comp"]["chem_comp"]["formula"])
      lig_name = str(npe["nonpolymer_comp"]["chem_comp"]["name"])
      for chain_id in npe["rcsb_nonpolymer_entity_container_identifiers"]["auth_asym_ids"]:
        c_id = str(chain_id)
        result.append([pdb_id, c_id, lig_id, lig_mw, lig_formula, lig_name, smiles])
  return result

def get_emdb_id_for_pdb_id(pdb_id):
  """ Find out EMDB ID given PDB ID by quering RCSB portal.

  Args:
      pdb_id (str): pdb id
  Returns:
    list of emdb ids, e.g. ['EMD-37438'] or None if X-ray or not defined
  """

  graphql_query = '''
query
{
  entry(entry_id:"%s") {
    exptl {
      method
    }
    rcsb_entry_container_identifiers {
      emdb_ids
    }
  }
}
''' % pdb_id
  r = requests.post(report_base_url, json={"query":graphql_query})
  data_entry = r.json()['data']['entry']
  if not data_entry:
    return None
  if data_entry['exptl'][0]['method'] != 'ELECTRON MICROSCOPY':
    return None
  emdb_ids = data_entry['rcsb_entry_container_identifiers']['emdb_ids']
  if len(emdb_ids)==0:
    return None
  return emdb_ids

def get_similar_ligands_via_smiles(smiles, match_type='sub-struct-graph-relaxed-stereo', **kwds):
  # graph-relaxed-stereo
  # graph-relaxed
  # fingerprint-similarity
  # sub-struct-graph-relaxed-stereo
  # sub-struct-graph-relaxed
  similar_ligand_query = '''
{
  "query": {
    "type": "terminal",
    "service": "chemical",
    "parameters": {
      "type": "descriptor",
      "value": "%s",
      "descriptor_type": "SMILES",
      "match_type": "%s"
    }
  },
  "return_type": "mol_definition",
  "request_options": {
    "results_verbosity":"minimal",
    "paginate": {
      "start": 0,
      "rows": 100
    },
    "sort": [
      {
        "sort_by": "score",
        "direction": "desc"
      }
    ],
    "scoring_strategy": "combined"
  }
}
'''
  assert (3 <= len(smiles)), 'short SMILES "%s" return too many results' % smiles
  sqr = similar_ligand_query % (smiles, match_type)
  jsq = json.loads(sqr)
  return post_query(query_json=jsq, xray_only=False, **kwds)
