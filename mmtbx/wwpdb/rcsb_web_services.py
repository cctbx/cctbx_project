
"""
Module for querying the RCSB web server using the REST API, as described here:
http://www.rcsb.org/pdb/software/rest.do

There is some overlap with iotbx.pdb.fetch, which really should have gone here
instead, but this module is intended to be used in higher-level automation
pipelines.
"""

from __future__ import division
from __future__ import print_function
import libtbx.utils
from xml.dom.minidom import parseString
import sys

url_base = "http://www.rcsb.org/pdb/rest"
url_search = url_base + "/search"

def post_query(query_xml, xray_only=True, d_max=None, d_min=None,
    protein_only=False, data_only=False, identity_cutoff=None, log=None,
    sort_by_resolution=False):
  """
  Generate the full XML for a multi-part query with generic search options,
  starting from the basic query passed by another function, post it to the
  RCSB's web service, and return a list of matching PDB IDs.

  Parameters
  ----------
  query_xml : str or None
  xray_only : bool, optional
  d_max : float, optional
  d_min : float, optional
  protein_only : bool, optional
  data_only : bool, optional
  identity_cutoff : int, optional
      Apply sequence filtering to remove structures with greater than this
      percentage sequence identity. Must be one of 100, 95, 90, 70, 50, 40, 30.
  log : file, optional
  sort_by_resolution : bool, optional
  """
  if (log is None):
    log = libtbx.utils.null_out()
  queries = []
  if (query_xml is not None):
    queries.append(query_xml)
  print("Setting up RCSB server query:", file=log)
  if (xray_only):
    print("  limiting to X-ray structures", file=log)
    xray_query = "<queryType>org.pdb.query.simple.ExpTypeQuery</queryType>\n"+\
      "<mvStructure.expMethod.value>X-RAY</mvStructure.expMethod.value>"
    if (data_only):
      xray_query += "\n<mvStructure.hasExperimentalData.value>Y</mvStructure.hasExperimentalData.value>"
    queries.append(xray_query)
  if (d_max is not None) or (d_min is not None):
    base_clause = "<queryType>org.pdb.query.simple.ResolutionQuery</queryType>"
    base_clause += "\n<refine.ls_d_res_high.comparator>between" + \
      "</refine.ls_d_res_high.comparator>"
    if (d_min is not None):
      assert (d_min >= 0)
      print("  applying resolution cutoff d_min > %f" % d_min, file=log)
      base_clause += \
        "\n<refine.ls_d_res_high.min>%f</refine.ls_d_res_high.min>" % d_min
    if (d_max is not None):
      assert (d_max >= 0)
      print("  applying resolution cutoff d_min < %f" % d_max, file=log)
      base_clause += \
        "\n<refine.ls_d_res_high.max>%f</refine.ls_d_res_high.max>" % d_max
    queries.append(base_clause)
  if (protein_only):
    print("  excluding non-protein models", file=log)
    queries.append(
      "<queryType>org.pdb.query.simple.ChainTypeQuery</queryType>\n" +
      "<containsProtein>Y</containsProtein>")
  if (identity_cutoff is not None):
    print("  filtering based on sequence identity < %d %%" % \
      int(identity_cutoff), file=log)
    assert (identity_cutoff > 0) and (identity_cutoff <= 100)
    queries.append(
      "<queryType>org.pdb.query.simple.HomologueReductionQuery</queryType>\n"+
      "<identityCutoff>%d</identityCutoff>" % int(identity_cutoff))
  if (len(queries) == 0):
    raise RuntimeError("No queries specified!")
  queries_string = ""
  k = 0
  for clause in queries :
    queries_string += """
        <queryRefinement>
          <queryRefinementLevel>%d</queryRefinementLevel>
          <conjunctionType>and</conjunctionType>
          <orgPdbQuery>
            %s
          </orgPdbQuery>
        </queryRefinement>""" % (k, clause)
    k += 1
  query_str = """\
<orgPdbCompositeQuery version="1.0">
  %s
</orgPdbCompositeQuery>
""" % (queries_string)
  parsed = parseString(query_str)
  url_final = url_search
  if (sort_by_resolution):
    url_final += "/?sortfield=Resolution"
    print("  will sort by resolution", file=log)
  print("  executing HTTP request...", file=log)
  result = libtbx.utils.urlopen(url_final, query_str).read()
  return result.splitlines()

def sequence_search(sequence, **kwds):
  search_type = kwds.pop("search_type", "blast")
  expect = kwds.pop("expect", 0.01)
  min_identity = kwds.pop("min_identity", 0.0)
  if (min_identity <= 1):
    min_identity *= 100
  """
  Homology search for an amino acid sequence.  The advantage of using this
  service over the NCBI/EBI BLAST servers (in iotbx.pdb.fetch) is the ability
  to exclude non-Xray structures.
  """
  assert (search_type in ["blast", "fasta", "psiblast"])
  query_str = """\
    <queryType>org.pdb.query.simple.SequenceQuery</queryType>
      <description>Sequence Search (expect = %g, search = %s)</description>
      <sequence>%s</sequence>
      <eCutOff>%g</eCutOff>
      <sequenceIdentityCutoff>%g</sequenceIdentityCutoff>
    <searchTool>%s</searchTool>""" % (expect, search_type, sequence, expect,
      min_identity, search_type)
  return post_query(query_str, **kwds)

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
  >>> len(chemical_id_search("ZN", data_only=True, identity_cutoff=70))
  2874
  """
  assert (1 <= len(resname) <= 3)
  polymeric_type = kwds.pop("polymeric_type", "Any")
  assert (polymeric_type in ["Any", "Free", "Polymeric"])
  polymer_limit = "<polymericType>%s</polymericType>" % polymeric_type
  query_str = """\
    <queryType>org.pdb.query.simple.ChemCompIdQuery</queryType>
    <description>Chemical ID: %s</description>
    <chemCompId>%s</chemCompId>
    %s""" % (resname, resname, polymer_limit)
  return post_query(query_str, **kwds)

def get_custom_report_table(pdb_ids, columns, log=sys.stdout,
    prefix="dimStructure", check_for_missing=True):
  """Given a list of PDB IDs and a list of attribute identifiers, returns a
  Python list of lists for the IDs and attributes."""
  assert (len(columns) > 0)
  if (len(pdb_ids) == 0) : return []
  url_base = "http://www.rcsb.org/pdb/rest/customReport?"
  url = url_base + "pdbids=%s" % ",".join(pdb_ids)
  all_columns = ["structureId"] + columns
  url += "&customReportColumns=%s" % ",".join(all_columns)
  result = libtbx.utils.urlopen(url).read()
  # The RCSB's custom report follows this format (using the high-resolution
  # limit as an example):
  #
  # <?xml version='1.0' standalone='no' ?>
  # <dataset>
  #   <record>
  #     <dimStructure.structureId>1A0I</dimStructure.structureId>
  #     <dimStructure.highResolutionLimit>2.6</dimStructure.highResolutionLimit>
  #   </record>
  # </dataset>
  #
  # XXX note that 'dimStructure' is substituted for by 'dimEntity' in the
  # ligand report (and probably others), and that this table will not
  # necessarily have one row per PDB ID.
  xmlrec = parseString(result)
  table = []
  report_ids = set([])
  records = xmlrec.getElementsByTagName("record")
  for record in records :
    pdb_id_nodes = record.getElementsByTagName("%s.structureId" % prefix)
    assert (len(pdb_id_nodes) == 1), record.toxml()
    pdb_id = pdb_id_nodes[0].childNodes[0].data
    report_ids.add(pdb_id)
    row = [ pdb_id ] # pdb_id
    for col_name in columns :
      row_col = record.getElementsByTagName("%s.%s" % (prefix, col_name))
      assert (len(row_col) == 1)
      row.append(row_col[0].childNodes[0].data)
    table.append(row)
  if (check_for_missing):
    missing_ids = set(pdb_ids) - report_ids
    if (len(missing_ids) > 0):
      print("WARNING: missing report info for %d IDs:"%len(missing_ids), file=log)
      print("  %s" % " ".join(sorted(list(missing_ids))), file=log)
  return table

def get_high_resolution_for_structures(pdb_ids):
  return get_custom_report_table(pdb_ids, columns=["highResolutionLimit"])

def get_high_resolution_and_residue_count_for_structures(pdb_ids):
  return get_custom_report_table(pdb_ids, columns=["highResolutionLimit", 'residueCount'])

def get_ligand_info_for_structures(pdb_ids):
  """Return a list of ligands in the specified structure(s), including the
  SMILES strings."""
  return get_custom_report_table(pdb_ids,
    columns=["chainId","ligandId","ligandMolecularWeight","ligandFormula",
      "ligandName","ligandSmiles"],
    prefix="dimEntity",
    check_for_missing=False)
