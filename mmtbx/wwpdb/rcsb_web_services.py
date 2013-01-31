
"""
Module for querying the RCSB web server using the REST API, as described here:
http://www.rcsb.org/pdb/software/rest.do

There is some overlap with iotbx.pdb.fetch, which really should have gone here
instead, but this module is intended to be used in higher-level automation
pipelines.
"""

from __future__ import division
import urllib

url_base = "http://www.rcsb.org/pdb/rest"
url_search = url_base + "/search"

def post_query (query_xml, xray_only=True) :
  other_query = ""
  if (xray_only) :
    other_query = \
      "<mvStructure.expMethod.value>X-RAY</mvStructure.expMethod.value>"
  query_str = """\
<orgPdbQuery>
%s
%s
</orgPdbQuery>
""" % (query_xml, other_query)
  result = urllib.urlopen(url_search, query_str).read()
  return result.splitlines()

def sequence_search (sequence,
    file_name=None,
    search_type="blast",
    expect=0.01,
    xray_only=True) :
  """
  Homology search for an amino acid sequence.  The advantage of using this
  service over the NCBI/EBI BLAST servers (in iotbx.pdb.fetch) is the ability
  to exclude non-Xray structures.
  """
  assert (search_type in ["blast", "fasta", "psiblast"])
  query_str = """\
<queryType>org.pdb.query.simple.SequenceQuery</queryType>
<description>Sequence Search (Expectation Value = %g, Search Tool = BLAST)</description>
<sequence>%s</sequence>
<eCutOff>%g</eCutOff>
<searchTool>%s</searchTool>
""" % (expect, sequence, expect, search_type)
  return post_query(query_str, xray_only=xray_only)
