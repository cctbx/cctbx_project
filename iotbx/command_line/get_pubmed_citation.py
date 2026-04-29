"""Fetch a citation for a PubMed article, output in either the format used
internally in Phenix, or BibText"""

from __future__ import absolute_import, division, print_function
from iotbx.bioinformatics import pubmed
import iotbx.phil
import sys

master_phil = """
pmid = None
  .type = int
  .multiple = True
bibtex = False
  .type = bool
"""

def run(args, out=sys.stdout):
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_phil,
    integer_def="pmid",
    usage_string="""\
iotbx.get_pubmed_citation PMID

Fetch a citation for a PubMed article, output in either the format used
internally in Phenix, or BibText (if bibtex=True or --bibtext specified).
""")
  params = cmdline.work.extract()
  assert (params.pmid is not None) and (len(params.pmid) > 0)
  for pmid in params.pmid :
    article = pubmed.get_pubmed_xml(pmid)
    if (params.bibtex):
      print(article.as_bibtex_citation(), file=out)
    else :
      print(article.as_phenix_citation(), file=out)

if (__name__ == "__main__"):
  run(sys.argv[1:])
