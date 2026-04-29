"""Given a PDB ID (and an email address, as required by the NCBI), will attempt
to identify the PubMed ID of the primary citation for the structure, download
the article XML from NCBI, and print a bibliography entry.
"""

# LIBTBX_SET_DISPATCHER_NAME iotbx.pdb.get_citation

from __future__ import absolute_import, division, print_function
from libtbx.utils import to_unicode, Sorry, Usage
from xml.dom.minidom import parseString
import sys

master_phil_string = """
pdb_id = None
  .type = str
email = None
  .type = str
"""

def run(args, out=sys.stdout):
  if (len(args) == 0) or ("--help" in args):
    raise Usage("""\
iotbx.pdb.get_citation ID EMAIL

Given a PDB ID (and an email address, as required by the NCBI), will attempt
to identify the PubMed ID of the primary citation for the structure, download
the article XML from NCBI, and print a bibliography entry.
""")
  import iotbx.pdb.fetch
  import iotbx.phil
  try :
    from Bio import Entrez
  except ImportError :
    raise Sorry("BioPython not installed.")
  class _cmdline(iotbx.phil.process_command_line_with_files):
    def process_other(self, arg):
      if (iotbx.pdb.fetch.valid_pdb_id(arg)):
        return iotbx.phil.parse("""pdb_id=%s""" % arg)
      elif (arg.count("@") == 1):
        return iotbx.phil.parse("""email=\"%s\"""" % arg)
      return False
  cmdline = _cmdline(
    args=args,
    master_phil_string=master_phil_string,
    integer_def="pmid")
  params = cmdline.work.extract()
  validate_params(params)
  # XXX should probably use mmCIF here
  pdb_data = iotbx.pdb.fetch.fetch(id=params.pdb_id)
  pmid = None
  for line in pdb_data.readlines():
    if (line.startswith("JRNL")):
      fields = line.split()
      if (fields[1] == "PMID"):
        pmid = int(fields[2].strip())
        break
  if (pmid is None):
    raise Sorry("Couldn't extract PMID from PDB header.")
  # This looks gross, but it's much safer than trying to do any more parsing
  # of unstructured text (even when using mmCIF)
  print("PMID = %s" % pmid, file=out)
  Entrez.email = params.email
  data = Entrez.efetch(db="pubmed", id=str(pmid), rettype="xml",
    retmode="text")
  def get_node_data(xml_node, node_name):
    child_nodes = xml_node.getElementsByTagName(node_name)
    return to_unicode(child_nodes[0].childNodes[0].data)
  xmlrec = parseString(data.read())
  articles = xmlrec.getElementsByTagName("PubmedArticle")
  assert (len(articles) == 1)
  article = articles[0]
  authors = article.getElementsByTagName("Author")
  author_list = []
  for author in authors :
    name = get_node_data(author, "LastName") + ", " + \
      ".".join(get_node_data(author, "Initials")) + "."
    author_list.append(name)
  people = ", ".join(author_list[:-1]) + " & " + author_list[-1]

  title = get_node_data(article, "ArticleTitle")
  if (not title.endswith(".")):
    title += "."
  date = article.getElementsByTagName("PubDate")[0]
  year = get_node_data(article, "Year")
  journal = get_node_data(article, "MedlineTA")
  volume = get_node_data(article, "Volume")
  pages = get_node_data(article, "MedlinePgn")
  # TODO other formats
  print("%s (%s).  %s %s %s, %s." % \
    (people.encode('ascii', 'ignore'), year, title.encode('ascii', 'ignore'),
     journal, volume, pages), file=out)

def validate_params(params):
  if (params.pdb_id is None):
    raise Sorry("PDB ID not specified!")
  if (params.email is None):
    raise Sorry("Please give a valid email address.")
  return True

if (__name__ == "__main__"):
  run(sys.argv[1:])
