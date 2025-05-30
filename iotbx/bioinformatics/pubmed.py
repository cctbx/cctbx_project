"""Get pubmed reference"""
from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry
import libtbx.utils
from libtbx import slots_getstate_setstate
from xml.dom.minidom import parseString
from six.moves import urllib

def get_node_data(xml_node, node_name):
  child_nodes = xml_node.getElementsByTagName(node_name)
  return child_nodes[0].childNodes[0].data

def get_pubmed_xml(pmid):
  if isinstance(pmid, str):
    try :
      pmid = int(pmid)
    except ValueError :
      raise Sorry("PubMed IDs must be integers.")
  url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
  params = urllib.parse.urlencode({
    "id" : str(pmid),
    "rettype" : "full",
    "db" : "PubMed",
    "retmode" : "xml",
  })
  result = libtbx.utils.urlopen(url, params).read()
  xmlrec = parseString(result)
  articles = xmlrec.getElementsByTagName("PubmedArticle")
  if (len(articles) == 0):
    raise Sorry("No article with PMID %s found.")
  assert len(articles) == 1
  return article(articles[0])

class author(slots_getstate_setstate):
  __slots__ = ["last_name", "initials"]
  def __init__(self, xml_node):
    self.last_name = get_node_data(xml_node, "LastName").encode("utf-8")
    self.initials = get_node_data(xml_node, "Initials").encode("utf-8")

  def __str__(self):
    return ("%s %s" % (self.last_name, self.initials)).strip()

  def bibtex(self):
    initials = ""
    for i in list(self.initials):
      initials += "%s." % i
    return ("%s, %s" % (self.last_name, initials)).strip()

class article(slots_getstate_setstate):
  __slots__ = ["authors", "title", "year", "journal", "volume", "pages",
               "pmid", "doi"]
  def __init__(self, xmlrec):
    authors = xmlrec.getElementsByTagName("Author")
    self.authors = []
    for author_xml_node in authors :
      self.authors.append(author(author_xml_node))
    self.title = get_node_data(xmlrec, "ArticleTitle").encode("utf-8")
    self.year = get_node_data(xmlrec, "Year").encode("utf-8")
    self.journal = get_node_data(xmlrec, "MedlineTA").encode("utf-8")
    self.volume = get_node_data(xmlrec, "Volume").encode("utf-8")
    self.pages = get_node_data(xmlrec, "MedlinePgn").encode("utf-8")
    self.pmid = get_node_data(xmlrec, "PMID").encode("utf-8")
    self.doi = None
    article_ids = xmlrec.getElementsByTagName("ArticleId")
    for xml_node in article_ids :
      id_type = xml_node.getAttribute("IdType")
      if (id_type == "doi"):
        self.doi = xml_node.childNodes[0].data.encode("utf-8")

  def as_phenix_citation(self):
    doi = "None"
    if (self.doi is not None):
      doi = "\"%s\"" % self.doi
    return "\n".join([
      "citation {",
      "  article_id = %s" % self.pmid,
      "  authors = %s" % ", ".join([ str(a) for a in self.authors ]),
      "  title = %s" % self.title,
      "  journal = %s" % self.journal,
      "  volume = %s" % self.volume,
      "  pages = %s" % self.pages,
      "  year = %s" % self.year,
      "  doi_id = %s" % doi,
      "  pmid = %s" % self.pmid,
      "}"])

  def as_bibtex_citation(self):
    return "\n".join([
      "@Article{%s," % self.pmid,
      "  Author = {%s}," % " and ".join([ a.bibtex() for a in self.authors ]),
      "  Year = {%s}," % self.year,
      "  Title = {%s}," % self.title,
      "  Journal = {%s}," % self.journal,
      "  Volume = {%s}," % self.volume,
      "  Pages = {%s}," % self.pages,
      "}"])
