
from libtbx.utils import Sorry
import xml.dom.minidom
import urllib2
import sys

def run (args) :
  try :
    pmid = int(args[0])
  except ValueError :
    raise Sorry("First argument must be a PubMed ID (i.e. integer).")
  url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
  query_args = "db=PubMed&retmode=xml&rettype=full&id=%d" % pmid
  data = urllib2.urlopen("%s?%s" % (url, query_args), timeout=2)
  #print data.read()
  xmlrec = xml.dom.minidom.parseString(data.read())
  article = xmlrec.getElementsByTagName("PubmedArticle")[0]
  authors = article.getElementsByTagName("Author")
  author_list = []
  for author in authors :
    initials = find_and_extract_string(author, "Initials")
    name = author.getElementsByTagName("LastName")[0].childNodes[0].data
    if (initials is not None) :
      name += " " + initials
    author_list.append(name)
  title = find_and_extract_string(article, "ArticleTitle")
  date = article.getElementsByTagName("PubDate")[0]
  year = find_and_extract_string(date, "Year")
  id_list = article.getElementsByTagName("ArticleIdList")
  doi = None
  id_list = article.getElementsByTagName("ArticleId")
  for id in id_list :
    if (id.getAttribute("IdType") == "doi") :
      doi = id.childNodes[0].data
  journal = find_and_extract_string(article, "MedlineTA")
  volume = find_and_extract_string(article, "Volume")
  pages = find_and_extract_string(article, "MedlinePgn")
  print """
# You should replace the article_id with a simple text string; the PMID will
# still be used to generate a link to PubMed.
citation {
  article_id = %d
  authors = %s
  title = %s
  year = %s
  journal = %s
  volume = %s
  pages = %s
  pmid = %d
  doi_id = %s
}
""" % (pmid, ", ".join(author_list), title, str(year), journal, volume, pages,
  pmid, doi)

def find_and_extract_string (node, name) :
  named_nodes = node.getElementsByTagName(name)
  if (len(named_nodes) > 0) :
    return named_nodes[0].childNodes[0].data
  return None

if __name__ == "__main__" :
  run(sys.argv[1:])
