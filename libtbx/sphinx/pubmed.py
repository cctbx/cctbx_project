from __future__ import division

from docutils.parsers.rst import Directive
from docutils.parsers.rst import directives
from docutils import statemachine


def setup(app):
  app.add_directive('pubmed', PubMedDirective)


class PubMedDirective(Directive):

  # this disables content in the directive
  has_content = False
  required_arguments = 1
  optional_arguments = 1
  final_argument_whitespace = True
  option_spec = {'reprint-url': directives.unicode_code}

  def run(self):

    PMID = self.arguments[0]
    description = self.arguments[1] # this isn't used at all
    reprint_url = self.options.get('reprint-url', None)

    from Bio import Entrez
    Entrez.email="NKSauter@lbl.gov"

    handle = Entrez.efetch(
      db="pubmed",id=PMID, # must be a quoted string of comma-separated PMIDs
      retmode="xml")
    XML = Entrez.read(handle)

    text = []

    for i in xrange(len(XML)):
      # Title/doi link:
      possible_doi = [ idx for idx in XML[i]["PubmedData"]["ArticleIdList"]
                       if idx.attributes["IdType"]=="doi" ]
      article = XML[i]["MedlineCitation"]["Article"]
      get_title = article["ArticleTitle"]
      get_title = get_title.strip(".")
      if len(possible_doi) > 0:
        text.append("| `%s <http://dx.doi.org/%s>`_" %(
          get_title.encode("ascii","xmlcharrefreplace"), possible_doi[0]))

      # Author list
      authors = [ " ".join([elem["LastName"],elem["Initials"]])
                  for elem in article["AuthorList"] ]
      text.append(
        "| %s." %(", ".join(authors).encode("ascii","xmlcharrefreplace")))

      # Journal reference
      journal = article["Journal"]
      journal_text = "| *%s*" %(journal["ISOAbbreviation"])
      issue = journal["JournalIssue"]
      if issue.has_key("Volume"):
        journal_text += " **%s**, %s" %(
          issue["Volume"], article["Pagination"]["MedlinePgn"] )
      date = issue["PubDate"]
      if date.has_key("Month"):
        journal_text += " (%s %s %s)."%(date.get("Day",1), date["Month"], date["Year"])
      else:
        journal_text += " (%s)"%(date["Year"])
      journal_text += " [PMID:%s]"%PMID
      possible_pmc = [ idx for idx in XML[i]["PubmedData"]["ArticleIdList"]
                       if idx.attributes["IdType"]=="pmc" ]
      if len(possible_pmc) > 0:
        journal_text += " [PMC reprint: `%s <http://ncbi.nlm.nih.gov/pmc/articles/%s/>`_]"%(
          str(possible_pmc[0]), str(possible_pmc[0]) )
      if reprint_url is not None:
        journal_text += " `(Reprint) <%s>`_" %(reprint_url)
      text.append(journal_text)

    print "\n".join(text)

    # insert rst
    source = self.state_machine.input_lines.source(
      self.lineno - self.state_machine.input_offset - 1)
    lines = statemachine.string2lines(
      "\n".join(text), self.state.document.settings.tab_width, convert_whitespace=True)
    self.state_machine.insert_input(lines, source)
    return []
