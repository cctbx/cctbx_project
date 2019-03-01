'''
Functionality for handling citations
'''
from __future__ import absolute_import, division, print_function

from libtbx.citations import journals_db

# =============================================================================
def citations_as_cif_block(articles, cif_block=None):
  import iotbx.cif.model
  if cif_block is None:
    cif_block = iotbx.cif.model.block()
  def replace_none_with_question_mark(s):
    if s is None: return '?'
    return s
  citation_loop = iotbx.cif.model.loop(header=(
    '_citation.id', '_citation.title', '_citation.journal_abbrev',
    '_citation.journal_volume', '_citation.page_first', '_citation.page_last',
    '_citation.year', '_citation.journal_id_ASTM', '_citation.journal_id_ISSN',
    '_citation.journal_id_CSD'))
  for article in articles:
    if article.pages is None:
      first_page, last_page = "?", "?"
    else:
      pages = article.pages.split('-')
      first_page = pages[0]
      if len(pages) == 1:
        last_page = '?'
      else:
        assert len(pages) == 2
        last_page = pages[1]
    journal = journals_db.get(article.journal)
    assert journal is not None
    citation_loop.add_row(
      {'_citation.id': article.article_id,
       '_citation.title': article.title,
       '_citation.journal_abbrev': journal.abbrev_CAS,
       '_citation.journal_volume': article.volume,
       '_citation.page_first': first_page,
       '_citation.page_last': last_page,
       '_citation.year': article.year,
       '_citation.journal_id_ASTM':
         replace_none_with_question_mark(journal.id_ASTM),
       '_citation.journal_id_ISSN':
         replace_none_with_question_mark(journal.id_ISSN),
       '_citation.journal_id_CSD':
         replace_none_with_question_mark(journal.id_CSD),
       })
  cif_block.add_loop(citation_loop)
  return cif_block

# =============================================================================
# end
