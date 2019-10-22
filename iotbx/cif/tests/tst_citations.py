from __future__ import absolute_import, division, print_function

from iotbx.cif.citations import citations_as_cif_block
from libtbx import citations

def test_format():
  article_ids = ['phenix2010', 'phenix.refine', 'phaser', 'molprobity']
  citation_list = [ citations.citations_db[article_id]
                    for article_id in article_ids ]
  cif_block = citations_as_cif_block(citation_list)
  assert list(cif_block['_citation.id']) == [
    'phenix2010', 'phenix.refine', 'phaser', 'molprobity']
  assert list(cif_block['_citation.journal_id_CSD']) == [
    '0766', '0766', '0228', '?']
  assert list(cif_block['_citation.journal_volume']) == ['66', '68', '40', '27']
  expected_keys = ('_citation.id', '_citation.title', '_citation.journal_abbrev',
                   '_citation.journal_volume', '_citation.page_first',
                   '_citation.page_last', '_citation.year',
                   '_citation.journal_id_ASTM', '_citation.journal_id_ISSN',
                   '_citation.journal_id_CSD')
  for key in expected_keys:
    assert key in cif_block

if (__name__ == '__main__'):
  test_format()
  print('OK')
