from __future__ import absolute_import, division, print_function

from io import StringIO

from libtbx import citations
from libtbx.test_utils import show_diff

def test_phil():
  citation = citations.citations_db['cctbx']
  assert(citation.authors ==
         'Grosse-Kunstleve RW, Sauter NK, Moriarty NW, Adams PD')

  journal = citations.journals_db['Nature']
  assert(journal.id_ASTM == 'NATUAS')

def test_sort():
  article_ids = ['phenix.polder', 'phenix2010', 'phenix.refine', 'phaser',
                 'molprobity']
  citation_list = [ citations.citations_db[article_id]
                    for article_id in article_ids ]

  # check unsorted first because sorting is in place
  with StringIO() as output:
    citations.show_citations(citation_list, out=output, sort_by_name=False)
    assert(output.getvalue().startswith('Liebschner'))

  with StringIO() as output:
    citations.show_citations(citation_list, out=output)
    assert(output.getvalue().startswith('Adams'))

def test_format():
  authors = citations.author_list_with_periods('Poon BK')
  assert(authors[0].count('B.K.') == 1)

  expected = """Adams PD, Afonine PV, Bunkoczi G, Chen VB, Davis IW, Echols N, Headd JJ, Hung LW, Kapral GJ, Grosse-Kunstleve RW, McCoy AJ, Moriarty NW, Oeffner R, Read RJ, Richardson DC, Richardson JS, Terwilliger TC, Zwart PH. (2010) PHENIX: a comprehensive Python-based system for macromolecular structure solution. Acta Cryst. D66:213-221."""
  citation = citations.citations_db['phenix2010']
  actual = citations.format_citation(citation)
  assert not show_diff(actual, expected)

  expected = """Grosse-Kunstleve R.W., Sauter N.K., Moriarty N.W., and Adams P.D. (2002). The Computational Crystallography Toolbox: crystallographic algorithms in a reusable software framework. J. Appl. Cryst. 35, 126-136."""
  citation = citations.citations_db['cctbx']
  actual = citations.format_citation_cell(citation)
  assert not show_diff(actual, expected)

  expected = """McCoy A.J., Grosse-Kunstleve R.W., Adams P.D., Winn M.D., Storoni L.C., & Read R.J. (2007). J Appl Crystallogr 40, 658-674."""
  citation = citations.citations_db['phaser']
  actual = citations.format_citation_iucr(citation)
  assert not show_diff(actual, expected)

  expected = """<b>Iterative-build OMIT maps: map improvement by iterative model building and refinement without model bias.</b> T.C. Terwilliger, R.W. Grosse-Kunstleve, P.V. Afonine, N.W. Moriarty, P.D. Adams, R.J. Read, P.H. Zwart, and L.-W. Hung. <a href="https://doi.org/doi:10.1107/S0907444908004319"><i>Acta Cryst.</i> D<b>64</b>, 515-524 (2008)</a>."""
  citation = citations.citations_db['autobuild_omit']
  actual = citations.format_citation_html(citation)
  assert not show_diff(actual, expected)

if (__name__ == '__main__'):
  test_phil()
  test_sort()
  test_format()
  print('OK')
