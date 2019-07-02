from __future__ import absolute_import, division, print_function
from libtbx import easy_pickle
from libtbx.utils import get_svn_revision, get_build_tag, plural_s
from operator import itemgetter
import glob, os, sys
from six.moves import cStringIO as StringIO

class analyse(object):

  def __init__(self):
    self.n_hkl = 0
    self.n_cif = 0
    self.n_hkl_cif_pairs = 0
    self.parsing_errors = {}
    self.build_errors = {}
    self.ignored_errors = {}
    self.skipped = set()
    g = glob.glob("result_*.pickle")
    for f in g:
      r = easy_pickle.load(f)
      self.n_hkl = r.n_hkl
      self.n_cif = r.n_cif
      self.n_hkl_cif_pairs = r.n_hkl_cif_pairs
      self.parsing_errors.update(r.parsing_errors)
      self.build_errors.update(r.build_errors)
      self.ignored_errors.update(r.ignored_errors)
      self.skipped.update(r.skipped)

  def show_summary(self, out=None):
    if out is None: out = sys.stdout
    n_total = self.n_cif + self.n_hkl
    print("Number of cif: %7i" % self.n_cif, file=out)
    print("Number of hkl: %7i" % self.n_hkl, file=out)
    print("Number of hkl+cif pairs: ", self.n_hkl_cif_pairs, file=out)
    print(file=out)

    print("%i parsing error%s" % plural_s(len(self.parsing_errors)), file=out)
    print("%i exception%s" % plural_s(len(self.build_errors)), file=out)
    print("%i ignored exception%s" % plural_s(len(self.ignored_errors)), file=out)
    print("%i skipping" % len(self.skipped), file=out)
    print(file=out)
    rev = get_svn_revision()
    cod_path = os.environ.get("COD_SVN_WORKING_COPY")
    cod_rev = None
    if cod_path is not None:
      cod_rev = get_svn_revision(path=cod_path)
    build_tag = get_build_tag()
    if cod_rev is not None:
      print("COD svn revision: %i" %cod_rev, file=out)
    if rev is not None:
      print("cctbx svn revision: %i" %rev, file=out)
    if build_tag is not None:
      print("cctbx build tag: ", build_tag, file=out)

  def show_all(self, out=None):
    if out is None: out = sys.stdout
    self.show_summary(out=out)
    print(file=out)
    self.show_skipping(out=out)
    print(file=out)
    self.show_exceptions(out=out)

  def show_parsing_errors(self, out=None):
    if out is None: out = sys.stdout
    for cod_id, error in sorted(
      sorted(self.parsing_errors.items()), key=itemgetter(1)):
      print("%s: %s" % (cod_id, error), file=out)

  def show_exceptions(self, out=None):
    if out is None: out = sys.stdout
    for cod_id, error in sorted(
      sorted(self.build_errors.items()), key=itemgetter(1)):
      print("%s: %s" % (cod_id, error), file=out)

  def show_skipping(self, out=None):
    if out is None: out = sys.stdout
    for s in sorted(self.skipped):
      print("SKIPPED: ", s, file=out)

  def as_html(self):
    with open("parsing_errors", "w") as out:
      self.show_parsing_errors(out=out)
    with open("exceptions", "w") as out:
      self.show_exceptions(out=out)
    with open("skipping", "w") as out:
      self.show_skipping(out=out)
    s = StringIO()
    self.show_summary(out=s)
    with open("summary.html", "w") as out:
      print(html % s.getvalue(), file=out)

html = """\
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">

<head>
<META http-equiv=Content-Type content="text/html; charset=utf-8">
<title>COD stats</title>
</head>

<body>

<hr>

<pre>
%s
</pre>

<hr>

<a href="parsing_errors">[parsing errors]</a>
<a href="exceptions">[exceptions]</a>
<a href="skipping">[skipping]</a>

<hr>

<a href="/cod_stats/current">[directory listing]</a>

<hr>

</body>
"""

def run():
  a = analyse()
  a.as_html()


if __name__ == '__main__':
  run()
