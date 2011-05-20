from libtbx import easy_pickle
from libtbx.utils import get_svn_revision, get_build_tag, plural_s
from operator import itemgetter
import glob, os, sys
from cStringIO import StringIO

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
    print >> out, "Number of cif: %7i" % self.n_cif
    print >> out, "Number of hkl: %7i" % self.n_hkl
    print >> out, "Number of hkl+cif pairs: ", self.n_hkl_cif_pairs
    print >> out

    print >> out, "%i parsing error%s" % plural_s(len(self.parsing_errors))
    print >> out, "%i exception%s" % plural_s(len(self.build_errors))
    print >> out, "%i ignored exception%s" % plural_s(len(self.ignored_errors))
    print >> out, "%i skipping" % len(self.skipped)
    print >> out
    rev = get_svn_revision()
    cod_path = os.environ.get("COD_SVN_WORKING_COPY")
    cod_rev = None
    if cod_path is not None:
      cod_rev = get_svn_revision(path=cod_path)
    build_tag = get_build_tag()
    if cod_rev is not None:
      print >> out, "COD svn revision: %i" %cod_rev
    if rev is not None:
      print >> out, "cctbx svn revision: %i" %rev
    if build_tag is not None:
      print >> out, "cctbx build tag: ", build_tag

  def show_all(self, out=None):
    if out is None: out = sys.stdout
    self.show_summary(out=out)
    print >> out
    self.show_skipping(out=out)
    print >> out
    self.show_exceptions(out=out)

  def show_parsing_errors(self, out=None):
    if out is None: out = sys.stdout
    for cod_id, error in sorted(
      sorted(self.parsing_errors.items()), key=itemgetter(1)):
      print >> out, "%s: %s" % (cod_id, error)

  def show_exceptions(self, out=None):
    if out is None: out = sys.stdout
    for cod_id, error in sorted(
      sorted(self.build_errors.items()), key=itemgetter(1)):
      print >> out, "%s: %s" % (cod_id, error)

  def show_skipping(self, out=None):
    if out is None: out = sys.stdout
    for s in sorted(self.skipped):
      print >> out, "SKIPPED: ", s

  def as_html(self):
    with open("parsing_errors", "wb") as out:
      self.show_parsing_errors(out=out)
    with open("exceptions", "wb") as out:
      self.show_exceptions(out=out)
    with open("skipping", "wb") as out:
      self.show_skipping(out=out)
    s = StringIO()
    self.show_summary(out=s)
    with open("summary.html", "wb") as out:
      print >> out, html % s.getvalue()

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
