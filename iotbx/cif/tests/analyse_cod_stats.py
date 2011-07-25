from libtbx import easy_pickle
from libtbx.utils import get_svn_revision, get_build_tag, plural_s, time_log
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
    self.structure_counts = {}
    self.skipped = set()
    self.timer = time_log("cif")
    self.max_delta = 0
    self.max_delta_id = ""
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
      self.timer.n += r.timer.n
      self.timer.accumulation += r.timer.accumulation
      self.max_delta = max(self.max_delta, r.max_delta)
      if self.max_delta == r.max_delta:
        self.max_delta_id = r.max_delta_id
      for n, count in r.structure_counts.items():
        if n in self.structure_counts:
          self.structure_counts[n] += count
        else:
          self.structure_counts[n] = count

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
    print >> out
    print >> out, self.timer.legend
    print >> out, self.timer.report()
    print >> out, "max delta: %.3g (%s)" %(self.max_delta, self.max_delta_id)

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
  a.show_summary()


if __name__ == '__main__':
  run()
