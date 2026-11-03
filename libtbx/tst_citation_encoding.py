from __future__ import absolute_import, division, print_function

"""
Regression test for citation output on narrow locale codecs (Windows cp874/cp1252).

Citation author names contain Latin diacritics (Bunkóczi, Brünger, Söding, Millán,
Pražnikar, Žídek, Iwaï, ...) and non-ASCII punctuation (en-dash U+2013, Unicode hyphen
U+2010, ellipsis U+2026). Writing these to a stream whose codec cannot encode them raised
UnicodeEncodeError and killed the job (observed in AutoSol_display_summary on a Thai/cp874
Windows machine: char '\\xf3' = 'ó' in "Bunkóczi").

citation_text_for_stream(text, out) contract (citation display only; lossy):
  * out is None -> inspect sys.stdout (where print(file=None) goes);
  * out has no .encoding -> unconstrained Unicode sink (e.g. StringIO): return unchanged;
  * out has a codec -> transliterate ONLY characters that codec cannot strictly encode;
  * result is always strictly encodable in the target codec;
  * raises no encoding-conversion exception for valid str input + codec metadata
    (it does NOT suppress errors from the stream's own write()).

All assertions pin the codec explicitly, so this runs on Mac/Linux CI without a Thai box.
"""

import io
import sys
from libtbx.citations import citation_text_for_stream as cts
from libtbx.citations import _downgrade_char

SAMPLES = [
  u"Bunk\u00f3czi G, Read RJ",   # ó  U+00F3  (the reported crash character)
  u"Br\u00fcnger AT",            # ü  U+00FC
  u"S\u00f6ding J",              # ö  U+00F6
  u"Mill\u00e1n C",              # á  U+00E1
  u"Pra\u017enikar J",           # ž  U+017E
  u"\u017d\u00eddek A",          # Ž U+017D + í U+00ED
  u"Iwa\u00ef PR",               # ï  U+00EF
  u"pages 679\u2013682",         # en-dash U+2013
  u"low\u2010resolution",        # Unicode hyphen U+2010
  u"see ref\u2026",              # ellipsis U+2026
]

class _Stream(io.TextIOBase):
  """Stream with a chosen .encoding whose write() encodes with it, like a real file.
  Only ever constructed with a real codec name (never None)."""
  def __init__(self, enc):
    assert enc is not None
    self._b = io.BytesIO(); self._e = enc
  @property
  def encoding(self): return self._e
  def write(self, s):
    self._b.write(s.encode(self._e)); return len(s)
  def getvalue(self): return self._b.getvalue()


# --------------------------------------------------------------------------- #
# EFFICACY: prove it fixes the reported problem                               #
# --------------------------------------------------------------------------- #

def test_efficacy():
  bunkoczi = u"Bunk\u00f3czi G, Read RJ"

  # F1 (control): the raw citation text is unencodable on cp874. PASSES before
  # AND after the patch -- it merely confirms the test models the real crash.
  try:
    bunkoczi.encode("cp874"); raise AssertionError("expected cp874 UnicodeEncodeError")
  except UnicodeEncodeError:
    pass

  # F2: through the helper, the same content becomes strictly encodable on the
  # narrow codec. FAILS before the patch (helper absent), passes after.
  for enc in ("cp874", "cp1252"):
    cts(bunkoczi, _Stream(enc)).encode(enc)   # must not raise

  # F3: specific downgrades on cp874 are correct.
  st = _Stream("cp874")
  assert cts(u"Bunk\u00f3czi", st)  == "Bunkoczi"
  assert cts(u"Br\u00fcnger",  st)  == "Brunger"
  assert cts(u"S\u00f6ding",   st)  == "Soding"
  assert cts(u"Mill\u00e1n",   st)  == "Millan"
  assert cts(u"Pra\u017enikar",st)  == "Praznikar"
  assert cts(u"\u017d\u00eddek",st) == "Zidek"
  assert cts(u"Iwa\u00ef",     st)  == "Iwai"

  # Punctuation must be REPLACED, never DELETED (anti-regression vs ascii/ignore):
  # a deleted dash would fuse digits ("679-682" -> "679682") and change meaning.
  for enc in ("ascii", "cp1252", "cp874"):
    st2 = _Stream(enc)
    for dash in (u"\u2010", u"\u2011", u"\u2012", u"\u2013", u"\u2014"):
      r = cts(u"679" + dash + u"682", st2)
      assert "679682" not in r, "dash deleted -> %r (enc=%s)" % (r, enc)
      r.encode(enc)
  # ellipsis maps to '...' on a codec that lacks it
  assert cts(u"ref\u2026", _Stream("ascii")) == "ref..."


# --------------------------------------------------------------------------- #
# SAFETY                                                                      #
# --------------------------------------------------------------------------- #

def test_safety_utf8_is_noop():
  # On UTF-8, byte-identical (diacritics preserved).
  st = _Stream("utf-8")
  for txt in SAMPLES:
    assert cts(txt, st) == txt, txt

def test_safety_stringio_is_unicode_sink():
  # StringIO has no .encoding -> unconstrained Unicode sink -> unchanged (Policy A),
  # independent of the host locale. This is the key correctness fix from review.
  buf = io.StringIO()
  text = u"Bunk\u00f3czi \u2013 \u017d\u00eddek"
  assert cts(text, buf) == text

def test_safety_none_means_stdout():
  # out=None must inspect sys.stdout, not the process locale.
  saved = sys.stdout
  try:
    sys.stdout = _Stream("cp874")
    assert cts(u"Bunk\u00f3czi", None) == "Bunkoczi"
    sys.stdout = _Stream("utf-8")
    assert cts(u"Bunk\u00f3czi", None) == u"Bunk\u00f3czi"
  finally:
    sys.stdout = saved

def test_safety_minimal_downgrade_cp1252():
  # cp1252 contains Ž/ž/í and the en-dash; only U+2010 is absent.
  st = _Stream("cp1252")
  assert cts(u"\u017d\u00eddek", st) == u"\u017d\u00eddek"   # kept
  assert cts(u"Bunk\u00f3czi", st)   == u"Bunk\u00f3czi"     # kept
  assert cts(u"a\u2010b", st)        == "a-b"                # U+2010 not in cp1252
  for txt in SAMPLES:
    cts(txt, st).encode("cp1252")

def test_safety_type_and_result_invariants():
  # Non-str raises TypeError (declared contract); None passes through.
  assert cts(None, _Stream("cp874")) is None
  raised = False
  try:
    cts(123, _Stream("cp874"))
  except TypeError:
    raised = True
  assert raised, "expected TypeError on non-str input"
  # Result always encodes in the target codec, for every codec/sample.
  for enc in ("utf-8", "cp874", "cp1252", "ascii"):
    st = _Stream(enc)
    for txt in SAMPLES:
      cts(txt, st).encode(enc)
  # Unmappable script (no ASCII base) -> '?', still encodable, no crash.
  assert cts(u"\u0410\u0411", _Stream("cp874")) == "??"   # Cyrillic AB


def test_downgrade_char_is_ascii():
  # Every non-ASCII character in the citation inventory must downgrade to ASCII
  # (so the transliterated result is always encodable in any ASCII-superset codec).
  inventory = [
    u"\u00e1", u"\u00e9", u"\u00ed", u"\u00ef", u"\u00f3", u"\u00f6",
    u"\u00fc", u"\u017e", u"\u017d",                         # letters
    u"\u2010", u"\u2011", u"\u2012", u"\u2013", u"\u2014",  # dashes
    u"\u2018", u"\u2019", u"\u201c", u"\u201d", u"\u2026",  # quotes/ellipsis
  ]
  for ch in inventory:
    r = _downgrade_char(ch)
    assert r.isascii() and r != "", "non-ASCII or empty downgrade for %r -> %r" % (ch, r)

def test_invalid_codec_name_surfaces():
  # An invalid out.encoding is a defective stream, not a locale condition: it
  # must NOT be silently returned through; LookupError should surface.
  class _Bad:
    encoding = "not-a-real-codec"
  raised = False
  try:
    cts(u"Bunk\u00f3czi", _Bad())
  except LookupError:
    raised = True
  assert raised, "invalid codec name should raise LookupError, not pass through"

def run():
  test_efficacy()
  test_safety_utf8_is_noop()
  test_safety_stringio_is_unicode_sink()
  test_safety_none_means_stdout()
  test_safety_minimal_downgrade_cp1252()
  test_safety_type_and_result_invariants()
  test_downgrade_char_is_ascii()
  test_invalid_codec_name_surfaces()
  print("OK")

if __name__ == "__main__":
  run()
