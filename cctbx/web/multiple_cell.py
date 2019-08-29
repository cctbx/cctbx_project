from __future__ import absolute_import, division, print_function
URL_cctbx_web = "cctbx_web.cgi"
URL_target_module = "explore_symmetry"

# this list was generated with generate_multiple_cell_stage2.py
multiple_cell_list = (
("C 4", (" C 4", " C 4")),
("C 41", (" C 4w", " C 4w")),
("C 42", (" C 4c", " C 4c")),
("C 43", (" C 4cw", " C 4cw")),
("F 4", (" F 4", " F 4")),
("F 41", (" F 4ad", " F 4d")),
("C -4", (" C -4", " C -4")),
("F -4", (" F -4", " F -4")),
("C 4/m", ("-C 4", "-C 4")),
("C 42/m", ("-C 4c", "-C 4c")),
("C 4/a :1", (" C 4a -1a", " C 4a -1a")),
("C 4/a :2", ("-C 4uv", "-C 4auv")),
("C 42/a :1", (" C 4ac -1ac", " C 4ac -1ac")),
("C 42/a :2", ("-C 4acuv", "-C 4wd")),
("F 4/m", ("-F 4", "-F 4")),
("F 41/d :1", (" F 4ad -1ad", " F 4d -1d")),
("F 41/d :2", ("-F 4auw", "-F 4vw")),
("C 4 2 2", (" C 4 2", " C 4 2")),
("C 4 2 21", (" C 4a 2", " C 4a 2")),
("C 41 2 2", (" C 4w 2w", " C 4w 2cw")),
("C 41 2 21", (" C 4aw 2c", " C 4aw 2")),
("C 42 2 2", (" C 4c 2c", " C 4c 2c")),
("C 42 2 21", (" C 4ac 2", " C 4ac 2")),
("C 43 2 2", (" C 4cw 2cw", " C 4cw 2w")),
("C 43 2 21", (" C 4acw 2c", " C 4acw 2")),
("F 4 2 2", (" F 4 2", " F 4 2")),
("F 41 2 2", (" F 4ad 2", " F 4d 2")),
("C 4 m m", (" C 4 -2", " C 4 -2")),
("C 4 m g1", (" C 4 -2a", " C 4 -2a")),
("C 42 m c", (" C 4c -2", " C 4c -2")),
("C 42 m g2", (" C 4ac -2", " C 4ac -2")),
("C 4 c c", (" C 4 -2c", " C 4 -2c")),
("C 4 c g2", (" C 4 -2ac", " C 4 -2ac")),
("C 42 c m", (" C 4c -2c", " C 4c -2c")),
("C 42 c g1", (" C 4c -2ac", " C 4c -2ac")),
("F 4 m m", (" F 4 -2", " F 4 -2")),
("F 4 m c", (" F 4 -2a", " F 4 -2a")),
("F 41 d m", (" F 4ad -2ad", " F 4d -2d")),
("F 41 d c", (" F 4ad -2d", " F 4d -2ad")),
("C -4 m 2", (" C -4 -2", " C -4 -2")),
("C -4 c 2", (" C -4 -2c", " C -4 -2c")),
("C -4 m 21", (" C -4 -2a", " C -4 -2a")),
("C -4 c 21", (" C -4 -2ac", " C -4 -2ac")),
("C -4 2 m", (" C -4 2", " C -4 2")),
("C -4 2 c", (" C -4 2c", " C -4 2c")),
("C -4 2 g1", (" C -4 2a", " C -4 2a")),
("C -4 2 g2", (" C -4 2ac", " C -4 2ac")),
("F -4 2 m", (" F -4 2", " F -4 2")),
("F -4 2 c", (" F -4 2a", " F -4 2a")),
("F -4 m 2", (" F -4 -2", " F -4 -2")),
("F -4 d 2", (" F -4 -2ad", " F -4 -2d")),
("C 4/m m m", ("-C 4 2", "-C 4 2")),
("C 4/m c c", ("-C 4 2c", "-C 4 2c")),
("C 4/a m g1 :1", (" C 4 2 -1a", " C 4 2 -1a")),
("C 4/a m g1 :2", ("-C 4uv 2a", "-C 4auv 2")),
("C 4/a c g2 :1", (" C 4 2 -1ac", " C 4 2 -1ac")),
("C 4/a c g2 :2", ("-C 4uv 2ac", "-C 4auv 2c")),
("C 4/m m g1", ("-C 4 2a", "-C 4 2a")),
("C 4/m c g2", ("-C 4 2ac", "-C 4 2ac")),
("C 4/a m m :1", (" C 4a 2 -1a", " C 4a 2 -1a")),
("C 4/a m m :2", ("-C 4uv 2", "-C 4auv 2a")),
("C 4/a c c :1", (" C 4a 2c -1a", " C 4a 2c -1a")),
("C 4/a c c :2", ("-C 4uv 2c", "-C 4auv 2ac")),
("C 42/m c m", ("-C 4c 2c", "-C 4c 2c")),
("C 42/m m c", ("-C 4c 2", "-C 4c 2")),
("C 42/a c g1 :1", (" C 4ac 2a -1ac", " C 4ac 2a -1ac")),
("C 42/a c g1 :2", ("-C 4wd 2ac", "-C 4acuv 2c")),
("C 42/a m g2 :1", (" C 4ac 2ac -1ac", " C 4ac 2ac -1ac")),
("C 42/a m g2 :2", ("-C 4wd 2a", "-C 4acuv 2")),
("C 42/m c g1", ("-C 4c 2ac", "-C 4c 2ac")),
("C 42/m m g2", ("-C 4ac 2", "-C 4ac 2")),
("C 42/a c m :1", (" C 4ac 2 -1ac", " C 4ac 2 -1ac")),
("C 42/a c m :2", ("-C 4wd 2c", "-C 4acuv 2ac")),
("C 42/a m c :1", (" C 4ac 2c -1ac", " C 4ac 2c -1ac")),
("C 42/a m c :2", ("-C 4wd 2", "-C 4acuv 2a")),
("F 4/m m m", ("-F 4 2", "-F 4 2")),
("F 4/m m c", ("-F 4 2a", "-F 4 2a")),
("F 41/d d m :1", (" F 4ad 2 -1ad", " F 4d 2 -1d")),
("F 41/d d m :2", ("-F 4uw 2ud", "-F 4ud 2ud")),
("F 41/d d c :1", (" F 4ad 2a -1ad", " F 4d 2a -1d")),
("F 41/d d c :2", ("-F 4uw 2vw", "-F 4ud 2vw")),
("H 3", (" H 3", " H 3", " H 3")),
("H 31", (" H 31", " H 31", " H 31")),
("H 32", (" H 32", " H 32", " H 32")),
("H -3", ("-H 3", "-H 3", "-H 3")),
("H 3 2 1", (" H 3 2\"", " H 3 2\"", " H 3 2\"")),
("H 3 1 2", (" H 3 2", " H 3 2", " H 3 2")),
("H 31 2 1", (" H 31 2\"", " H 31 2\" (0 0 2)", " H 31 2\" (0 0 4)")),
("H 31 1 2", (" H 31 2 (0 0 2)", " H 31 2 (0 0 4)", " H 31 2")),
("H 32 2 1", (" H 32 2\"", " H 32 2\" (0 0 4)", " H 32 2\" (0 0 2)")),
("H 32 1 2", (" H 32 2 (0 0 4)", " H 32 2 (0 0 2)", " H 32 2")),
("H 3 1 m", (" H 3 -2", " H 3 -2", " H 3 -2")),
("H 3 m 1", (" H 3 -2\"", " H 3 -2\"", " H 3 -2\"")),
("H 3 1 c", (" H 3 -2c", " H 3 -2c", " H 3 -2c")),
("H 3 c 1", (" H 3 -2\"c", " H 3 -2\"c", " H 3 -2\"c")),
("H -3 m 1", ("-H 3 2\"", "-H 3 2\"", "-H 3 2\"")),
("H -3 c 1", ("-H 3 2\"c", "-H 3 2\"c", "-H 3 2\"c")),
("H -3 1 m", ("-H 3 2", "-H 3 2", "-H 3 2")),
("H -3 1 c", ("-H 3 2c", "-H 3 2c", "-H 3 2c")),
("H 6", (" H 6", " H 6", " H 6")),
("H 61", (" H 61", " H 61", " H 61")),
("H 65", (" H 65", " H 65", " H 65")),
("H 62", (" H 62", " H 62", " H 62")),
("H 64", (" H 64", " H 64", " H 64")),
("H 63", (" H 6c", " H 6c", " H 6c")),
("H -6", (" H -6", " H -6", " H -6")),
("H 6/m", ("-H 6", "-H 6", "-H 6")),
("H 63/m", ("-H 6c", "-H 6c", "-H 6c")),
("H 6 2 2", (" H 6 2", " H 6 2", " H 6 2")),
("H 61 2 2", (" H 61 2 (0 0 4)", " H 61 2", " H 61 2 (0 0 2)")),
("H 65 2 2", (" H 65 2 (0 0 2)", " H 65 2", " H 65 2 (0 0 4)")),
("H 62 2 2", (" H 62 2 (0 0 2)", " H 62 2", " H 62 2 (0 0 4)")),
("H 64 2 2", (" H 64 2 (0 0 4)", " H 64 2", " H 64 2 (0 0 2)")),
("H 63 2 2", (" H 6c 2", " H 6c 2", " H 6c 2")),
("H 6 m m", (" H 6 -2", " H 6 -2", " H 6 -2")),
("H 6 c c", (" H 6 -2c", " H 6 -2c", " H 6 -2c")),
("H 63 m c", (" H 6c -2c", " H 6c -2c", " H 6c -2c")),
("H 63 c m", (" H 6c -2", " H 6c -2", " H 6c -2")),
("H -6 2 m", (" H -6 -2", " H -6 -2", " H -6 -2")),
("H -6 2 c", (" H -6c -2c", " H -6c -2c", " H -6c -2c")),
("H -6 m 2", (" H -6 2", " H -6 2", " H -6 2")),
("H -6 c 2", (" H -6c 2", " H -6c 2", " H -6c 2")),
("H 6/m m m", ("-H 6 2", "-H 6 2", "-H 6 2")),
("H 6/m c c", ("-H 6 2c", "-H 6 2c", "-H 6 2c")),
("H 63/m m c", ("-H 6c 2c", "-H 6c 2c", "-H 6c 2c")),
("H 63/m c m", ("-H 6c 2", "-H 6c 2", "-H 6c 2")),
)

import cgi
from six.moves import urllib

def escape(s):
  return cgi.escape(s, 1)

def LinkExploreSymmetry(Hall):
  return \
    "<a href=\"%s?target_module=%s&convention=Hall;sgsymbol=%s\">%s</a>" % (
     URL_cctbx_web,
     URL_target_module,
     urllib.parse.quote_plus(Hall),
     Hall)

def run():
  tetragonal_list = []
  trigonal_list = []
  hexagonal_list = []
  for HM, HallSymbols in multiple_cell_list:
    if (HM[0] != "H"):
      tetragonal_list.append((HM, HallSymbols))
    elif (not "6" in HM):
      trigonal_list.append((HM, HallSymbols))
    else:
      hexagonal_list.append((HM, HallSymbols))
  print(r"""<html>
<head>
<title>cctbx - Multiple cell C or F and triple cell H settings</title>
</head>
<body bgcolor=#ffffff>
<hr>
<h1>cctbx - Multiple cell C or F and triple cell H settings</h1>
<hr>
The Hermann-Mauguin symbols listed in the tables below are
not recognized by the class sgtbx::space_group_symbols,
but you may click on the Hall symbols in the right
columns to explore the symmetries.
<p>
The multiple cell and triple cell symbols were first
introduced in "Internationale Tabellen zur Bestimmung von
Kristallstrukturen" in 1935 (IT-1935). They were dropped in
the International Tables for Crystallography, Volume I from
1952 (IT-I1952) but kept in the comparative tables. In the
International Tables for Crystallography, Volume A from
1983 (IT-A1983), the multiple cell and triple cell symbols
reappear in the sub- and supergroup data, and are also
included in the <i>Synoptic Tables of Space-Group Symbols</i>
(Table 4.3.1).
<p>
IT-A1983 defines two transformation matrices for
transforming from a P- or I-representation in the
tetragonal system to C or F respectively. The matrices are
shown in the table headers below. For most settings the
choice of the transformation matrix does not effect the
result of the transformation. However, there are a number
of settings where the result of the transformation depends
on the choice of the matrix. For these settings IT-A1983
does not unambiguously define a specific space group
representation for the given Hermann-Mauguin symbol (see
e.g. F41 below). The situation is similar for the
transformation from a primitive trigonal or hexagonal
setting to a triple cell H setting. Unfortunately there is
no established notation for resolving these ambigous
definitions of IT-1983.
<p>
IT-1935 defines unambigous transformations for the multiple
cell and triple cell settings (p. 30):
<p>
<i>The equation for transforming from a P- or I-representation
(in the tetragonal system) to C or F respectively are:</i>
<ul>
<li><b>a</b>(C,F) = <b>a</b>(P,I) + <b>b</b>(P,I);
    <b>b</b>(C,F) = -<b>a</b>(P,I) + <b>b</b>(P,I)
</ul>
<p>
<i>For the change from C to H in the hexagonal system we have:</i>
<ul>
<li><b>a</b>(H) = <b>a</b>(C) + 2<b>b</b>(C);
    <b>b</b>(H) = -2<b>a</b>(C) - <b>b</b>(C)
</ul>
<p>
<b>c</b>(C,F) = <b>c</b>(P,I) and <b>c</b>(H) = <b>c</b>(C) are implied.<br>
<b>Note that IT-1935 uses the symbol C for the <i>primitive</i>
settings in the trigonal and hexagonal systems.</b>
<p>
Unfortunately, some primitive C settings in IT-1935 are
different from the primitive settings in IT-I1952 and
IT-A1983 (Schoenflies symbols D3^3 and D3^5). Therefore it
is not possible to derive all triple cell H settings of
IT-1935 by applying a consistent transformation to the
primitive settings in IT-A1983. It is currently not known
to me (rwgk) if a similar problem exists for the tetragonal
space groups.
<p>
It is beyond the scope of the sgtbx to provide a table of
historic space group symbols that are only referred to in
very old papers. A related problem is that the Protein Data
Bank (PDB) sometimes uses the symbol H to distinguish
between the hexagonal and the rhombohedral setting of
trigonal space groups (e.g. the symbol <i>H 3</i> in the
PDB corresponds to <i>R 3 with hexagonal axes</i> in
IT-A1983). It could therefore cause significant confusion
if the sgtbx interprets the PDB symbols according to the
IT-1935 conventions.
<p>
Considering the ambiguities and conflicting definitions, a
direct support for multiple cell and triple cell settings
is not included in the sgtbx. The tables below are
provided for the rare cases where information about these
symbols is required.
""")
  print("<hr>")
  print("<h2>Tetragonal multiple cell C or F</h2>")
  print("See section 4.3.4 in the International Tables for Crystallography,")
  print("Volume A, 1983 or later.")
  print("<p>")
  print("The Hall symbols in this table were generated with the program")
  print('<a href="http://xtal.crystal.uwa.edu.au/">Xtal 3.7.0</a>.')
  print("<p>")
  print("<table border=2 cellpadding=2>")
  print("<tr>")
  print("<th>Hermann-Mauguin<br>symbol")
  print("<th>Hall symbol<br>a-b,a+b,c")
  print("<th>Hall symbol<br>a+b,-a+b,c")
  for HM, HallSymbols in tetragonal_list:
    print("<tr>")
    print("<td>" + escape(HM))
    if (HallSymbols[0] == HallSymbols[1]):
      print("<td align=center colspan=2>" + LinkExploreSymmetry(HallSymbols[0]))
    else:
      for Hall in HallSymbols:
        print("<td align=center>" + LinkExploreSymmetry(Hall))
  print("</table>")
  for system, section, list in (
    ("Trigonal", "4.3.5", trigonal_list),
    ("Hexagonal", "4.3.5", hexagonal_list)):
    print("<hr>")
    print("<h2>%s triple cell H</h2>" % system)
    print("See section %s" % section)
    print("in the International Tables for Crystallography,")
    print("Volume A, 1983 or later.")
    print("<p>")
    print("The Hall symbols in this table were generated manually.")
    print("<p>")
    print("<table border=2 cellpadding=2>")
    print("<tr>")
    print("<th>Hermann-Mauguin<br>symbol")
    print("<th>Hall symbol<br>a-b,a+2b,c")
    print("<th>Hall symbol<br>2a+b,-a+b,c")
    print("<th>Hall symbol<br>a+2b,-2a-b,c")
    for HM, HallSymbols in list:
      print("<tr>")
      print("<td>" + escape(HM))
      if (    HallSymbols[0] == HallSymbols[1]
          and HallSymbols[0] == HallSymbols[2]):
        print("<td align=center colspan=3>" + LinkExploreSymmetry(
          HallSymbols[0]))
      else:
        for Hall in HallSymbols:
          print("<td align=center>" + LinkExploreSymmetry(Hall))
    print("</table>")
  print("""
<hr>
R.W. Grosse-Kunstleve, October 2001
<hr>
</body>
</html>""")

if (__name__ == "__main__"):
  run()
