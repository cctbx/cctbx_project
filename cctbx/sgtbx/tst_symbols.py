# $Id$

from cctbx_boost.sgtbx import SpaceGroupSymbols

def show(sgsymbols):
  print sgsymbols.SgNumber()
  print sgsymbols.Schoenflies()
  print sgsymbols.Qualifier()
  print sgsymbols.Hermann_Mauguin()
  print sgsymbols.Extension()
  print sgsymbols.Hall()

import symbols

def loop_table(table, TableId, Label, Quiet = 0):
  for hm in table:
    if (hm in symbols.table_kriber_1_1_exclude): continue
    if (not Quiet): print ">>>", hm
    try:
      syms = SpaceGroupSymbols(hm, TableId)
      if (not Quiet): show(syms)
    except RuntimeError, e:
      print hm, Label, e.args

def loop_tables(Quiet = 0):
  for table in symbols.tables.keys():
    print 'Next Table:', symbols.tables[table][0]
    loop_table(table,
               symbols.tables[table][1],
               symbols.tables[table][0],
               Quiet)

def check_icsd_names():
  for icsd_pair in symbols.table_icsd_names:
    print icsd_pair
    try:
      print ' ', icsd_pair[0]
      syms1 = SpaceGroupSymbols(icsd_pair[0])
      print ' ', icsd_pair[1]
      syms2 = SpaceGroupSymbols(icsd_pair[1])
      if (syms1.SgNumber() != syms2.SgNumber()):
        print 'Error:', syms1.SgNumber(), syms2.SgNumber()
    except RuntimeError, e:
      print icsd_pair, e.args

def weird_symbols():
  symbols = (
    "  hall p 3  ",
    " hall:   p 4",
    " hall :   p 6",
    " p    1 21/n     1",
    "C  6  v  #  1",
    "1  2",
  )
  loop_table(symbols, "", "weird")

def compare_syms(syms1, syms2):
  return (    syms1.SgNumber() == syms2.SgNumber()
          and syms1.Schoenflies() == syms2.Schoenflies()
          and syms1.Qualifier() == syms2.Qualifier()
          and syms1.Hermann_Mauguin() == syms2.Hermann_Mauguin()
          and syms1.Extension() == syms2.Extension()
          and syms1.Hall() == syms2.Hall())

def recycle(TableId = ""):
  for SgNumber in xrange(1, 231):
    syms1 = SpaceGroupSymbols(SgNumber, "", TableId)
    syms2 = SpaceGroupSymbols(str(SgNumber), TableId)
    if (not compare_syms(syms1, syms2)):
      print 'Error:', syms1.SgNumber(), syms2.SgNumber()
    syms2 = SpaceGroupSymbols(syms1.Schoenflies(), TableId)
    if (not compare_syms(syms1, syms2)):
      print 'Error:', syms1.SgNumber(), syms2.SgNumber()
    syms2 = SpaceGroupSymbols(syms1.Hermann_Mauguin(), TableId)
    if (not compare_syms(syms1, syms2)):
      print 'Error:', syms1.SgNumber(), syms2.SgNumber()
    if (syms1.Extension() != ""):
      syms2 = SpaceGroupSymbols(SgNumber, syms1.Extension(), TableId)
      if (not compare_syms(syms1, syms2)):
        print 'Error:', syms1.SgNumber(), syms2.SgNumber()
      syms2 = SpaceGroupSymbols(str(SgNumber) + ":" + syms1.Extension(), TableId)
      if (not compare_syms(syms1, syms2)):
        print 'Error:', syms1.SgNumber(), syms2.SgNumber()

#for i in xrange(10000):
#  loop_tables(1)

#check_icsd_names()
#weird_symbols()
recycle()
recycle("A")
recycle("I")
