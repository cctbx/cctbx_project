from cctbx_boost import uctbx
from cctbx_boost import sgtbx
from cctbx.development import debug_utils

import symbols

def loop_table(table, TableId, verbose = 0):
  for hm in table:
    if (hm in symbols.table_kriber_1_1_exclude): continue
    if (verbose): print ">>>", hm
    syms = sgtbx.SpaceGroupSymbols(hm, TableId)

def test_symbols_from_tables(verbose = 0):
  for table in symbols.tables.keys():
    if (verbose): print 'Next Table:', symbols.tables[table][0]
    loop_table(table,
               symbols.tables[table][1],
               verbose)

def test_weird_symbols(verbose = 0):
  symbols = (
    "  hall p 3  ",
    " hall:   p 4",
    " hall :   p 6",
    " p    1 21/n     1",
    "C  6  v  #  1",
    "1  2",
  )
  loop_table(symbols, "", verbose)

def test_icsd_symbols(verbose = 0):
  for icsd_pair in symbols.table_icsd_names:
    if (verbose): print icsd_pair
    syms = sgtbx.SpaceGroupSymbols(icsd_pair[0])
    other = str(syms.SgNumber())
    if (syms.Qualifier()):
      other += ":" + syms.Qualifier()
    elif (syms.Extension()):
      other += ":" + syms.Extension()
    assert other == icsd_pair[1], (icsd_pair, other)

def compare_syms(syms1, syms2):
  return (    syms1.SgNumber() == syms2.SgNumber()
          and syms1.Schoenflies() == syms2.Schoenflies()
          and syms1.Qualifier() == syms2.Qualifier()
          and syms1.Hermann_Mauguin() == syms2.Hermann_Mauguin()
          and syms1.Extension() == syms2.Extension()
          and syms1.Hall() == syms2.Hall())

def assert_syms(syms1, syms2):
  assert compare_syms(syms1, syms2), (
    syms1.ExtendedHermann_Mauguin(), syms2.ExtendedHermann_Mauguin())

def test_recycle_space_group_symbols(TableId = ""):
  for SgNumber in xrange(1, 231):
    syms1 = sgtbx.SpaceGroupSymbols(SgNumber, "", TableId)
    syms2 = sgtbx.SpaceGroupSymbols(str(SgNumber), TableId)
    assert_syms(syms1, syms2)
    syms2 = sgtbx.SpaceGroupSymbols(syms1.Schoenflies(), TableId)
    assert_syms(syms1, syms2)
    syms2 = sgtbx.SpaceGroupSymbols(syms1.Hermann_Mauguin(), TableId)
    assert_syms(syms1, syms2)
    if (syms1.Extension() != ""):
      syms2 = sgtbx.SpaceGroupSymbols(SgNumber, syms1.Extension(), TableId)
      assert_syms(syms1, syms2)
      syms2 = sgtbx.SpaceGroupSymbols(
        str(SgNumber) + ":" + syms1.Extension(), TableId)
      assert_syms(syms1, syms2)

def test_rtmx():
  print sgtbx.STBF
  print sgtbx.CRBF
  print sgtbx.CTBF
  s = sgtbx.parse_string("x+1/2,y,z")
  print s.string()
  print s.where()
  m = sgtbx.RTMx(s)
  print m.as_xyz(0, 0, "xyz", ",")
  print m.as_xyz(0, 1, "xyz", ",")
  print m.as_tuple()
  print m.as_tuple(12, 72)
  print m.modPositive()
  mm = m * m
  print mm.as_xyz()
  mi = m.inverse()
  print mi.as_xyz()
  print mi.as_xyz(0, 1, "xyz", ",")
  print (m * mi).as_xyz()
  m = sgtbx.RTMx("2*x,2*y,2*z", "", 72, 12)
  mi = m.inverse()
  print mi.isValid()
  print mi.as_xyz()
  m = sgtbx.RTMx("-x,-y,z")
  print (m * m).as_xyz()
  m = sgtbx.RTMx("y+z-1/3,x+z-1/3,x+y", "", 2, 3)
  print m
  print m.inverse()
  assert (m * m.inverse()).isUnit()
  assert (m.multiply(m.inverse())).isUnit()
  mi = m.inverse_with_cancel()
  print mi
  mmi = m.multiply(mi)
  print mmi
  assert mmi.isUnit()

def test_cancel_and_tidy():
  SgOps = sgtbx.SpaceGroup("P 61")
  for m in SgOps:
    print m.as_xyz()
    print m.cancel().TBF()
  SgOps.makeTidy()
  for m in SgOps:
    print m.as_xyz()

def test_space_group_comparison():
  SgOps1 = sgtbx.SpaceGroup("P 3")
  SgOps2 = sgtbx.SpaceGroup("P 4")
  assert (SgOps1 == SgOps2) == 0
  assert (SgOps1 != SgOps2) == 1
  SgOps2 = sgtbx.SpaceGroup("P 3")
  SgOps2.makeTidy()
  assert (SgOps1 == SgOps2) == 1
  assert (SgOps1 != SgOps2) == 0

def test_Z2P():
  SgOps = sgtbx.SpaceGroup("F 4 2 3")
  Z2POp = SgOps.getZ2POp()
  print Z2POp.M().as_xyz()
  print Z2POp.InvM().as_xyz()
  P2ZOp = Z2POp.swap()
  print P2ZOp.M().as_xyz()
  print P2ZOp.InvM().as_xyz()
  SgOpsP = SgOps.ChangeBasis(Z2POp)
  SgOpsZ = SgOpsP.ChangeBasis(P2ZOp)
  assert SgOps == SgOpsZ
  x = (0.123, 0.456, 0.789)
  px = Z2POp(x)
  zx = P2ZOp(px)
  print "x (%.6g, %.6g, %.6g)" % x
  print "px (%.6g, %.6g, %.6g)" % px
  print "zx (%.6g, %.6g, %.6g)" % zx
  assert ("(%.6g, %.6g, %.6g)" % x) == ("(%.6g, %.6g, %.6g)" % zx)
  SgOps = sgtbx.SpaceGroup("P 4c -2ab -1b -1ac -1a -1bc")
  Z2POp = SgOps.getZ2POp()
  print Z2POp.M().as_xyz(), Z2POp.InvM().as_xyz()

def test_expand_with_hall():
  SgOps = sgtbx.SpaceGroup("P 2")
  SgOps.ParseHallSymbol("P 2x")
  assert SgOps == sgtbx.SpaceGroup("P 2 2")

def test_check_metrical_matrix():
  uc = uctbx.UnitCell((1, 1, 1, 90, 90, 120))
  G = uc.getMetricalMatrix()
  SgOps = sgtbx.SpaceGroup(sgtbx.parse_string("P 3"))
  SgOps.CheckMetricalMatrix(G)
  SgOps = sgtbx.SpaceGroup(sgtbx.parse_string("P 4 3*"))
  ok = 0
  try:
    SgOps.CheckMetricalMatrix(G)
  except RuntimeError, e:
    ok = 1
  assert ok
  uc = uctbx.UnitCell((1, 1, 1, 90, 90, 90))
  G = uc.getMetricalMatrix()
  SgOps.CheckMetricalMatrix(G)

def test_pickle():
  SgOps = sgtbx.SpaceGroup("P 63")
  import pickle
  pstr = pickle.dumps(SgOps)
  up = pickle.loads(pstr)
  assert SgOps.Info().BuildHallSymbol(1) == up.Info().BuildHallSymbol(1)

def show_position(p):
  print "%.6g %.6g %.6g" % tuple(p)

def show_distance(d):
  print "%.6g" % (d,)

def test_site_symmetry():
  uc = uctbx.UnitCell([])
  sg = sgtbx.SpaceGroup("F 4 2")
  sginfo = sg.Info()
  print sginfo.BuildLookupSymbol()
  WTab = sgtbx.WyckoffTable(sginfo)
  WTab.expand(sg)
  for p in WTab:
    print p.M(), p.Letter(), p.SpecialOp()
    for o in p: print "   ", o
  p = WTab("a")
  print p.M(), p.Letter(), p.SpecialOp()
  X = (-1.45, 2.45, 0.2)
  MinMateDistance = 0.3
  SnapParameters = sgtbx.SpecialPositionSnapParameters(
    uc, sg, 0, MinMateDistance)
  SS = sgtbx.SiteSymmetry(SnapParameters, X, 1)
  show_position(SS.OriginalPosition())
  show_position(SS.SnapPosition())
  show_distance(SS.DistanceMoved())
  show_distance(SS.ShortestDistance())
  print SS.isWellBehaved()
  print SS.M()
  print SS.SpecialOp()
  print SS.PointGroupType()
  for M in SS:
    print M
  SpecialPositionTolerances = sgtbx.SpecialPositionTolerances(
    uc, sg, MinMateDistance, MinMateDistance)
  SE = sgtbx.SymEquivCoordinates(SpecialPositionTolerances, X)
  print SE.M()
  for sx in SE:
    show_position(sx)
  SE = sgtbx.SymEquivCoordinates(SS)
  print SE.M()
  for sx in SE:
    show_position(sx)
  p = WTab.getWyckoffMapping(SS).WP()
  print p.M(), p.Letter(), p.SpecialOp()

def sub_test_random_recycle_xyz(SgOps):
  Decimal = 0
  TrFirst = 0
  LettersXYZ = "xyz"
  Separator = ","
  newSgOps = sgtbx.SpaceGroup()
  matrix_indices = range(SgOps.OrderZ())
  debug_utils.random.shuffle(matrix_indices)
  for i in matrix_indices:
    xyz = SgOps(i).as_xyz(Decimal, TrFirst, LettersXYZ, Separator)
    s = sgtbx.parse_string(xyz)
    rtmx = sgtbx.RTMx(s)
    assert xyz == rtmx.as_xyz(Decimal, TrFirst, LettersXYZ, Separator)
    newSgOps.expandSMx(rtmx)
  assert SgOps.nLTr() == newSgOps.nLTr()
  assert SgOps.fInv() == newSgOps.fInv()
  assert SgOps.nSMx() == newSgOps.nSMx()
  assert SgOps == newSgOps

def test_loop_internal_symbol_table():
  enantiomorphic_pairs = {}
  cbop = sgtbx.ChOfBasisOp(sgtbx.RTMx("x+1/12,y+1/12,z+1/12"))
  iter = sgtbx.SpaceGroupSymbolIterator()
  for s in iter:
    SgOps = sgtbx.SpaceGroup(s)#ChangeBasis(cbop)
    sub_test_random_recycle_xyz(SgOps)
    SgInfo = SgOps.Info()
    print "OrderZ:", SgOps.OrderZ()
    print "isChiral:", SgOps.isChiral()
    ts = SgOps.MatchTabulatedSettings()
    symbols = sgtbx.SpaceGroupSymbols(ts.ExtendedHermann_Mauguin())
    assert symbols.Hall() == s.Hall()
    ch1 = SgInfo.getChangeOfHandOp()
    print s.ExtendedHermann_Mauguin(), ch1.M()
    b = sgtbx.Brick(SgInfo)
    print b
    for i in xrange(3):
      for j in xrange(2):
        p = b(i, j)
        print p.Point(), p.Off()
    e = SgOps.ChangeBasis(ch1)
    eInfo = e.Info()
    if (e != SgOps):
      assert SgInfo.isEnantiomorphic()
      enantiomorphic_pairs[s.SgNumber()] = eInfo.SgNumber()
    else:
      assert not SgInfo.isEnantiomorphic()
    ch2 = eInfo.getChangeOfHandOp()
    ch1ch2 = ch1.M().multiply(ch2.M()).modPositive()
    assert ch1ch2.isUnit() # this is sufficient, but not necessary.
    SECg = sgtbx.SymEquivCoordinates(SgOps, (0.123,0.456,0.789))
    flippedX = ch1(SECg(0))
    SECf = sgtbx.SymEquivCoordinates(e, flippedX)
    Fg = SECg.StructureFactor((3,5,7))
    Ff = SECf.StructureFactor((3,5,7))
    d = abs(Fg) - abs(Ff)
    assert d < 1.e-5
    grid_sg = SgOps.refine_gridding()
    grid_ss = sgtbx.StructureSeminvariant(SgOps).refine_gridding()
    grid_sg_ss = SgOps.refine_gridding(grid_ss)
    eucl = SgInfo.expandAddlGeneratorsOfEuclideanNormalizer(1, 1)
    grid_eucl = eucl.refine_gridding(grid_ss)
    print "Gridding SpaceGroup", grid_sg
    print "Gridding StructureSeminvariant", grid_ss
    print "Gridding SpaceGroup+StructureSeminvariant", grid_sg_ss
    print "Gridding Euclidean normalizer", grid_eucl
  keys = enantiomorphic_pairs.keys()
  keys.sort()
  for p in keys:
    print sgtbx.SpaceGroupSymbols(p).Hermann_Mauguin(),
    print sgtbx.SpaceGroupSymbols(enantiomorphic_pairs[p]).Hermann_Mauguin()
  assert len(enantiomorphic_pairs) == 22

def run():
  import sys, os
  Flags = debug_utils.command_line_options(sys.argv[1:], (
    "RandomSeed",
    "ShortCut",
  ))
  if (not Flags.RandomSeed): debug_utils.set_random_seed(0)
  test_symbols_from_tables()
  test_weird_symbols()
  test_icsd_symbols()
  test_recycle_space_group_symbols()
  test_recycle_space_group_symbols("A1983")
  test_recycle_space_group_symbols("I1952")
  test_rtmx()
  test_cancel_and_tidy()
  test_space_group_comparison()
  test_Z2P()
  test_expand_with_hall()
  test_check_metrical_matrix()
  test_pickle()
  test_site_symmetry()
  test_loop_internal_symbol_table()
  import tst6, tst7, tst9, tst10, tst11, tst12, tst13, tst14
  tst6.run(0)
  tst7.run(0)
  tst9.run(0)
  tst10.run(0)
  tst11.run(0)
  tst12.run(0)
  tst13.run(0)
  tst14.run(0)
  t = os.times()
  print "u+s,u,s:", t[0] + t[1], t[0], t[1]

def run_other(args, timing, function_object, shortcut_settings=None):
  import os
  Flags = debug_utils.command_line_options(args, (
    "RandomSeed",
    "ShortCut",
    "Endless",
    "AllSpaceGroups",
    "AllTabulatedSpaceGroups",
    "AllSettings",
  ))
  if (not Flags.RandomSeed): debug_utils.set_random_seed(0)
  if (Flags.ShortCut):
    settings = shortcut_settings
  elif (Flags.AllSettings):
    from settings import settings # see examples/python/make_settings.py
  elif (Flags.AllTabulatedSpaceGroups):
    settings = []
    for i in sgtbx.SpaceGroupSymbolIterator():
      settings.append(i.ExtendedHermann_Mauguin())
  else:
    settings = debug_utils.get_test_space_group_symbols(Flags.AllSpaceGroups)
  while 1:
    function_object(settings)
    if (timing and not Flags.ShortCut):
      t = os.times()
      print "u+s,u,s:", t[0] + t[1], t[0], t[1]
    if (not Flags.Endless): break

if (__name__ == "__main__"):
  run()
