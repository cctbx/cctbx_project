from scitbx import matrix

def vectors(handle,all):
  from scitbx.array_family import flex
  x = flex.double()
  y = flex.double()
  xpred = flex.double()
  ypred = flex.double()
  for line in handle.readlines():
    if line.find("CV ")!=0: continue
    tokens = line.split()
    obscen = matrix.col((float(tokens[2]),float(tokens[3])))
    refcen = matrix.col((float(tokens[5]),float(tokens[6])))
    obsspo = matrix.col((float(tokens[8]),float(tokens[9])))
    predspo = matrix.col((float(tokens[11]),float(tokens[12])))
    xpred.append(predspo[0])
    ypred.append(predspo[1])
    prediction = predspo-refcen
    observation = obsspo-obscen
    cv = prediction-observation
    x.append(cv[0])
    y.append(cv[1])

  if all:
    print "Plotting all %d spots one graph."%len(x)
    from matplotlib import pyplot as plt
    plt.plot(x,y,"r.")
    plt.show()
    return

  print len(x),len(y)
  #from spotfinder.applications.xfel.cxi_run3 import get_initial_cxi_scope
  #print "ONLY FOR RUN 3!!!"
  #params = get_initial_cxi_scope()
  #tiling = params.distl.detector_tiling

  print "ONLY FOR RUN 4!!!"
  tiling = [518, 439, 712, 624, 715, 439, 909, 624, 519, 652, 713, 837, 716, 652, 910, 837, 510, 19, 695, 213, 510, 216, 695, 410, 721, 19, 906, 213, 721, 216, 906, 410, 87, 233, 281, 418, 284, 233, 478, 418, 88, 20, 282, 205, 285, 20, 479, 205, 108, 447, 293, 641, 108, 644, 293, 838, 321, 445, 506, 639, 321, 642, 506, 836, 437, 853, 622, 1047, 437, 1050, 622, 1244, 649, 853, 834, 1047, 649, 1050, 834, 1244, 19, 1069, 213, 1254, 216, 1069, 410, 1254, 18, 856, 212, 1041, 215, 856, 409, 1041, 230, 1282, 415, 1476, 230, 1479, 415, 1673, 16, 1282, 201, 1476, 16, 1479, 201, 1673, 442, 1469, 636, 1654, 639, 1469, 833, 1654, 443, 1257, 637, 1442, 640, 1257, 834, 1442, 852, 1137, 1046, 1322, 1049, 1137, 1243, 1322, 852, 925, 1046, 1110, 1049, 925, 1243, 1110, 1067, 1350, 1252, 1544, 1067, 1547, 1252, 1741, 854, 1352, 1039, 1546, 854, 1549, 1039, 1743, 1280, 1342, 1474, 1527, 1477, 1342, 1671, 1527, 1282, 1554, 1476, 1739, 1479, 1554, 1673, 1739, 1467, 924, 1652, 1118, 1467, 1121, 1652, 1315, 1255, 925, 1440, 1119, 1255, 1122, 1440, 1316, 1142, 521, 1327, 715, 1142, 718, 1327, 912, 930, 521, 1115, 715, 930, 718, 1115, 912, 1359, 514, 1553, 699, 1556, 514, 1750, 699, 1358, 727, 1552, 912, 1555, 727, 1749, 912, 1353, 92, 1538, 286, 1353, 289, 1538, 483, 1565, 91, 1750, 285, 1565, 288, 1750, 482, 932, 111, 1126, 296, 1129, 111, 1323, 296, 931, 323, 1125, 508, 1128, 323, 1322, 508]

  for itile in xrange(len(tiling)//4):
    print "tile",itile,
    print "(%4d %4d)-(%4d %4d)"%tuple(tiling[4*itile:4*itile+4]),
    selection = flex.bool()
    for i in xrange(len(x)):
      selection.append(
         tiling[4*itile+0]<xpred[i]<tiling[4*itile+2] and
         tiling[4*itile+1]<ypred[i]<tiling[4*itile+3]
     )
    print "in selection of %d/%d"%(selection.count(True),len(selection))
    if selection.count(True)<10:continue
    from matplotlib import pyplot as plt
    plt.plot(x.select(selection),y.select(selection),"r.")
    plt.plot([flex.mean(x.select(selection))],[flex.mean(y.select(selection))],"go")
    print "Delta x=%.1f"%flex.mean(x.select(selection)),"Delta y=%.1f"%flex.mean(y.select(selection))
    plt.show()

if __name__=="__main__":
  import sys
  handle = open(sys.argv[1],"r")
  vectors(handle,all = "all" in sys.argv)
