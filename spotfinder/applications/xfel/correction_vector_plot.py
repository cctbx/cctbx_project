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
  from spotfinder.applications.xfel.cxi_run3 import get_initial_cxi_scope
  print "ONLY FOR RUN 3!!!"
  params = get_initial_cxi_scope()
  tiling = params.distl.detector_tiling
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
    plt.show()

if __name__=="__main__":
  import sys
  handle = open(sys.argv[1],"r")
  vectors(handle,all = "all" in sys.argv)
