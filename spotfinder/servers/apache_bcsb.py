import StringIO

def run(args, verbose=False):
  from libtbx.utils import Sorry
  import os
  from spotfinder.command_line.signal_strength import master_params

  #For the Apache server version, do not allow site, user, or dataset preferences
  #all parameters are to be passed in through the http: query line

  phil_objects = []
  argument_interpreter = master_params.command_line_argument_interpreter(
    home_scope="distl")

  from spotfinder.servers.apache_utils import LongLineSimpleNode as SimpleNode
  from spotfinder.applications import signal_strength

  logger = StringIO.StringIO()
  top = SimpleNode("spotfinder")

  try:
    for key in args.keys():
        arg = "%s=%s"%(key,args.get(key,""))
        command_line_params = argument_interpreter.process(arg=arg)
        phil_objects.append(command_line_params)

    working_params = master_params.fetch(sources=phil_objects)
    params = working_params.extract()

    top.child(SimpleNode(tag="file_name",contents=params.distl.image))

    if not os.path.isfile(params.distl.image):
      raise Sorry("%s not a readable file"%params.distl.image)

    Org = signal_strength.run_signal_strength(params)
    assert len(Org.S.images.keys())==1 # there is only one image
    key = Org.S.images.keys()[0]

    # List of spots between specified high- and low-resolution limits
    if Org.S.images[key].has_key('lo_pass_resolution_spots'):
      spots = Org.S.images[key]['lo_pass_resolution_spots']
    elif Org.S.images[key].has_key('inlier_spots'):
      spots = Org.S.images[key]['inlier_spots']
    else:
      spots = []

    if Org.S.images[key].has_key('N_spots_total'):
      total = "%d"%Org.S.images[key]["N_spots_total"]
    else:
      total = "0"

    if Org.S.images[key].has_key('resolution'):
      resolution = Org.S.images[key]['resolution']
    elif Org.S.images[key].has_key('distl_resolution'):
      resolution = Org.S.images[key]['distl_resolution']
    else:
      resolution = 0.0

    top.child(SimpleNode(tag="total_spots",contents=total))
    top.child(SimpleNode(tag="good_spots",contents="%d"%len(spots)))
    top.child(SimpleNode(tag="resolution",contents="%.3f"%resolution))

    if len(Org.S.reporters[key])==0:
      top.child(SimpleNode(tag="total_integrated",contents="0"))
      top.child(SimpleNode(tag="mean_isigi",contents="0"))
    else:
      reporter = Org.S.reporters[key][-1]
      normalizer = reporter.weights.sum()
      summation = 0;
      for x in xrange(reporter.S_table_rows):
        summation += reporter.weights[x] * reporter.MeanIsigI[x]
      top.child(SimpleNode(tag="mean_isigi",contents="%.3f"%(summation/normalizer)))
      integrated = reporter.Integrated.sum()
      top.child(SimpleNode(tag="integrated",contents="%.3f"%integrated))

  except Exception,e:
    top.child(SimpleNode(tag="status",contents=repr(e)))
    top.emit(logger)
    return logger.getvalue()

  top.child(SimpleNode(tag="status",contents="OK"))
  top.emit(logger)

  return logger.getvalue()
