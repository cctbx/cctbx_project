import sys, os
pythonw = os.path.join(sys.exec_prefix, 'pythonw')
if not os.path.exists(pythonw):
  pythonw = "%s.%s" % (pythonw, 'exe')
  if not os.path.exists(pythonw):
    print "Error: can't find pythonw"
    sys.exit(1)
exit(os.system("%s %s" % (pythonw, " ".join(sys.argv[1:]))))