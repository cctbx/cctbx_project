import fftbx

def initcseq(N):
  cseq = fftbx.vector_of_double()
  for i in xrange(1, N + 1):
    cseq.append(float(37-i)/(N+11))
    cseq.append(float(i-73)/(N+13))
  return cseq

def showcseq(cseq):
  for i in xrange(len(cseq) / 2):
    print "%12.4f %12.4f" % (cseq[2*i], cseq[2*i+1])
  sys.stdout.flush()

def initrseq(N):
  rseq = fftbx.vector_of_double()
  for i in xrange(1, N + 1):
    rseq.append(float(37-i)/(N+11))
  return rseq

def showrseq(rseq):
  for i in xrange(len(rseq)):
    print "%12.4f" % (rseq[i],)
  sys.stdout.flush()

def showfactors(xfft):
  print N,
  for f in xfft.Factors(): print f,
  print
  #for w in xfft.WA(): print "%8.5f" % (w,)
  sys.stdout.flush()

def test_complex(Nmax, transform_type):
  for N in xrange(2, Nmax + 1):
    print "N %4d" % (N,)
    cseq = initcseq(N)
    cfft = fftbx.complex_to_complex(N)
    if (transform_type == "f"):
      cfft.forward(cseq)
    else:
      cfft.backward(cseq)
    showcseq(cseq)

def test_real(Nmax, transform_type):
  for N in xrange(2, Nmax + 1):
    print "N %4d" % (N,)
    rseq = initrseq(N)
    rfft = fftbx.real_to_complex(N)
    if (transform_type == "f"):
      rfft.forward(rseq)
    else:
      rfft.backward(rseq)
    showrseq(rseq)

def trace_func(frame, event, arg):
  #if (event == "call"):
  print frame.f_code.co_filename, frame.f_lineno
  return trace_func

if (__name__ == "__main__"):
  import sys
  #sys.settrace(trace_func)
  assert sys.argv[1][0] in "cr"
  assert sys.argv[1][1] in "bf"
  Nmax = float(sys.argv[2])
  if (sys.argv[1][0] == "c"):
    test_complex(Nmax, sys.argv[1][1])
  else:
    test_real(Nmax, sys.argv[1][1])
