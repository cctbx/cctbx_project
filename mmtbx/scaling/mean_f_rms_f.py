from __future__ import absolute_import, division, print_function
import sys,math

def mean_fh(d_min,b):
  if b <1.0: b=1.0
  pi=3.14159
  return 6.0 * pi**0.5 * d_min**3 * math.erf(b**0.5/(2.*d_min))/b**1.5 \
    - 6.0 * d_min**2 * math.exp(-1.0*b/(4.0*d_min**2))/b
def mean_fh_square(d_min,b):
  if b <1.0: b=1.0
  pi=3.14159
  return 3.0 * (pi/2.)**0.5 * d_min**3 * math.erf((b/2.)**0.5/d_min)/b**1.5 \
    - 3.0 * d_min**2 * math.exp(-1.0*b/(2.0*d_min**2))/b

def ratio_mean_f_to_rms_f(d_min,b):
  fh=mean_fh(d_min,b)
  fh2=mean_fh_square(d_min,b)
  return fh/fh2**0.5

def test(out=sys.stdout):
  b=90
  d_min=3.5
  ratio=ratio_mean_f_to_rms_f(d_min,b)
  text= "B: %7.2f   D_min: %7.2f  <f>/<f**2>**0.5: %7.3f" %( b,d_min,ratio)
  print(text, file=out)
  expected_text="B:   90.00   D_min:    3.50  <f>/<f**2>**0.5:   0.889"
  o=float(text.split()[-1])
  e=float(expected_text.split()[-1])
  if abs(o-e) > 0.02:
    raise AssertionError("test of mean_f_rms_f failed: Expected %7.3f got %7.3f: diff is %7.3f" %(
      o,e,abs(o-e)))

def exercise_mean_fh():
  i=-1
  for b in [0,30,60,90,200]:
    for d_min in [2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8]:
       i+=1
       ratio=ratio_mean_f_to_rms_f(d_min,b)
       print("B: %7.2f   D_min: %7.2f  <f>/<f**2>**0.5: %7.3f" %(
         b,d_min,ratio))

if __name__=="__main__":
  if 'exercise' in sys.argv:
    exercise_mean_fh()
    sys.exit(0)

  if 'tst' in sys.argv:
    test()
    sys.exit(0)

  sum_n=0.
  sum_f=0.
  sum_f2=0.
  d_min=2.5
  b=float(sys.argv[1].split("_")[2])
  fh=mean_fh(d_min,b)
  fh2=mean_fh_square(d_min,b)
  ratio=fh/math.sqrt(fh2)

  import math

  for line in open(sys.argv[1]).readlines():
    spl=line.lstrip().rstrip().split()
    if len(spl) != 5: continue
    f2=float(spl[3])
    f=math.sqrt(f2)
    sum_n+=1
    sum_f+=f
    sum_f2+=f2
  print("Mean f  %7.2f  rms f: %7.2f  Mean/rms: %7.2f  N: %d  Target: %7.2f" %(
     sum_f/sum_n, math.sqrt(sum_f2/sum_n),(sum_f/sum_n)/math.sqrt(sum_f2/sum_n),int(sum_n),ratio))
