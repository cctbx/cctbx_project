def gss(f,a,b,eps=1e-7,N=1000):
  c = 0.61803398875 #golden ratio
  x1 = c*a + (1-c)*b
  fx1 = f(x1)
  x2 = (1-c)*a + c*b
  fx2 = f(x2)
  for i in range(1,N-2):
     if fx1 < fx2:
        b = x2
        x2 = x1
        fx2 = fx1
        x1 = c*a + (1-c)*b
        fx1 = f(x1)
     else:
        a = x1
        x1 = x2
        fx1 = fx2
        x2 = (1-c)*a + c*b
        fx2 = f(x2)
     if (abs(b-a) < eps):
       if fx1 < fx2:
         return x1
       else:
        return x2
  if fx1 < fx2:
    return x1
  else:
    return x2

def function(x):
  y = (x-3.5)**2.0
  return y

def tst_it():
  opti = gss(function,-10,10)
  assert abs(opti-3.5)<1e-4

if __name__ == "__main__":
  tst_it()
  print "OK"
