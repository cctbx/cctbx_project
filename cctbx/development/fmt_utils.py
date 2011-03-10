def dot6g_list(l):
  s = ""
  for x in l: s += (" %.6g" % (x,))
  return s[1:]

def dot6gdot(x):
  s = ("%.6g" % (x,)).upper()
  if ("E" in s): return "%.10f" % (x,)
  if ("." in s): return s
  return s + "."

def dot6gdot_list(l):
  s = ""
  for x in l: s += " " + dot6gdot(x)
  return s[1:]

def dot6fdot(x):
  s = "%.6f" % (x,)
  if ("." in s): return s
  return "." + s

def dot6fdot_list(l):
  s = ""
  for x in l: s += " " + dot6fdot(x)
  return s[1:]

def xtal6g(x):
  import string
  s = "%.6g" % (x,)
  i = string.find(s, "e")
  if (i < 0): return s
  return s[:i] + s[i+1:]

def xtal6g_list(l):
  s = ""
  for x in l: s += " " + xtal6g(x)
  return s[1:]

def MathFormVec3(V):
  return ("{%.6g,%.6g,%.6g}" % tuple(V)).replace("e", "*^")

def MathFormMx33(M):
  return ("{{%.6g,%.6g,%.6g},{%.6g,%.6g,%.6g},{%.6g,%.6g,%.6g}}"
          % tuple(M)).replace("e", "*^")

def MathForm_adp(adp):
  return MathFormMx33((adp[0], adp[3], adp[4],
                       adp[3], adp[1], adp[5],
                       adp[4], adp[5], adp[2]))
