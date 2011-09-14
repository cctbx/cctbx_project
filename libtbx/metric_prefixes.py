"""http://en.wikipedia.org/wiki/SI_prefix

An SI prefix (also known as a metric prefix) is a name or associated
symbol that precedes a basic unit of measure (or its symbol) to form a
decimal multiple or submultiple.
"""

yotta = 1e24
zetta = 1e21
exa = 1e18
peta = 1e15
tera = 1e12
giga = 1e9
mega = 1e6
kilo = 1e3
hecto = 1e2
deca = 1e1
deci = 1e-1
centi = 1e-2
milli = 1e-3
micro = 1e-6
nano = 1e-9
pico = 1e-12
femto = 1e-15
atto = 1e-18
zepto = 1e-21
yocto = 1e-24

if (__name__ == "__main__"):
  p = yotta*zetta*exa*peta*tera*giga*mega*kilo*hecto*deca \
      *deci*centi*milli*micro*nano*pico*femto*atto*zepto*yocto
  assert abs(p-1) < 1.e-15
  print "OK"
