import math
from cctbx_boost.arraytbx import shared
a = shared.int(3)
assert a.size() == 3
assert tuple(a) == (0, 0, 0)
a = shared.int()
assert a.size() == 0
a.assign(3, 1)
assert tuple(a) == (1, 1, 1)
a.push_back(2)
assert tuple(a) == (1, 1, 1, 2)
a.insert(0, 3)
assert tuple(a) == (3, 1, 1, 1, 2)
a.insert(2, 3, 4)
assert tuple(a) == (3, 1, 4, 4, 4, 1, 1, 2)
a.pop_back()
assert tuple(a) == (3, 1, 4, 4, 4, 1, 1)
a.erase(2)
assert tuple(a) == (3, 1, 4, 4, 1, 1)
a.erase(3, 5)
assert tuple(a) == (3, 1, 4, 1)
a.clear()
assert a.size() == 0
a = shared.double((0, 1, 2, 3))
b = shared.double((4, 5, 6))
a.append(b)
assert tuple(a) == (0, 1, 2, 3, 4, 5, 6)
assert tuple(a.indices()) == tuple(xrange(a.size()))
assert list(a.items()) == zip(xrange(a.size()), a)
x = shared.double((0, 1, 2, 3))
y = shared.double((1, 3, 5, 7))
r = shared.linear_regression(x, y)
assert r.is_well_defined()
assert abs(r.b() - 1) < 1.e-6
assert abs(r.m() - 2) < 1.e-6
assert abs(r.cc() - 1) < 1.e-6
s = shared.statistics(x)
assert (s.min() - 0) < 1.e-6
assert (s.max() - 3) < 1.e-6
assert (s.mean() - 6./4.) < 1.e-6
assert (s.mean2() - 14./4.) < 1.e-6
assert (s.sigma() - math.sqrt(14./4. - 36./16.)) < 1.e-6
c = 1+2j
x = shared.complex_double((c,))
a = shared.abs(x)
assert abs(a[0] - abs(c)) < 1.e-6
p = shared.arg(x)
y = shared.polar(a, p)
d = y[0]
assert abs(d.real - c.real) < 1.e-6
assert abs(d.imag - c.imag) < 1.e-6
p = shared.arg(x, 0)
y = shared.polar(a, p, 0)
d = y[0]
assert abs(d.real - c.real) < 1.e-6
assert abs(d.imag - c.imag) < 1.e-6
p = shared.arg(x, 1)
y = shared.polar(a, p, 1)
d = y[0]
assert abs(d.real - c.real) < 1.e-6
assert abs(d.imag - c.imag) < 1.e-6
y = shared.polar(a, p, 0)
d = y[0]
assert abs(d.real - c.real) > 1.e-6
assert abs(d.imag - c.imag) > 1.e-6
print "OK"
