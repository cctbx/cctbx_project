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
a = shared.bool((0, 1, 0, 1))
b = shared.bool((0, 1, 1, 0))
assert tuple(~a) == (1, 0, 1, 0)
assert tuple(a & b) == (0, 1, 0, 0)
assert tuple(a | b) == (0, 1, 1, 1)
a &= b
assert tuple(a) == (0, 1, 0, 0)
a |= shared.bool((1, 0, 1, 0))
assert tuple(a) == (1, 1, 1, 0)
assert a.count(0) == 1
assert a.count(1) == 3
a &= 1
assert tuple(a) == (1, 1, 1, 0)
a &= 0
assert tuple(a) == (0, 0, 0, 0)
a |= 1
assert tuple(a) == (1, 1, 1, 1)
a |= 0
assert tuple(a) == (1, 1, 1, 1)
a = shared.int((4, 9))
b = shared.int((2, 3))
assert tuple(-a) == (-4, -9)
assert tuple(a + b) == (6, 12)
assert tuple(a - b) == (2, 6)
assert tuple(a * b) == (8, 27)
assert tuple(a / b) == (2, 3)
assert tuple(a % b) == (0, 0)
assert tuple(a.add(3)) == (7, 12)
assert tuple(a.sub(4)) == (0, 5)
assert tuple(a.mul(5)) == (20, 45)
assert tuple(a.div(2)) == (2, 4)
assert tuple(a.mod(2)) == (0, 1)
assert shared.sum(a) == 13
a += 3
assert tuple(a) == (7, 12)
a -= 3
assert tuple(a) == (4, 9)
a *= 3
assert tuple(a) == (12, 27)
a /= 3
assert tuple(a) == (4, 9)
a %= 3
assert tuple(a) == (1, 0)
a = shared.int((4, 9))
assert tuple(a.as_double()) == (4, 9)
b = shared.int((2, 12))
assert tuple(a == b) == (0, 0)
assert tuple(a != b) == (1, 1)
assert tuple(a < b) == (0, 1)
assert tuple(a > b) == (1, 0)
assert tuple(a <= b) == (0, 1)
assert tuple(a >= b) == (1, 0)
assert tuple(a == 9) == (0, 1)
a = shared.int((-1, 0, 1))
assert tuple(shared.abs(a)) == (1, 0, 1)
a = shared.double((-2, 0, 3))
assert tuple(shared.pow(a, 2)) == (4, 0, 9)
a = shared.double((1, 0, 3, 2))
b = shared.double((4, 5, 6))
assert shared.min_index(a) == 1
assert shared.min_index(b) == 0
assert shared.max_index(a) == 2
assert shared.max_index(b) == 2
assert shared.min(a) == 0
assert shared.min(b) == 4
assert shared.max(a) == 3
assert shared.max(b) == 6
assert shared.sum(a) == 6
assert shared.sum(b) == 15
assert shared.product(a) == 0
assert shared.product(b) == 120
assert shared.mean(a) == 6. / 4
assert shared.mean(b) == 15. / 3
assert shared.rms(a) == math.sqrt((1+3.*3.+2.*2.) / 4)
assert shared.rms(b) == math.sqrt((4.*4.+5.*5.+6.*6.) / 3)
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
import pickle
a = shared.bool((1,0,1))
p = pickle.dumps(a)
b = pickle.loads(p)
assert b.size() == 3
assert tuple(b) == (1,0,1)
a = shared.double((1,2,3))
a = shared.double(())
p = pickle.dumps(a)
b = pickle.loads(p)
assert b.size() == 0
a = shared.double((1,2,3))
p = pickle.dumps(a)
b = pickle.loads(p)
assert b.size() == 3
assert tuple(b) == (1,2,3)
a = shared.int((1,2,3))
p = pickle.dumps(a)
b = pickle.loads(p)
assert b.size() == 3
assert tuple(b) == (1,2,3)
a = shared.float((1,2,3))
p = pickle.dumps(a)
b = pickle.loads(p)
assert b.size() == 3
assert tuple(b) == (1,2,3)
a = shared.complex_double((1+2j, 2+3j, 4+5j))
p = pickle.dumps(a)
b = pickle.loads(p)
assert tuple(b) == (1+2j, 2+3j, 4+5j)
assert b.size() == 3
print "OK"
