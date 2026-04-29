# Converted from the NRUBS.java file in the following repository: https://github.com/rlabduke/javadev
# Comments from the original NRUBS.java file.

from __future__ import absolute_import, division, print_function
import numpy as np

'''
<code>NRUBS</code> implements non-rational, uniform B-splines, a much simpler
thing than NURBS (non-uniform rational B-splines).

<p>My reference was Computer Graphics: Principles and Practice, 2nd ed.
in the Systems Programming Series by Foley, van Dam, Feiner, and Hughes (1993).
With some work, I could probably use it to implement NURBS (which are a
generalization of NRUBS, on two counts), but since Prekin uses the plain
B-splines for ribbons, I don't see the point.

<p>Here's my simplified layout of the math of splines. <!-- {{{ -->
B-splines require 4 control points per segment, and adjacent segments
share 3 of the 4 control points.
For simplicity, I deal with only one segment at a time here,
with control points P0, P1, P2, P3.
Segments are C2 continuous at join points, meaning their second derivatives
are equal. Doubling or tripling guide points may break this.
The points along the actual spline are given as Q(t), where t is on [0,1].
Although the Q(t) don't actually pass through any of the control points,
they do lie within the convex hull of those points.
You can make the spline go closer by doubling up a control point, but if
three of them are the same, then that spline segment will be a straight line.

<p><pre>
t                               the parameter for spline coordinates

Q(t) = [x(t) y(t) z(t)]         the functions of t that give spline coordinates

    [P0]   [P0_x P0_y P0_z]
G = [P1] = [P1_x P1_y P1_z]     the control points (guide points)
    [P2]   [P2_x P2_y P2_z]
    [P3]   [P3_x P3_y P3_z]

            [-1  3 -3  1]
M = (1/6) * [ 3 -6  3  0]       the basis matrix, which determines blending functions
            [-3  0  3  0]
            [ 1  4  1  0]

T = [t^3  t^2  t  1]            powers of t -- makes this a CUBIC spline function

Q(t) = T * M * G                the spline function itself
     = T * C                        alternate representation: C = M * G
     = B * G                        alternate representation: B = T * M
</pre>

<p>Notice that if you're going to do a bunch of segments with the same
subdivisions of t (e.g. t = [0.00  0.25  0.50  0.75  1.00]), it makes sense
to precalculate all the B(t) you need. This gives:
<pre>
Q(t) = (1/6) * {  [(1-t)^3]*P0 + [3t^3 - 6t^2 + 4]*P1 + [-3t^3 + 3t^2 + 3t + 1]*P2 + [t^3]*P3  }
</pre>
<!-- }}} -->

<p>Copyright (C) 2006 by Ian W. Davis. All rights reserved.
<br>Begun on Thu Jan 19 16:56:14 EST 2006
'''

class Triple:
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    def setXYZ(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def addMult(self, scalar, other):
        self.x += scalar * other.x
        self.y += scalar * other.y
        self.z += scalar * other.z

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def getZ(self):
        return self.z

class NRUBS:
    def __init__(self):
        self.Bcache = {}
        self.work = Triple()

    def spline(self, guidePts, nIntervals):
        '''Generates interpolated points between two guide points (guideStart and guideStart+1),
        and stores them in splineOut at splineStart through splineStart+nIntervals.
        (Yes, that's right: it calculates nIntervals+1 points.)
        '''
        B = self.getB(nIntervals)
        out = [Triple() for _ in range(nIntervals * (len(guidePts) - 3) + 1)]

        for i in range(len(guidePts) - 3):
            self._spline(guidePts, i, out, i * nIntervals, B)

        return out

    def _spline(self, guidePts, guideStart, splineOut, splineStart, B):
        '''Calculates the N+1 points along the spline for a particular segment.
        Using this repeatedly along a series of guidepoints will result in small
        amount of duplicate calculation (each join calc'd twice).
        '''
        P0 = guidePts[guideStart + 0]
        P1 = guidePts[guideStart + 1]
        P2 = guidePts[guideStart + 2]
        P3 = guidePts[guideStart + 3]

        for i in range(len(B)):
            self.work.setXYZ(0, 0, 0)
            self.work.addMult(B[i][0], P0)
            self.work.addMult(B[i][1], P1)
            self.work.addMult(B[i][2], P2)
            self.work.addMult(B[i][3], P3)
            splineOut[splineStart + i].setXYZ(self.work.getX(), self.work.getY(), self.work.getZ())

    def getB(self, nIntervals):
        '''Returns the "B" coefficients for a particular number N of intervals
        along a spline segment. Because N+1 points are generated (fencepost problem)
        this function returns a double[N+1][4].

        If possible, the value will be retrieved from cache instead of recalculated.
        '''
        if nIntervals not in self.Bcache:
            self.Bcache[nIntervals] = self.calculateB(nIntervals)
        return self.Bcache[nIntervals]

    @staticmethod
    def calculateB(nIntervals):
        '''Returns the "B" coefficients for a particular number N of intervals
        along a spline segment. Because N+1 points are generated (fencepost problem)
        this function returns a double[N+1][4].
        '''
        B = np.zeros((nIntervals + 1, 4))
        for i in range(nIntervals + 1):
            t = i / nIntervals
            t2 = t * t
            t3 = t2 * t
            _1_t = 1 - t
            B[i][0] = (_1_t * _1_t * _1_t) / 6
            B[i][1] = (3 * t3 - 6 * t2 + 4) / 6
            B[i][2] = (-3 * t3 + 3 * t2 + 3 * t + 1) / 6
            B[i][3] = t3 / 6
        return B

# Example usage
if __name__ == "__main__":
    guidePts = [Triple(0, 0, 0), Triple(1, 2, 0), Triple(2, 3, 0), Triple(4, 0, 0)]
    nrubs = NRUBS()
    result = nrubs.spline(guidePts, 10)
    for point in result:
        print("({}, {}, {})".format(point.getX(), point.getY(), point.getZ()))

