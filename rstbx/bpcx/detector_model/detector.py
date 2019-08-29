from __future__ import absolute_import, division, print_function
#!/usr/bin/env python

from scitbx import matrix
from rstbx.bpcx import sensor
from six.moves import range

class detector:
    '''An abstract class definition for X-ray detectors which are assumed to
    be composed of one or more flat rectangular sensor areas. Will initially
    assume that sensor areas to not shadow one another.'''

    def __init__(self, sensors):
        self._sensors = sensors

        return

    def sensors(self):
        return self._sensors

class reflection_predictor:
    '''A class to convert supplied diffraction vectors and detector surfaces
    to impact positions as sensor number and mm in the coordinate frame
    attached to that sensor. Need to return number of sensor and coordinate.
    The units of mm are not enforced - it is whatever the detector matrix d
    is set up with'''

    def __init__(self, rays, sensors):

        self._rays = rays
        self._sensors = sensors
        return

    def intersections(self):
        '''Calculate intersections for each diffracted ray in the list'''

        intersections = []
        for ray in self._rays:
            x = self.intersect_all_sensors(ray)
            if x:
                intersections.append(x)

        return intersections

    def intersect_all_sensors(self, ray):
        # basic way - we have to check the vector against each sensor in turn.
        # there may be more efficient short cuts

        hit = None
        for j, s in enumerate(self._sensors):
            x = self.intersect(ray, s)
            if x:
                if hit is None:
                    dist = (matrix.col(s.origin) + \
                            x[0] * matrix.col(s.dir1) + \
                            x[1] * matrix.col(s.dir2)).length()
                    hit = (j, x[0], x[1], dist)
                else:
                    dist = (matrix.col(s.origin) + \
                            x[0] * matrix.col(s.dir1) + \
                            x[1] * matrix.col(s.dir2)).length()
                    if dist < hit[3]:
                        hit = (j, x[0], x[1], dist)
        return hit

    def intersect(self, ray, sensor):
        '''Compute intersection of sensor with ray from frame origin, returning
        none if intersection not within limits.'''
        raise RuntimeError('overload me')

class reflection_predictor_thomas(reflection_predictor):
    '''Implementation of reflection_predictor using David Thomas' matrix
    formalism'''

    def intersect(self, ray, sensor):
        v = matrix.sqr(sensor.D) * ray
        x1 = v[0] / v[2]
        x2 = v[1] / v[2]

        if x1 < sensor.lim1[0] or x1 > sensor.lim1[1]:
            return None
        if x2 < sensor.lim2[0] or x2 > sensor.lim2[1]:
            return None

        return (x1, x2)

class reflection_predictor_trivial(reflection_predictor):
    '''Implementation of reflection_predictor using the vectors directly'''

    def intersect(self, ray, sensor):
        scale = ray.dot(matrix.col(sensor.normal))
        r = (ray * sensor.distance / scale) - matrix.col(sensor.origin)
        x1 = r.dot(matrix.col(sensor.dir1))
        x2 = r.dot(matrix.col(sensor.dir2))
        if x1 < sensor.lim1[0] or x1 > sensor.lim1[1]:
            return None
        if x2 < sensor.lim2[0] or x2 > sensor.lim2[1]:
            return None

        return (x1, x2)

def predictor_factory(selection):
    '''Factory function to choose a suitable implementation of the
    reflection_predictor class. This is a basic version to flesh out the API.
    A more realistic version might use e.g. PHIL parameters to make the choice
    and would not hard code the alternatives'''

    class_choice = {
        'matrix': reflection_predictor_thomas,
        'vector': reflection_predictor_trivial
    }.get(selection)

    return class_choice

def test_all():
    '''Test all features of the module using a detector composed of four
    overlapping square sensors.'''

    d1 = matrix.col((1, 0, 0))
    d2 = matrix.col((0, 1, 0))
    lim = (0,50)
    panel0 = sensor(matrix.col((-45, -45, 100)), d1, d2, lim, lim)
    panel1 = sensor(matrix.col((-5, -45, 99)), d1, d2, lim, lim)
    panel2 = sensor(matrix.col((-45, -5, 98)), d1, d2, lim, lim)
    panel3 = sensor(matrix.col((-5, -5, 97)), d1, d2, lim, lim)

    d = detector([panel0, panel1, panel2, panel3])
    rays = [matrix.col((0, 0, 1))]

    for ver in ('matrix', 'vector'):
        rp = predictor_factory(ver)
        rp = rp(rays, d.sensors())

        # only one ray to consider
        rpi = rp.intersections()

        # the closest intersection should be on panel3
        assert rpi[0][0] == 3

    print('OK')

def test_work():
    '''Test all features of the module using a detector composed of four
    overlapping square sensors, with lots of rays.'''

    import random, time

    nrays = 10000

    d1 = matrix.col((1, 0, 0))
    d2 = matrix.col((0, 1, 0))
    lim = (0,50)
    panel0 = sensor(matrix.col((-45, -45, 100)), d1, d2, lim, lim)
    panel1 = sensor(matrix.col((-5, -45, 99)), d1, d2, lim, lim)
    panel2 = sensor(matrix.col((-45, -5, 98)), d1, d2, lim, lim)
    panel3 = sensor(matrix.col((-5, -5, 97)), d1, d2, lim, lim)

    d = detector([panel0, panel1, panel2, panel3])

    rays = [matrix.col((random.random() - 0.5,
                        random.random() - 0.5,
                        random.random() - 0.5)) for j in range(nrays)]

    rpi = { }

    for ver in ('matrix', 'vector'):
        t0 = time.time()
        rp = predictor_factory(ver)(rays, d.sensors())

        rpi[ver] = rp.intersections()

    for j, intersection_matrix in enumerate(rpi['matrix']):
        intersection_vector = rpi['vector'][j]
        assert(intersection_matrix[0] == intersection_vector[0])

    print('OK')

if __name__ == '__main__':
    test_work()
