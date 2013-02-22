def plot_shoebox(shoebox):
    '''
    Render a shoebox as a sequence of png images using matplotlib.
    '''

    frame, slow, fast = shoebox.accessor().all()

    max_value = 100

    import wx
    import warnings
    warnings.filterwarnings("ignore", category=wx.wxPyDeprecationWarning)

    for j in range(frame):
        from scitbx.array_family import flex
        import numpy
        import matplotlib
        from matplotlib import pyplot
        grid = flex.grid((slow, fast))
        slice = flex.double(grid)
        for s in range(slow):
            for f in range(fast):
                if shoebox[(j, s, f)] < 0:
                    slice[(s, f)] = max_value
                else:
                    slice[(s, f)] = shoebox[(j, s, f)]

        image = numpy.reshape(slice.as_numpy_array(), (slow, fast))
        pyplot.imshow(image, cmap = matplotlib.cm.Greys, vmax = max_value,
                      interpolation = 'nearest', norm = None, filternorm = 0.0)
        pyplot.savefig('frame%02d.png' % (j))

def tst_cube():

    import os
    import sys

    from dxtbx.model.cube import cube

    import os

    template = os.path.join('/Users/graeme/svn/cctbx/sources/cctbx_project/bpcx_regression/use_case_xds_method',
                            'thau2_O0_K0_P0_1_####.cbf')

    c = cube(template)

    s = c.get(1, 5, 0, 100, 0, 100)

    # plot_shoebox(s)

    print 'OK'

if __name__ == '__main__':
    tst_cube()
