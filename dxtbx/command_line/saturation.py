from __future__ import division

def saturation(image_file):
    from dxtbx import load
    i = load(image_file)
    d = i.get_detector()
    return i.get_scan().get_image_range()[0], \
           max(i.get_raw_data()) / d.get_trusted_range()[1]

if __name__ == '__main__':
    import sys

    for image_file in sys.argv[1:]:
        i, s = saturation(image_file)
        print '%6d %.6f' % (i, s)

