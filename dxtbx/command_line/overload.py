from __future__ import division

def overload(image_file):
    from dxtbx import load
    i = load(image_file)
    d = i.get_detector()
    return max(i.get_raw_data()) > d.get_trusted_range()[1]

if __name__ == '__main__':
    import sys

    for image_file in sys.argv[1:]:
        if overload(image_file):
            print image_file

