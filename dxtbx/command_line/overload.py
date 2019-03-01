from __future__ import absolute_import, division, print_function


def overload(image_file):
    from dxtbx import load

    i = load(image_file)
    data = i.get_raw_data()
    if not isinstance(data, tuple):
        data = (data,)
    detector = i.get_detector()
    for pid, (d, p) in enumerate(zip(data, detector)):
        if max(d) > p.get_trusted_range()[1]:
            return True
    return False


if __name__ == "__main__":
    import sys

    for image_file in sys.argv[1:]:
        if overload(image_file):
            print(image_file)
