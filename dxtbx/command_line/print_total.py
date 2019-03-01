from __future__ import absolute_import, division, print_function


def print_total():
    import sys
    from dxtbx.format.Registry import Registry

    # this will do the lookup for every frame - this is strictly not needed
    # if all frames are from the same instrument

    for arg in sys.argv[1:]:
        print("=== %s ===" % arg)
        format_class = Registry.find(arg)
        print("Using header reader: %s" % format_class.__name__)
        i = format_class(arg)
        image_size = i.get_detector()[0].get_image_size()
        d = i.get_raw_data()
        if not isinstance(d, tuple):
            d = (d,)
        d = [p.as_1d() for p in d]
        total = sum([sum(p.select(p >= 0)) for p in d])
        print("Total Counts: %d" % total)
        print("Average Counts: %.2f" % (total / (image_size[0] * image_size[1])))


if __name__ == "__main__":
    print_total()
