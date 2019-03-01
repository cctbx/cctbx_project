# LIBTBX_SET_DISPATCHER_NAME dev.dxtbx.debug_memory
from __future__ import absolute_import, division, print_function

import resource
import sys

import dxtbx


def main():
    frame = sys.argv[1]
    powers = [2 ** n for n in range(20)]

    for j in range(powers[-1] + 1):
        i = dxtbx.load(frame)
        mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if j in powers:
            print(j, mem)


if __name__ == "__main__":
    main()
