from __future__ import absolute_import, division, print_function
from dxtbx.format.Registry import Registry

Registry.setup()
import sys


def print_class(classobj, filename=None, depth=1):
    if filename is None:
        print("% 5d" % depth, "  " * depth, classobj.__name__)
    else:
        try:
            ok = classobj.understand(filename)
        except Exception:
            ok = False
        if ok:
            print("% 5d" % depth, "  " * depth, classobj.__name__)
        else:
            return
    for child in classobj._children:
        print_class(child, filename, depth + 1)


def show_registry(filename=None):
    if filename is None:
        extrabit = ""
    else:
        extrabit = " that understand image %s" % filename
    print(
        "Showing hierarchy of classes in the dxtbx registry%s. The root classes are shown with depth of 1, and subclasses are shown indented and with a higher depth number."
        % extrabit
    )
    print()
    print("Depth  Class name")
    print("    0  Format")

    for classobj in Registry._formats:
        print_class(classobj, filename)


if __name__ == "__main__":
    if len(sys.argv) == 2:
        show_registry(sys.argv[1])
    else:
        show_registry()
