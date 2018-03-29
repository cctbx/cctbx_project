from __future__ import absolute_import, division
from dxtbx.format.Registry import Registry
Registry.setup()

print "Showing hierarchy of classes in the dxtbx registry. The root classes are shown with depth of 1, and subclasses are shown indented and with a higher depth number."
print
print "Depth  Class name"
print "    0  Format"

def print_class(classobj, depth = 1):
    print "% 5d"%depth, "  " * depth, classobj.__name__
    for child in classobj._children:
        print_class(child, depth + 1)

for classobj in Registry._formats: print_class(classobj)
