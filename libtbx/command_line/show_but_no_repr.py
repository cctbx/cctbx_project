from __future__ import division, print_function
import sys
import ast

def class_has_show_but_not_repr(node):
    found_repr = False
    found_show = False

    for obj in node.body:
        if isinstance(obj, ast.FunctionDef):
            if obj.name == '__repr__':
                found_repr = True
            if obj.name == '__str__':
                found_repr = True
            if obj.name == 'show':
                found_show = True
    return found_show and not found_repr

class rummage(ast.NodeVisitor):
    def __init__(self, node):
        self.__offenders = []
        self.visit(node)

    def visit_ClassDef(self, node):
        if class_has_show_but_not_repr(node):
            self.__offenders.append(node.name)
        self.generic_visit(node)

    def offenders(self):
        return self.__offenders

if __name__ == '__main__':
    for arg in sys.argv:
        with open(arg) as f:
            tree = ast.parse(f.read())

        r = rummage(tree)
        off = r.offenders()
        if off:
            print('File: %s' % arg)
            for o in off:
                print('Class: %s' % o)
