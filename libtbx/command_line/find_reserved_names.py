from __future__ import absolute_import, division, print_function
import ast
import builtins
import keyword
import os
import sys

reserved_words = keyword.kwlist + [b for b in dir(builtins) if not b.startswith("_")]
reserved_words.remove("copyright")


def check(filename):
    root = ast.parse(open(filename).read())
    text = open(filename).readlines()

    for node in ast.walk(root):
        if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
            if node.id in reserved_words:
                print('%s:%d "%s"' % (filename, node.lineno, node.id))
                print(text[node.lineno - 1])


def main(start):
    for path, dnames, fnames in os.walk(start):
        for f in fnames:
            if f.endswith(".py"):
                check(os.path.join(path, f))


if __name__ == "__main__":
    if len(sys.argv) == 1:
        args = ["."]
    else:
        args = sys.argv[1:]
    for arg in args:
        print("Checking %s" % arg)
        main(arg)
