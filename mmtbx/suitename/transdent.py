from __future__ import nested_scopes, generators, division, absolute_import
from __future__ import  with_statement, print_function
import re, os, sys

N = 2   # indentation standard
white = "([ \t]*)"
simple = re.compile(".*(#|:)")
stringy = re.compile("[^\"#:]*(\".*\"|'.*'|#|:[ \t]*(#.*)?$)")

class Line(object):
    content = ""
    level = 0
    actual = 0  # the actual whitespace at the front


def analyzeWhitespace(white):
    count = 0
    for char in white:
        if char == "\t":
            count += N
        else:
            count += 1
    return count


def analyzeLine(raw):
    # discover indentation of raw line and whether it ends with a colon
    content = raw.strip()
    match1 = re.match(white, raw)
    mark = match1.end()
    spaces = analyzeWhitespace(match1[1])
    if content == "":  # blank line
        spaces = 0
    while True:
        # this loop is exited only by returning
        match2 = stringy.match(raw, mark)
        if match2 is None:
            return content, spaces, False

        what = match2.group(1)
        char = what[0]
        if char == "#":
            return content, spaces, False
        elif char == ":":
            return content, spaces, True
        elif char == '"' or char == "'":
            mark = match2.end()


def main():
    file = sys.argv[1] + ".py"
    #root, ext = os.path.splitext(file)
    outFile = sys.argv[1] + ".ind.py"
    # file = "target.py"
    stream = open(file)
    lines = stream.readlines()
    oStream = open(outFile, "w")
    #oStream = sys.stderr

    indents = []
    prevSpaces = 0
    level = 0
    newLevel = False
    freeze = False
    count=-1
    for rawLine in lines:
        count += 1
        content, spaces, indenting = analyzeLine(rawLine)
        if spaces == 0 and rawLine.startswith("def "):
          # sanity reset: no matter what, restart at left margin
          # when defining a top level function
          level = 0; newLevel = False
        elif newLevel:
          if spaces > prevSpaces:
            # A true syntactic indent (else a stray colon, ignore)
            level += 1
            indents.append(prevSpaces)
          newLevel = False
        elif spaces > prevSpaces:
          # indent without : is a line-continuation or data-formatting situation
          freeze = True
          indents.append(prevSpaces)
        if content == "":  # blank lines influence nothing
          print("", file=oStream)
          continue
#        if spaces < prevSpaces and level > 0:
        if spaces < prevSpaces and spaces <= indents[-1]:
            oldSpaces = indents.pop()
            while spaces < oldSpaces:
                oldSpaces = indents.pop()
            level = len(indents)
            freeze = False
        prevSpaces = spaces
        if freeze:
          print((level * N + (spaces-indents[-1])) * " " + content, file=oStream)
        else:
          print(level * N * " " + content, file=oStream)
        if indenting:
            newLevel = True
    oStream.close()

main()

