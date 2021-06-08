import re, os, sys

white = "(\s*)"
simple = re.compile(".*(#|:)")
stringy = re.compile(".*(\".*\"|'.*'|#|:)")


class Line:
    content = ""
    level = 0
    actual = 0  # the actual whitespace at the front


def analyzeLine(raw):
    # discover indentation of raw line and whether it ends with a colon
    content = raw.strip()
    match1 = re.match(white, raw)
    mark = match1.end()
    spaces = match1[1]
    while True:
        # this loop is exited only by returning
        match2 = stringy.match(raw, mark)
        what = match2.group(1)
        if what is None:
            return content, spaces, False

        char = what[0]
        if char == "#":
            return content, spaces, False
        elif char == ":":
            return content, spaces, True
        elif char == '"' or char == "'":
            mark = match2.end()


def main():
    file = sys.argv[1]
    stream = open(file)
    lines = stream.readlines()

    prevSpaces = 0
    level = 0
    newLevel = True
    for rawLine in lines:
        if newLevel:
            # Determined by PREVIOUS line
            level += 1
            newLevel = False
        content, spaces, indenting = analyzeLine(rawLine)
        if spaces < prevSpaces and level > 0:
            level -= 1
            prevSpaces = spaces
        pass
        print(level * N * " " + content)
        if indenting:
            newLevel = level + 1
