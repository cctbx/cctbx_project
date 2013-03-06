from __future__ import division
import re

# filename template code stolen from xia2...

# N.B. these are reversed patterns...

patterns = [r'([0-9]{2,12})\.(.*)',
            r'(.*)\.([0-9]{2,12})_(.*)',
            r'(.*)\.([0-9]{2,12})(.*)']

joiners = ['.', '_', '']

compiled_patterns = [re.compile(pattern) for pattern in patterns]

def template_regex(filename):
    '''Try a bunch of templates to work out the most sensible. N.B. assumes
    that the image index will be the last digits found in the file name.'''

    rfilename = filename[::-1]

    global patterns, compiled_patterns

    for j, cp in enumerate(compiled_patterns):
        match = cp.match(rfilename)
        if not match:
            continue
        groups = match.groups()

        if len(groups) == 3:
            exten = '.' + groups[0][::-1]
            digits = groups[1][::-1]
            prefix = groups[2][::-1] + joiners[j]
        else:
            exten = ''
            digits = groups[0][::-1]
            prefix = groups[1][::-1] + joiners[j]

        template = prefix + ''.join(['#' for d in digits]) + exten
        break

    return template, int(digits)
