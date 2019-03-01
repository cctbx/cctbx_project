from __future__ import absolute_import, division, print_function
from libtbx import code_analysis
from libtbx.option_parser import option_parser

def run(args, debug=False):
  comments = code_analysis.comments(args, debug=debug)
  print(comments.commented_lines, comments.lines,\
        round(comments.commented_lines/comments.lines * 100, 1))

if __name__ == '__main__':
  import sys
  command_line = (option_parser(
    usage="",
    description="")
    .option(None, "--debug",
            dest='debug',
            action="store_true",
            default=False)
  ).process(args=sys.argv[1:])
  run(command_line.args, **command_line.options.__dict__)
