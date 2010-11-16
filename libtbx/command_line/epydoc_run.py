import epydoc
import epydoc.docwriter
import epydoc.docwriter.xlink
import epydoc.cli

help_text = """Synopsis

libtbx.epydoc [options]

Description

Generate the epydoc documentation specified by the given configuration file.
Any option understood by epydoc may be passed to this script.
If the configuration file is not specified, the file epydoc.conf in the current directory is used.

Requirements

Epydoc 3.0 and docutils.

Instructions for documentation writers

When a pure Python class A uses the Boost.Python wrapping of a C++ class B,
the docstring of A should feature a link to the doxygen-generated
documentation of B. That link shall be written as e.g.
  :doxyclass:`scitbx::lbfgs::drop_convergence_test`
and will give
<a href="../c_plus_plus/classscitbx_1_1lbfgs_1_1drop__convergence__test.html">
class scitbx::lbfgs::drop_convergence_test</a>
The other magic keyword is "doxystruct" for struct instead of class.

This relies on an essential requirement on the part of the C++ documentation writers: the C++ documentation root directory and the Python documentation
root directory shall be in the same directory on the server.
"""

def help():
  print help_text
  exit(1)

# create our own custom external link classes
class DoxygenUrlGenerator(epydoc.docwriter.xlink.UrlGenerator):

  def get_url(self, name):
    url = name.replace('_', '__')     \
              .replace('::', '_1_1')
    url = "../c_plus_plus/%s%s.html" % (self.what_is_documented, url)
    return url

class DoxygenClassUrlGenerator(DoxygenUrlGenerator):
  what_is_documented = 'class'

class DoxygenStructUrlGenerator(DoxygenUrlGenerator):
  what_is_documented = 'struct'


def run():
  # register our custom external link classes
  epydoc.docwriter.xlink.register_api('doxyclass', DoxygenClassUrlGenerator())
  epydoc.docwriter.xlink.create_api_role('doxyclass', problematic=False)
  epydoc.docwriter.xlink.register_api('doxystruct', DoxygenClassUrlGenerator())
  epydoc.docwriter.xlink.create_api_role('doxystruct', problematic=False)

  # let epydoc handle the command-line arguments
  import sys
  if '--conf=' not in sys.argv[1:]:
    sys.argv.append('--conf=epydoc.conf')
  options = epydoc.cli.parse_arguments()

  # run it
  epydoc.cli.main(options)

if __name__ == "__main__":
  run()
