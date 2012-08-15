"""optik

A powerful, extensible, and easy-to-use command-line parser for Python.

By Greg Ward <gward@python.net>

See http://optik.sourceforge.net/
"""
from __future__ import division

# Copyright (c) 2001-2006 Gregory P. Ward.  All rights reserved.
# See the README.txt distributed with Optik for licensing terms.

__version__ = "1.5.1"


# Re-import these for convenience
from optik.option import Option
from optik.option_parser import *
from optik.help import *
from optik.errors import *

from optik import option, option_parser, help, errors
__all__ = (option.__all__ +
           option_parser.__all__ +
           help.__all__ +
           errors.__all__)


# Some day, there might be many Option classes.  As of Optik 1.3, the
# preferred way to instantiate Options is indirectly, via make_option(),
# which will become a factory function when there are many Option
# classes.
make_option = Option
