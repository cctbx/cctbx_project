from __future__ import absolute_import, division, print_function

try:
  import argparse

except ImportError:
  from libtbx.utils import Sorry
  raise Sorry("The argparse module is not available. Check Python version >= 2.7")

import libtbx.phil


def phil_from_string(string):

  try:
    return libtbx.phil.parse( string )

  except RuntimeError:
    raise ValueError("Incorrect PHIL format")


def phil_from_file(file_name):

  with open( file_name ) as f:
    string = f.read()

  return phil_from_string( string = string )


class PhilArgumentFactory(object):

  def __init__(self, master_phil):

    from libtbx.phil import command_line
    self.interpreter = command_line.argument_interpreter( master_phil = master_phil )
    self.file_handlers = []


  def add_file_handler(self, handler):

    self.file_handlers.append( handler )


  def __call__(self, argument):

    if "=" in argument:
      try:
        return self.interpreter.process( arg = argument )

      except Exception as e:
        raise argparse.ArgumentTypeError( "Cannot interpret assignment '%s': %s" % ( argument, e ) )

    import os.path

    if os.path.isfile( argument ):
      for handler in self.file_handlers:
        try:
          return handler( argument = argument )

        except ValueError:
          continue

      try:
        return phil_from_file( file_name = argument )

      except ValueError:
        pass

      raise argparse.ArgumentTypeError( "Does not recognize file: %s" % argument )

    raise argparse.ArgumentTypeError( "Cannot interpret argument: %s" % argument )


# File recognition
class ExtensionFileHandler(object):

  def __init__(self, extensions, philpath):

    self.extensions = set( extensions )
    self.philpath = philpath


  def __call__(self, argument):

    import os.path

    if os.path.splitext( argument )[1] in self.extensions:
      return phil_from_string( string = "%s=%s" % ( self.philpath, argument ) )

    raise ValueError("Unrecognized file: %s" % argument)


# Parser actions
class PhilPrintAction(argparse.Action):
  """
  Action to print phil
  """

  def __init__(
    self,
    option_strings,
    dest = argparse.SUPPRESS,
    default = argparse.SUPPRESS,
    help = None,
    ):

    super( PhilPrintAction, self ).__init__(
      option_strings = option_strings,
      dest = dest,
      default = default,
      nargs = 0,
      help = help,
      )


  def __call__(self, parser, namespace, values, option_string = None):

    parser.print_phil()
    parser.exit()


class CaptureStdinPhilAction(argparse.Action):
  """
  Action to capture stdin
  """

  def __init__(self, option_strings, dest, default = argparse.SUPPRESS, help = None):

    super( CaptureStdinPhilAction, self ).__init__(
      option_strings = option_strings,
      dest = dest,
      default = default,
      nargs = 0,
      help = help,
      )


  def __call__(self, parser, namespace, values, option_string = None):

    import sys

    string = sys.stdin.read()

    try:
      value = libtbx.phil.parse( string )

    except RuntimeError:
      parser.error( "Could not interpret input from stdin: %s\n" % string )

    setattr( namespace, self.dest, value )


class Parser(argparse.ArgumentParser):
  """
  Extended parser that handles help and phil
  """

  def __init__(self, master_phil, *args, **kwargs):

    super( Parser, self ).__init__( *args, **kwargs )
    self.master_phil = master_phil
    self.factory = PhilArgumentFactory( master_phil = self.master_phil )

    self.add_argument(
      "phils",
      nargs = "+",
      metavar = "PHIL",
      action = "store",
      default = [],
      type = self.factory,
      help = "PHIL argument (file name or PHIL command line assignment)"
      )

    prefix = self.prefix_chars[0]
    self.add_argument(
      prefix * 2 + "show-defaults",
      action = PhilPrintAction,
      help= "print PHIL and exit",
      )
    self.add_argument(
      prefix + "i", prefix * 2 + "stdin",
      action = CaptureStdinPhilAction,
      help= "read PHIL from stdin as well",
      )


  def recognize_file_type_as(self, extensions, philpath):

    self.factory.add_file_handler(
      ExtensionFileHandler( extensions = extensions, philpath = philpath ),
      )


  def print_phil(self, file = None):

    super( Parser, self )._print_message(
      self.master_phil.as_str() + "\n",
      file,
      )


  def parse_args(self, args = None, namespace = None):

    args = super( Parser, self ).parse_args(
      args = args,
      namespace = namespace
      )

    if hasattr( args, "stdin" ):
      args.phils.append( args.stdin )
      delattr( args, "stdin" )

    try:
      merged = self.master_phil.fetch( sources = args.phils )

    except Exception as e:
      self.error( "Error while merging arguments: %s" % e )

    delattr( args, "phils" )
    args.phil = merged

    try:
      args.phil_extract = args.phil.extract()

    except RuntimeError as e:
      self.error( str( e ) )

    return args

