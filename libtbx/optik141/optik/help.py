import textwrap

__revision__ = "$Id$"

__all__ = ['HelpFormatter', 'IndentedHelpFormatter', 'TitledHelpFormatter']


class HelpFormatter:

    """
    Abstract base class for formatting option help.  OptionParser
    instances should use one of the HelpFormatter subclasses for
    formatting help; by default IndentedHelpFormatter is used.

    Instance attributes:
      indent_increment : int
        the number of columns to indent per nesting level
      max_help_position : int
        the maximum starting column for option help text
      help_position : int
        the calculated starting column for option help text;
        initially the same as the maximum
      width : int
        total number of columns for output
      level : int
        current indentation level
      current_indent : int
        current indentation level (in columns)
      help_width : int
        number of columns available for option help text (calculated)
    """

    def __init__ (self,
                  indent_increment,
                  max_help_position,
                  width,
                  short_first):
        self.indent_increment = indent_increment
        self.help_position = self.max_help_position = max_help_position
        self.width = width
        self.current_indent = 0
        self.level = 0
        self.help_width = width - max_help_position
        self.short_first = short_first

    def indent (self):
        self.current_indent += self.indent_increment
        self.level += 1

    def dedent (self):
        self.current_indent -= self.indent_increment
        assert self.current_indent >= 0, "Indent decreased below 0."
        self.level -= 1

    def format_usage (self, usage):
        raise NotImplementedError, "subclasses must implement"

    def format_heading (self, heading):
        raise NotImplementedError, "subclasses must implement"

    def format_description (self, description):
        desc_width = self.width - self.current_indent
        indent = " "*self.current_indent
        return textwrap.fill(description, desc_width,
                             initial_indent=indent,
                             subsequent_indent=indent)

    def format_option (self, option):
        # The help for each option consists of two parts:
        #   * the opt strings and metavars
        #     eg. ("-x", or "-fFILENAME, --file=FILENAME")
        #   * the user-supplied help string
        #     eg. ("turn on expert mode", "read data from FILENAME")
        #
        # If possible, we write both of these on the same line:
        #   -x      turn on expert mode
        #
        # But if the opt string list is too long, we put the help
        # string on a second line, indented to the same column it would
        # start in if it fit on the first line.
        #   -fFILENAME, --file=FILENAME
        #           read data from FILENAME
        result = []
        opts = option.option_strings
        opt_width = self.help_position - self.current_indent - 2
        if len(opts) > opt_width:
            opts = "%*s%s\n" % (self.current_indent, "", opts)
            indent_first = self.help_position
        else:                       # start help on same line as opts
            opts = "%*s%-*s  " % (self.current_indent, "", opt_width, opts)
            indent_first = 0
        result.append(opts)
        if option.help:
            help_lines = textwrap.wrap(option.help, self.help_width)
            result.append("%*s%s\n" % (indent_first, "", help_lines[0]))
            result.extend(["%*s%s\n" % (self.help_position, "", line)
                           for line in help_lines[1:]])
        elif opts[-1] != "\n":
            result.append("\n")
        return "".join(result)

    def store_option_strings (self, parser):
        self.indent()
        max_len = 0
        for opt in parser.option_list:
            strings = self.format_option_strings(opt)
            opt.option_strings = strings
            max_len = max(max_len, len(strings) + self.current_indent)
        self.indent()
        for group in parser.option_groups:
            for opt in group.option_list:
                strings = self.format_option_strings(opt)
                opt.option_strings = strings
                max_len = max(max_len, len(strings) + self.current_indent)
        self.dedent()
        self.dedent()
        self.help_position = min(max_len + 2, self.max_help_position)

    def format_option_strings (self, option):
        """Return a comma-separated list of option strings & metavariables."""
        if option.takes_value():
            metavar = option.metavar or option.dest.upper()
            short_opts = [sopt + metavar for sopt in option._short_opts]
            long_opts = [lopt + "=" + metavar for lopt in option._long_opts]
        else:
            short_opts = option._short_opts
            long_opts = option._long_opts

        if self.short_first:
            opts = short_opts + long_opts
        else:
            opts = long_opts + short_opts

        return ", ".join(opts)

class IndentedHelpFormatter (HelpFormatter):
    """Format help with indented section bodies.
    """

    def __init__ (self,
                  indent_increment=2,
                  max_help_position=24,
                  width=80,
                  short_first=1):
        HelpFormatter.__init__(
            self, indent_increment, max_help_position, width, short_first)

    def format_usage (self, usage):
        return "usage: %s\n" % usage

    def format_heading (self, heading):
        return "%*s%s:\n" % (self.current_indent, "", heading)


class TitledHelpFormatter (HelpFormatter):
    """Format help with underlined section headers.
    """

    def __init__ (self,
                  indent_increment=0,
                  max_help_position=24,
                  width=80,
                  short_first=0):
        HelpFormatter.__init__ (
            self, indent_increment, max_help_position, width, short_first)

    def format_usage (self, usage):
        return "%s  %s\n" % (self.format_heading("Usage"), usage)

    def format_heading (self, heading):
        return "%s\n%s\n" % (heading, "=-"[self.level] * len(heading))
