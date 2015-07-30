from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 06/18/2015
Description : Progress bar. Adapted from DIALS progress bar (nearly a copy 
              thereof with small alterations for IOTA purposes).
'''

import sys

class ProgressBarTimer:
  """ A simple timer for the progress bar. """

  def __init__(self):
    """ Init the progress bar timer. """
    from time import time
    self._last_time = time()
    self._last_perc = 0
    self._update_period = 0.5
    self._n_seconds_left = -1

  def update(self, percent):
    """ Update the timer. """
    from time import time

    # Get the current time diff between last time
    curr_time = time()
    diff_time = curr_time - self._last_time

    # Only update after certain period or at 100%
    if percent < 0: percent = 0
    if percent > 100: percent = 100
    if diff_time >= self._update_period or percent >= 100:

      # Check the difference in percentage and calculate
      # number of seconds remaining
      diff_perc = percent - self._last_perc
      if (diff_perc == 0):
        self._n_seconds_left = 0
      else:
        self._n_seconds_left = diff_time * (100 - percent) / diff_perc

    # Return number of seconds
    return self._n_seconds_left


class ProgressBar:
  """ A command line progress bar. """

  def __init__(self, title=None, ticker=True, bar=True, estimate_time=False,
               indent=0, length=80):
    """ Init the progress bar parameters. """

    # Set the parameters
    self._title = title
    self._indent = indent
    self._ticker = ticker
    self._estimate_time = estimate_time
    self._bar = bar
    self._length = length

    self._timer = ProgressBarTimer()
    self._start_time = self._timer._last_time

    # Print 0 percent
    self.update(0, 0)

  def update(self, fpercent, spin_update):
    """ Update the progress bar with a percentage. """
    import sys
    from math import ceil

    # do not update if not a tty
    if not sys.stdout.isatty():
      return

    # Get integer percentage
    percent = int(fpercent)
    if percent < 0: percent = 0
    if percent > 100: percent = 100

    # Add a percentage counter
    right_str = ''
    left_str = ''
    if sys.stdout.isatty():
      left_str = '\r'
    left_str += ' ' * self._indent

    # Add a title if given
    if self._title:
      left_str += self._title + ': '

    left_str += '{0: >3}%'.format(percent)

    # Add a ticker
    if self._ticker:
      tick = ['-   ', '--  ', '--- ', '----', '--- ', '--  ']
      left_str += ' '
      left_str += '[ {0} ]'.format(tick[spin_update % 6])

    # Add a timer
    if self._estimate_time:
      n_seconds_left = self._timer.update(fpercent)
      if n_seconds_left < 0:
        n_seconds_left = '?'
      else:
        n_seconds_left = int(ceil(n_seconds_left))
      right_str = ' ' + 'est: {0}s'.format(n_seconds_left) + right_str

    # Add a bar
    if self._bar:
      bar_length = self._length - (len(left_str) + len(right_str)) - 5
      n_char = int(percent * bar_length / 100)
      n_space = bar_length - n_char
      left_str += ' '
      left_str += '[ {0}>{1} ]'.format('=' * n_char, ' ' * n_space)

    # Append strings
    progress_str = left_str + right_str

    # Print progress string to stdout
    sys.stdout.write(progress_str)
    sys.stdout.flush()

  def finished(self, string=None):
    """ The progress bar is finished. """
    if string:
      self._title = string
    else:
      string = ''

    ''' Print the 'end of comand' string.'''
    from sys import stdout
    from time import time


    if self._estimate_time:
      # Get the time string
      time_string = '{0:.2f}s'.format(time() - self._start_time)

      # Truncate the string
      max_length = self._length - self._indent - len(time_string) - 1
      string = string[:max_length]

      # Add an indent and a load of dots and then the time string
      dot_length = 1 + max_length - len(string)
      string = (' ' * self._indent) + string
      string = string + '.' * (dot_length)
      string = string + time_string

    else:

      # Truncaet the string
      max_length = self._length - self._indent
      string = string[:max_length]

      # Add a load of dots
      dot_length = max_length - len(string)
      string = (' ' * self._indent) + string
      string = string + '.' * (dot_length)

    # Write the string to stdout
    if stdout.isatty():
      string = '\r' + string + '\n'
    else:
      string = string + '\n'
    stdout.write(string)
    stdout.flush()

class Command(object):
  '''Class to nicely print out a command with timing info.'''

  # Variables available in class methods
  indent = 0
  max_length = 80
  print_time = True

  @classmethod
  def start(self, string):
    ''' Print the 'start command' string.'''
    from sys import stdout
    from time import time
    # from termcolor import colored

    # Get the command start time
    self._start_time = time()

    # do not output if not a tty
    if not stdout.isatty():
      return

    # Truncate the string to the maximum length
    max_length = self.max_length - self.indent - 3
    string = string[:max_length]
    string = (' ' * self.indent) + string + '...'

    # Write the string to stdout
    stdout.write(string)
    stdout.flush()

  @classmethod
  def end(self, string):
    ''' Print the 'end of comand' string.'''
    from sys import stdout
    from time import time
    #from termcolor import colored

    # Check if we want to print the time or not
    if self.print_time:

      # Get the time string
      time_string = '{0:.2f}s'.format(time() - self._start_time)

      # Truncate the string
      max_length = self.max_length - self.indent - len(time_string) - 1
      string = string[:max_length]

      # Add an indent and a load of dots and then the time string
      dot_length = 1 + max_length - len(string)
      string = (' ' * self.indent) + string
      string = string + '.' * (dot_length)
      string = string + time_string

    else:

      # Truncaet the string
      max_length = self.max_length - self.indent
      string = string[:max_length]

      # Add a load of dots
      dot_length = max_length - len(string)
      string = (' ' * self.indent) + string
      string = string + '.' * (dot_length)

    # Write the string to stdout
    if stdout.isatty():
      string = '\r' + string + '\n'
    else:
      string = string + '\n'
    stdout.write(string)
    stdout.flush()

try:
  import termcolor
except ImportError:
  termcolor = None

def coloured(text, *args, **kwargs):
  import sys
  if not sys.stdout.isatty() or termcolor is None:
    return text
  return termcolor.colored(text, *args, **kwargs)

def heading(text):
  return coloured(text, attrs=['bold'])

if __name__ == '__main__':
  import time

  p = ProgressBar()
  h = 200
  k = 100 / h

  for j in range(h):
    i = j * k
    p.update(i, j)
    time.sleep(0.25)

  p.finished()

  Command.start("Starting to do a command")
  time.sleep(1)
  Command.end("Ending the command")
