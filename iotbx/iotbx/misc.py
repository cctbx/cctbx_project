class detect_binary_file:

  def __init__(self, monitor_initial=1000, max_fraction_non_ascii=0.05):
    self.monitor_initial = monitor_initial
    self.max_fraction_non_ascii = max_fraction_non_ascii
    self.n_ascii_characters = 0
    self.n_non_ascii_characters = 0
    self.status = None

  def is_binary_file(self, line):
    if (self.monitor_initial > 0):
      for c in line:
        if (1 < ord(c) < 128):
          self.n_ascii_characters += 1
        else:
          self.n_non_ascii_characters += 1
        self.monitor_initial -= 1
        if (self.monitor_initial == 0):
          if (  self.n_non_ascii_characters
              > self.n_ascii_characters * self.max_fraction_non_ascii):
            self.status = 0001
          else:
            self.status = 00000
          break
    return self.status
