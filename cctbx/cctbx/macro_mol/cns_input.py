""" Transfer of CNS reflection files to shared arrays.
"""

__version__="$Revision$"[11:-2]

import exceptions
from cctbx_boost.arraytbx import shared
from cctbx_boost import miller

# <xray-reflection-statement> :==
#   nreflection=<integer>
#   anomalous=<logical>
#   declare <xray-declare-statement> end
#   group <xray-group-statement> end
#   index <integer> <integer> <integer>
#   <word>={<real> <real>|<real>|<integer>}
#
#   <xray-declare-statement> :==
#     name=<word>
#     domain=reciprocal|real
#     type=complex|real|integer
#
#   <xray-group-statement> :==
#     type=hl
#     object=<word>

class CNS_input_Error(exceptions.Exception):
  pass

class CNS_input:

  def __init__(self, file):
    self._readline = file.readline
    self._buffer = []
    self._LineNo = 0
    self._LastWord = ""
    self.level = 0

  def getNextWord(self, word_len = 0):
    while (len(self._buffer) == 0):
      line = self._readline()
      if (line == ""): raise EOFError
      self._LineNo = self._LineNo + 1
      # XXX take care of quotes
      i = line.find("!")
      if (i >= 0): line = line[:i]
      while 1:
        i = line.find("{")
        if (i < 0): break
        self.level = self.level + 1
        while 1:
          j = line.find("}")
          if (j >= 0):
            line = line[:i] + line[j + 1:]
            break
          next_line = self._readline()
          if (next_line == ""): raise EOFError
          line = line + next_line
        self.level = self.level - 1
      line = line.replace("=", " ")
      self._buffer = line.upper().split()
      self._buffer.reverse()
    self._LastWord = self._buffer.pop()
    if (word_len): return self._LastWord[:word_len]
    return self._LastWord

  def getLastWord(self, word_len = 0):
    if (word_len): return self._LastWord[:word_len]
    return self._LastWord

  def raiseError(self, message):
    raise CNS_input_Error, \
      "line %d, word \"%s\": " % (self._LineNo,
                                  self._LastWord) + message

def _ampl_phase_to_complex(ampl, phase):
  import math
  phase = phase * math.pi / 180. # convert to radians
  return complex(ampl * math.cos(phase), ampl * math.sin(phase))

class cns_reciprocal_space_object:

  def __init__(self, name, type):
    self.name = name
    self.type = type
    self.H = shared.miller_Index()
    if   (type == "real"):
      self.data = shared.double()
    elif (type == "complex"):
      self.data = shared.complex_double()
    elif (type == "integer"):
      self.data = shared.int()
    else:
      raise RuntimeError, "Internal Error."

  def __repr__(self):
    return "name=%s type=%s len(data)=%d" % (
      self.name, self.type, self.data.size())

  def append(self, H, value):
    self.H.append(H)
    self.data.append(value)


class cns_reflection_file:

  def __init__(self, nreflections, anomalous,
                     reciprocal_space_objects, groups):
    self.nreflections = nreflections
    self.anomalous = anomalous
    self.reciprocal_space_objects = reciprocal_space_objects
    self.groups = groups

  def __repr__(self):
    result =  "nreflections=%d\n" % (self.nreflections,)
    result += "anomalous=%d\n" % (self.anomalous,)
    for rso in self.reciprocal_space_objects.values():
      result += str(rso) + "\n"
    for g in self.groups:
      result += "group: " + str(g) + "\n"
    return result[:-1]

  def optimize(self):
    rsos = self.reciprocal_space_objects.values()
    for i in xrange(len(rsos)-1):
      h_i = rsos[i].H
      for j in xrange(i+1, len(rsos)):
        h_j = rsos[j].H
        if (h_i is h_j): continue # already optimized
        if (h_i.size() != h_j.size()): continue
        if ((h_i == h_j).count(0) == 0):
          rsos[j].H = h_i

  def join_hl_group(self, group_index=None):
    if (group_index == None):
      assert len(self.groups) == 1
      group_index = 0
    selected_group = self.groups[group_index]
    assert len(selected_group) == 4
    miller_indices = 0
    rsos = []
    joined_sets = []
    for name in selected_group:
      rso = self.reciprocal_space_objects[name]
      assert rso.type == "real"
      rsos.append(rso)
      if (miller_indices == 0): miller_indices = rso.H
      js = miller.join_sets(miller_indices, rso.H)
      assert not js.have_singles()
      joined_sets.append(js)
    hl = shared.hendrickson_lattman()
    for ih in xrange(miller_indices.size()):
      coeff = []
      for ic in xrange(4):
        ih0, ih1 = joined_sets[ic].pairs()[ih]
        assert ih0 == ih
        coeff.append(rsos[ic].data[ih1])
      hl.append(coeff)
    return miller_indices, hl

class CNS_xray_reflection_Reader(CNS_input):

  def __init__(self, file):
    CNS_input.__init__(self, file)

  def _read_nreflections(self):
    self.level = self.level + 1
    word = self.getNextWord()
    try:
      nreflections = int(word)
    except ValueError:
      self.raiseError("integer value expected")
    self.level = self.level - 1
    return nreflections

  def _read_anomalous(self):
    self.level = self.level + 1
    word = self.getNextWord(4)
    if   (word == "TRUE"):
      anomalous = 1
    elif (word == "FALS"):
      anomalous = 0
    else:
      self.raiseError("TRUE or FALSe expected")
    self.level = self.level - 1
    return anomalous

  def _read_declare(self, xray_objects):
    name = None
    domain = None
    type = None
    self.level = self.level + 1
    while 1:
      word = self.getNextWord(4)
      if   (word == "NAME"):
        self.level = self.level + 1
        name = self.getNextWord()
        self.level = self.level - 1
      elif (word == "DOMA"):
        self.level = self.level + 1
        word = self.getNextWord(4)
        if   (word == "RECI"):
          domain = "reciprocal"
        elif (word == "REAL"):
          domain = "real"
        else:
          self.raiseError("unrecognized keyword")
        self.level = self.level - 1
      elif (word == "TYPE"):
        self.level = self.level + 1
        word = self.getNextWord(4)
        if   (word == "REAL"):
          type = "real"
        elif (word == "COMP"):
          type = "complex"
        elif (word == "INTE"):
          type = "integer"
        else:
          self.raiseError("unrecognized keyword")
        self.level = self.level - 1
      elif (word == "END"):
        if (xray_objects.has_key(name)):
          self.raiseError("duplicate declaration of NAME=" + name)
        if (domain != "reciprocal"):
          self.raiseError("real space objects are not supported")
        xray_objects[name] = cns_reciprocal_space_object(name, type)
        break
      else:
        self.raiseError("unrecognized keyword")
    self.level = self.level - 1

  def _read_group(self, xray_objects, groups):
    self.level = self.level + 1
    this_group = []
    while 1:
      word = self.getNextWord(4)
      if   (word == "TYPE"):
        self.level = self.level + 1
        word = self.getNextWord(4)
        if (word != "HL"):
          self.raiseError("TYPE=HL expected")
        self.level = self.level - 1
      elif (word == "OBJE"):
        self.level = self.level + 1
        name = self.getNextWord()
        if (not xray_objects.has_key(name)):
          self.raiseError("reciprocal space object " + name
                          + " does not exist")
        if (xray_objects[name].type != "real"):
          self.raiseError(
            "Hendrickson-Lattman coefficients must be of type real")
        if (xray_objects[name].H.size()):
          self.raiseError(
            "GROUp statement must appear before reflection data")
        this_group.append(name)
        self.level = self.level - 1
      elif (word == "END"):
        if (len(this_group) > 0):
          if (len(this_group) != 4):
            self.raiseError(
              "there must be exactly four Hendrickson-Lattman coefficients")
          groups.append(this_group)
        break
      else:
        self.raiseError("unrecognized keyword")
    self.level = self.level - 1

  def load(self):
    gNW = self.getNextWord
    gLW = self.getLastWord

    nreflections = None
    anomalous = None
    reciprocal_space_objects = {}
    groups = []

    current_hkl = None
    reuse_word = 0

    try:
      while 1:
        if (not reuse_word): word = gNW(4)
        reuse_word = 0
        if (word == "INDE"):
          self.level = self.level + 1
          H = [None] * 3
          for i in xrange(3):
            word = gNW()
            try:
              H[i] = int(word)
            except ValueError:
              self.raiseError("integer values expected for hkl")
          current_hkl = tuple(H)
          self.level = self.level - 1
        elif (word == "NREF"):
          nreflections = self._read_nreflections()
        elif (word == "ANOM"):
          anomalous = self._read_anomalous()
        elif (word == "DECL"):
          self._read_declare(reciprocal_space_objects)
        elif (word == "GROU"):
          self._read_group(reciprocal_space_objects, groups)
        else:
          word = gLW()
          try:
            type = reciprocal_space_objects[word].type
          except KeyError:
            self.raiseError("unrecognized keyword")
          self.level = self.level + 1
          name = word
          n = 1
          if (type == "complex"): n = 2
          for i in xrange(n):
            word = gNW()
            if (type == "integer"):
              try:
                value = int(word)
              except ValueError:
                self.raiseError("integer value expected for array " + name)
            else:
              try:
                if (i == 0): value = float(word)
                else:        value = _ampl_phase_to_complex(value,float(word))
              except ValueError:
                if (i == 0):
                  self.raiseError("floating-point value expected for array "
                                  + name)
                # declared complex but only real part given
                reuse_word = 1
          reciprocal_space_objects[name].append(current_hkl, value)
          self.level = self.level - 1

    except EOFError:
      if (self.level != 0): raise CNS_input_Error, "premature end-of-file"

    return cns_reflection_file(nreflections, anomalous,
                               reciprocal_space_objects, groups)


def summary(cns_reflection_file):
  reader = CNS_xray_reflection_Reader(cns_reflection_file)
  reflection_file = reader.load()
  print reflection_file

def run():
  import sys, os
  if (len(sys.argv) == 1):
    summary(sys.stdin)
  else:
    for file in sys.argv[1:]:
      print file + ":"
      f = open(file, "r")
      summary(f)
      f.close()
      print
  t = os.times()
  print "u+s,u,s:", t[0] + t[1], t[0], t[1]

if (__name__ == "__main__"):
  run()
