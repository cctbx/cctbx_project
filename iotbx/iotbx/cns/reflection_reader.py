"Transfer of CNS reflection files to flex arrays."

import exceptions
from cctbx import crystal
from cctbx import miller
from cctbx.array_family import flex
from scitbx.python_utils import complex_math
from scitbx.python_utils import easy_pickle
import sys

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

class cns_reciprocal_space_object:

  def __init__(self, name, type):
    self.name = name
    self.type = type
    self.indices = flex.miller_index()
    if   (type == "real"):
      self.data = flex.double()
    elif (type == "complex"):
      self.data = flex.complex_double()
    elif (type == "integer"):
      self.data = flex.int()
    else:
      raise RuntimeError, "Internal Error."

  def show_summary(self, f=sys.stdout):
    print >> f, "name=%s type=%s len(data)=%d" % (
      self.name, self.type, self.data.size())

  def append(self, h, value):
    self.indices.append(h)
    self.data.append(value)

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
        if (xray_objects[name].indices.size()):
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

  def load(self, to_obj):
    gNW = self.getNextWord
    gLW = self.getLastWord
    to_obj.nreflections = None
    to_obj.anomalous = None
    to_obj.reciprocal_space_objects = {}
    to_obj.groups = []
    current_hkl = None
    reuse_word = 0
    try:
      while 1:
        if (not reuse_word): word = gNW(4)
        reuse_word = 0
        if (word == "INDE"):
          self.level = self.level + 1
          h = [None] * 3
          for i in xrange(3):
            word = gNW()
            try:
              h[i] = int(word)
            except ValueError:
              self.raiseError("integer values expected for hkl")
          current_hkl = tuple(h)
          self.level = self.level - 1
        elif (word == "NREF"):
          to_obj.nreflections = self._read_nreflections()
        elif (word == "ANOM"):
          to_obj.anomalous = self._read_anomalous()
        elif (word == "DECL"):
          self._read_declare(to_obj.reciprocal_space_objects)
        elif (word == "GROU"):
          self._read_group(to_obj.reciprocal_space_objects, to_obj.groups)
        else:
          word = gLW()
          try:
            type = to_obj.reciprocal_space_objects[word].type
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
                if (i == 0):
                  value = float(word)
                else:
                  value = complex_math.polar((value, float(word)), deg=0001)
              except ValueError:
                if (i == 0):
                  self.raiseError("floating-point value expected for array "
                                  + name)
                # declared complex but only real part given
                reuse_word = 1
          to_obj.reciprocal_space_objects[name].append(current_hkl, value)
          self.level = self.level - 1
    except EOFError:
      if (self.level != 0): raise CNS_input_Error, "premature end-of-file"

def as_miller_array(crystal_symmetry, anomalous_flag,
                    miller_indices, data, info):
  return (miller.array(
    miller_set=miller.set(
      crystal_symmetry=crystal_symmetry,
      indices=miller_indices,
      anomalous_flag=anomalous_flag),
    data=data)
    .set_info(info))

class cns_reflection_file:

  def __init__(self, file_handle):
    reader = CNS_xray_reflection_Reader(file_handle)
    reader.load(self)
    self.optimize()

  def show_summary(self, f=sys.stdout):
    print >> f, "nreflections=%d" % (self.nreflections,)
    print >> f, "anomalous=%d" % (self.anomalous,)
    for rso in self.reciprocal_space_objects.values():
      rso.show_summary(f)
    for g in self.groups:
      print >> f, "group: " + str(g)

  def optimize(self):
    rsos = self.reciprocal_space_objects.values()
    for i in xrange(len(rsos)-1):
      h_i = rsos[i].indices
      for j in xrange(i+1, len(rsos)):
        h_j = rsos[j].indices
        if (flex.order(h_i, h_j) == 0):
          rsos[j].indices = h_i

  def join_hl_group(self, group_index=None):
    if (group_index is None):
      assert len(self.groups) == 1
      group_index = 0
    selected_group = self.groups[group_index]
    assert len(selected_group) == 4
    miller_indices = 0
    rsos = []
    matches = []
    for name in selected_group:
      rso = self.reciprocal_space_objects[name]
      assert rso.type == "real"
      rsos.append(rso)
      if (type(miller_indices) == type(0)): miller_indices = rso.indices
      match = miller.match_indices(miller_indices, rso.indices)
      assert not match.have_singles()
      matches.append(match)
    hl = flex.hendrickson_lattman()
    for ih in miller_indices.indices():
      coeff = []
      for ic in xrange(4):
        ih0, ih1 = matches[ic].pairs()[ih]
        assert ih0 == ih
        coeff.append(rsos[ic].data[ih1])
      hl.append(coeff)
    return miller_indices, hl

  def as_miller_arrays(self, crystal_symmetry, force_symmetry=00000,
                             info_prefix=""):
    result = []
    done = {}
    for group_index in xrange(len(self.groups)):
      miller_indices, hl = self.join_hl_group(group_index)
      info = info_prefix + "hl_group_%d" % (group_index+1,)
      result.append(as_miller_array(
        crystal_symmetry, self.anomalous, miller_indices, hl, info))
      for name in self.groups[group_index]:
        done[name] = 1
    for rso in self.reciprocal_space_objects.values():
      if rso.name in done: continue
      info = info_prefix + rso.name.lower()
      result.append(as_miller_array(
        crystal_symmetry, self.anomalous, rso.indices, rso.data, info))
      done[rso.name] = 1
    return result

def run(args):
  import os
  to_pickle = "--pickle" in args
  for file_name in args:
    if (file_name.startswith("--")): continue
    print file_name + ":"
    f = open(file_name, "r")
    t0 = os.times()
    reflection_file = cns_reflection_file(f)
    tn = os.times()
    t_parse = tn[0]+tn[1]-t0[0]-t0[1]
    f.close()
    reflection_file.show_summary()
    print
    crystal_symmetry = crystal.symmetry((), "P 1")
    miller_arrays = reflection_file.as_miller_arrays(crystal_symmetry)
    for miller_array in miller_arrays:
      miller_array.show_summary()
      print
    if (to_pickle):
      pickle_file_name = os.path.split(file_name)[1] + ".pickle"
      t0 = os.times()
      easy_pickle.dump(pickle_file_name, reflection_file)
      tn = os.times()
      t_dump = tn[0]+tn[1]-t0[0]-t0[1]
      t0 = os.times()
      easy_pickle.load(pickle_file_name)
      tn = os.times()
      t_load = tn[0]+tn[1]-t0[0]-t0[1]
      print "parse: %.2f, dump: %.2f, load: %.2f" % (t_parse, t_dump, t_load)
    print
  t = os.times()
  print "u+s,u,s: %.2f %.2f %.2f" % (t[0] + t[1], t[0], t[1])
