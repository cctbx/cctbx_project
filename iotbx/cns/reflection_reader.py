"Transfer of CNS reflection files to flex arrays."

from iotbx.cns.crystal_symmetry_utils import \
  re_sg_uc, re_uc_sg, crystal_symmetry_from_re_match
from cctbx import crystal
from cctbx import miller
from cctbx.array_family import flex
from libtbx import complex_math
from libtbx import easy_pickle
import re
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

class CNS_input_Error(Exception):
  pass

class CNS_input(object):

  def __init__(self, file):
    self._readline = file.readline
    self._buffer = []
    self._LineNo = 0
    self._LastWord = ""
    self.level = 0
    self.comments = []
    self.remarks = []
    self.cryst1s = []

  def getNextWord(self, word_len = 0):
    while (len(self._buffer) == 0):
      line = self._readline()
      if (line == ""): raise EOFError
      self._LineNo += 1
      if (line.lstrip().upper().startswith("REMARK ")):
        self.remarks.append(line)
        continue
      if (line.startswith("CRYST1")
            and len(line) >= 57
              and " +-.0123456789".find(line[6]) >= 0):
        self.cryst1s.append(line)
        continue
      # XXX take care of quotes
      i = line.find("!")
      if (i >= 0): line = line[:i]
      while 1:
        i = line.find("{")
        if (i < 0): break
        self.level = self.level + 1
        while 1:
          j = line.find("}", i+1)
          if (j >= 0):
            self.comments.append(line[i:j+1])
            line = line[:i] + line[j + 1:]
            break
          next_line = self._readline()
          if (next_line == ""): raise EOFError
          line += next_line
        self.level -= 1
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
    raise CNS_input_Error(
      "line %d, word \"%s\": " % (self._LineNo, self._LastWord) + message)

  def raiseError_floating_point(self, name):
    self.raiseError("floating-point value expected for array " + name)

class cns_reciprocal_space_object(object):

  def __init__(self, name, type):
    self.name = name
    self.type = type
    self.indices = flex.miller_index()
    self.has_non_zero_phases = None
    if   (type == "real"):
      self.data = flex.double()
    elif (type == "complex"):
      self.data = flex.complex_double()
      self.has_non_zero_phases = False
    elif (type == "integer"):
      self.data = flex.int()
    else:
      raise RuntimeError("Internal Error.")

  def show_summary(self, f=None, prefix=""):
    if (f is None): f = sys.stdout
    print >> f, prefix + "name=%s type=%s len(data)=%d" % (
      self.name, self.type, self.data.size())

  def append(self, h, value):
    self.indices.append(h)
    self.data.append(value)

  def is_real(self, use_name_as_hint=True):
    if (self.type == "real"): return True
    if (self.type != "complex"): return False
    if (not self.has_non_zero_phases): return True
    if (self.name.lower() in [
      "fobs", "f_obs", "iobs", "i_obs", "obs"]): return True
    return False

  def real_data(self, use_name_as_hint=True):
    assert self.is_real(use_name_as_hint=use_name_as_hint)
    if (self.type == "real"): return self.data
    return flex.abs(self.data)

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
      anomalous = True
    elif (word == "FALS"):
      anomalous = False
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
    reuse_word = False
    n_words_processed = 0
    try:
      while True:
        if (not reuse_word):
          word = gNW(4)
          n_words_processed += 1
        reuse_word = False
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
          rso = to_obj.reciprocal_space_objects.get(word)
          if (rso is None):
            self.raiseError("unrecognized keyword")
          type = rso.type
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
              if (i == 0):
                try:
                  value = float(word)
                except ValueError:
                  self.raiseError_floating_point(name)
              else:
                try:
                  phase = float(word)
                except ValueError:
                  reuse_word = True # declared complex but only real part given
                else:
                  value = complex_math.polar((value, phase), deg=True)
                  if (phase != 0):
                    rso.has_non_zero_phases = True
          rso.append(current_hkl, value)
          self.level = self.level - 1
    except EOFError:
      if (self.level != 0): raise CNS_input_Error("premature end-of-file")
    if (n_words_processed == 0):
      raise CNS_input_Error("empty file")

class cns_reflection_file(object):

  def __init__(self, file_handle):
    reader = CNS_xray_reflection_Reader(file_handle)
    reader.load(self)
    self.optimize()
    self.remarks = reader.remarks
    self.cryst1s = reader.cryst1s
    self.comments = reader.comments

  def space_group_from_remark_symop(self):
    from cctbx import sgtbx
    result = None
    for remark in self.remarks:
      remark = remark.lstrip()[6:].strip().replace(" ", "").lower()
      if (    remark.startswith("symop(")
          and remark.endswith(")")):
        s = remark[6:-1]
        try:
          s = sgtbx.rt_mx(s)
        except RuntimeError:
          pass
        else:
          if (result is None):
            result = sgtbx.space_group()
          result.expand_smx(s)
    return result

  def crystal_symmetry_from_remark_uc_sg(self):
    sg = self.space_group_from_remark_symop()
    for remark in self.remarks:
      remark = remark.lstrip()[6:].strip()
      m = re.match(re_uc_sg, remark)
      if (m is None): continue
      result = crystal_symmetry_from_re_match(m=m, i_uc=1, i_sg=7)
      if (result is not None):
        if (sg is not None):
          result = crystal.symmetry(
            unit_cell=result.unit_cell(),
            space_group=sg)
        return result
    return None

  def crystal_symmetry_from_cryst1(self):
    for record in self.cryst1s:
      from iotbx.pdb import cryst1_interpretation
      result = cryst1_interpretation.crystal_symmetry(cryst1_record=record)
      if (result is not None
            and result.unit_cell() is not None
            and result.space_group_info() is not None):
        return result
    return None

  def crystal_symmetry_from_comments(self):
    for comment in self.comments:
      m = re.match(r'\{\s+' + re_sg_uc, comment)
      if (m is None): continue
      result = crystal_symmetry_from_re_match(m=m)
      if (result is not None): return result
    return None

  def crystal_symmetry(self,
        crystal_symmetry=None,
        force_symmetry=False):
    self_symmetry = self.crystal_symmetry_from_remark_uc_sg()
    if (self_symmetry is None):
      self_symmetry = self.crystal_symmetry_from_cryst1()
    if (self_symmetry is None):
      self_symmetry = self.crystal_symmetry_from_comments()
    if (crystal_symmetry is None):
      return self_symmetry
    if (self_symmetry is None):
      return crystal_symmetry
    return self_symmetry.join_symmetry(
      other_symmetry=crystal_symmetry,
      force=force_symmetry)

  def show_summary(self, f=None, prefix=""):
    if (f is None): f = sys.stdout
    print >> f, prefix + "nreflections=%d" % self.nreflections
    print >> f, prefix + "anomalous=" + str(self.anomalous)
    for rso in self.reciprocal_space_objects.values():
      rso.show_summary(f=f, prefix=prefix)
    for g in self.groups:
      print >> f, prefix + "group: " + str(g)

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
    names = []
    miller_indices = 0
    rsos = []
    matches = []
    for name in selected_group:
      names.append(name)
      rso = self.reciprocal_space_objects[name]
      assert rso.type == "real"
      rsos.append(rso)
      if (type(miller_indices) == type(0)): miller_indices = rso.indices
      match = miller.match_indices(miller_indices, rso.indices)
      assert not match.have_singles()
      matches.append(match)
    hl = flex.hendrickson_lattman()
    for ih in xrange(miller_indices.size()):
      coeff = []
      for ic in xrange(4):
        ih0, ih1 = matches[ic].pairs()[ih]
        assert ih0 == ih
        coeff.append(rsos[ic].data[ih1])
      hl.append(coeff)
    return names, miller_indices, hl

  def _as_miller_array(self, crystal_symmetry, miller_indices,
                             data, sigmas=None, obs_type=None):
    result = miller.set(
      crystal_symmetry=crystal_symmetry,
      indices=miller_indices,
      anomalous_flag=self.anomalous)
    if (self.anomalous is None):
      result = result.auto_anomalous()
    result = result.array(data=data, sigmas=sigmas)
    if (obs_type is not None):
      assert obs_type in ("f", "i")
      if (obs_type == "f"):
        result.set_observation_type_xray_amplitude()
      else:
        result.set_observation_type_xray_intensity()
    return result

  def as_miller_arrays(self,
        crystal_symmetry=None,
        force_symmetry=False,
        merge_equivalents=True,
        base_array_info=None):
    crystal_symmetry = self.crystal_symmetry(
      crystal_symmetry=crystal_symmetry,
      force_symmetry=force_symmetry)
    if (crystal_symmetry is None):
      crystal_symmetry = crystal.symmetry(
        unit_cell=None,
        space_group_info=None)
    if (base_array_info is None):
      base_array_info = miller.array_info(source_type="cns_reflection_file")
    result = []
    done = set()
    for group_index in xrange(len(self.groups)):
      names, miller_indices, hl = self.join_hl_group(group_index)
      result.append(self._as_miller_array(
        crystal_symmetry, miller_indices, hl).set_info(
          base_array_info.customized_copy(labels=names)))
      for name in names:
        done.add(name)
    real_arrays = {}
    for rso in self.reciprocal_space_objects.values():
      if (rso.name in done): continue
      if (not rso.is_real()): continue
      real_arrays[rso.name.lower()] = rso
    for obs,sigma,obs_type in group_obs_sigma(real_arrays):
      result.append(self._as_miller_array(
        crystal_symmetry, obs.indices,
        obs.real_data(), sigma.real_data(), obs_type).set_info(
          base_array_info.customized_copy(labels=[obs.name, sigma.name])))
      done.add(obs.name)
      done.add(sigma.name)
    for rso in self.reciprocal_space_objects.values():
      if (rso.name in done): continue
      result.append(self._as_miller_array(
        crystal_symmetry, rso.indices, rso.data).set_info(
          base_array_info.customized_copy(labels=[rso.name])))
      done.add(rso.name)
    return result

def group_obs_sigma(real_arrays):
  result = []
  done = set()
  i = real_arrays.get("iobs")
  if (i is not None):
    for name in ["sigi", "sigmai", "sigiobs"]:
      s = real_arrays.get(name)
      if (s is not None and i.indices.all_eq(s.indices)):
        result.append((i,s,"i"))
        done.add(i.name)
        done.add(s.name)
        break
  for name,f in real_arrays.items():
    if (f.name in done): continue
    rest = None
    for prefix in ("fobs", "f_obs", "f", ""):
      if (name.startswith(prefix)):
        rest = name[len(prefix):]
        break
    if (rest is None):
      continue
    for prefix in ("sigma", "sig", "s"):
      s = real_arrays.get(prefix+rest)
      if (s is not None):
        break
    if (s is None):
      for prefix in ("sigma", "sig"):
        s = real_arrays.get(prefix+name)
        if (s is not None):
          break
    if (s is not None and f.indices.all_eq(s.indices)):
      result.append((f,s,"f"))
      done.add(f.name)
      done.add(s.name)
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
