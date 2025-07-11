"""Interpret atom names"""
from __future__ import absolute_import, division, print_function
import sys

class dict_with_add(dict):

  def __add__(self, other):
    result = dict_with_add(self)
    result.update(other)
    return result

def alternative_hydrogen_pattern(pattern):
  if (len(pattern) > 1
      and pattern[1] == "h"
      and pattern[0] in "123456789"):
    return pattern[1:]+pattern[0]
  return None

class interpreter(object):

  def __init__(self,
        expected_patterns,
        synonym_patterns, non_hydrogens, hydrogens,
        mutually_exclusive_pairs=[]):
    expected_patterns_set = set()
    for expected_pattern in expected_patterns:
      expected_patterns_set.add(expected_pattern)
    assert len(expected_patterns_set) == len(expected_patterns)
    for synonym_pattern,expected_pattern in synonym_patterns.items():
      assert expected_pattern in expected_patterns_set
    for mep in mutually_exclusive_pairs:
      assert len(mep) == 3
      for expected_pattern in mep:
        if (expected_pattern not in expected_patterns_set):
          raise RuntimeError(
            "Inconsistent mutually_exclusive_pairs:\n"
            + "  given: %s\n" % str(mep)
            + "  pattern %s not in expected_patterns" % expected_pattern)
    expected = {}
    for expected_pattern in expected_patterns:
      assert expected_pattern.strip() == expected_pattern
      for h in ["H", "D"]:
        name = expected_pattern.replace("h", h)
        if (name in expected):
          raise RuntimeError(
            "Duplicate name %s given expected_patterns: %s" % (
              show_string(name), ", ".join([show_string(p)
                for p in expected_patterns])))
        expected[name] = expected_pattern
        if (name == expected_pattern): break
    synonym_patterns = dict(synonym_patterns)
    for expected_pattern in expected_patterns:
      alt = alternative_hydrogen_pattern(expected_pattern)
      if (alt is not None): synonym_patterns[alt] = expected_pattern
    for synonym_pattern,expected_pattern in list(synonym_patterns.items()):
      alt = alternative_hydrogen_pattern(synonym_pattern)
      if (alt is not None): synonym_patterns[alt] = expected_pattern
    synonyms = {}
    for synonym_pattern,expected_pattern in list(synonym_patterns.items()):
      for h in ["H", "D"]:
        name = synonym_pattern.replace("h", h)
        if (name in synonyms):
          raise RuntimeError(
            "Inconsistent synonym_patterns:\n"
            + "  synonym_pattern: %s\n" % synonym_pattern
            + "  name derived from synonym_pattern: %s\n" % name
            + "  Another synonym_pattern already lead to the same name.")
        synonyms[name] = expected_pattern.replace("h", h)
        if (name == synonym_pattern): break
    self.expected_patterns = expected_patterns
    self.synonym_patterns = synonym_patterns
    self.mutually_exclusive_pairs = mutually_exclusive_pairs
    self.expected = expected
    self.synonyms = synonyms
    self.non_hydrogens = non_hydrogens
    self.hydrogens = hydrogens

  def match_atom_names(self, atom_names):
    expected = {}
    unexpected = []
    for atom_name in atom_names:
      name = atom_name.strip().upper()
      expected_pattern = self.expected.get(self.synonyms.get(name, name))
      if (expected_pattern is None):
        unexpected.append(atom_name)
      else:
        expected.setdefault(expected_pattern, []).append(atom_name)
    return matched_atom_names(
      interpreter=self,
      atom_names=atom_names,
      expected=expected,
      unexpected=unexpected)

class matched_atom_names(object):

  def __init__(self, interpreter, atom_names, expected, unexpected):
    self.interpreter = interpreter
    self.atom_names = atom_names
    self.expected = expected
    self.unexpected = unexpected

  def __repr__(self):
    outl = "\n%s" % self.__class__.__name__
    for attr in ["atom_names",
                 "expected",
                 "unexpected",
                 ]:
      outl += "\n  %s" % attr
      for name in getattr(self, attr, []):
        outl += " %s" % name
    return outl

  def expected_patterns_with_multiple_matches(self):
    result = {}
    for expected_pattern,names in self.expected.items():
      if (len(names) != 1):
        result[expected_pattern] = names
    return result

  def mutually_exclusive_pairs(self):
    result = []
    for mep in self.interpreter.mutually_exclusive_pairs:
      if (    mep[0] in self.expected
          and mep[2] in self.expected):
        result.append((mep[0], mep[2]))
    return result

  def show_problems(self, out=None, prefix=""):
    result = 0
    if (out is None): out = sys.stdout
    if (len(self.unexpected) != 0):
      print(prefix+"unexpected atom names:", ", ".join(['"'+name+'"'
        for name in self.unexpected]), file=out)
      result += 1
    for expected_pattern,names in \
          self.expected_patterns_with_multiple_matches().items():
      print(prefix+"multiple matches: expected pattern=%s  names=%s" \
        % (expected_pattern, ", ".join(['"'+name+'"' for name in names])), file=out)
      result += 1
    for pair in self.mutually_exclusive_pairs():
      print(prefix+"mutually exclusive: %s" % " ".join(pair), file=out)
      result += 1
    return result

  def mon_lib_names(self):
    result = [None] * len(self.atom_names)
    name_indices = {}
    for i,name in enumerate(self.atom_names):
      name_indices.setdefault(name, []).append(i)
    mep_transl = {}
    for mep in self.interpreter.mutually_exclusive_pairs:
      if (mep[2] in self.expected):
        mep_transl[mep[1]] = mep[0]
        mep_transl[mep[2]] = mep[1]
      else:
        mep_transl[mep[0]] = mep[0]
        mep_transl[mep[1]] = mep[1]
    for expected_pattern,names in self.expected.items():
      expected_pattern = mep_transl.get(expected_pattern, expected_pattern)
      mon_lib_name = expected_pattern.upper()
      if (mon_lib_name[0] in "123456789"):
        mon_lib_name = mon_lib_name[1:] + mon_lib_name[0]
      for name in names:
        for i in name_indices[name]:
          result[i] = mon_lib_name
    return result

  def missing_atom_names(self, ignore_hydrogen=False):
    if ignore_hydrogen:
      return set(self.interpreter.non_hydrogens).difference(
        set(self.mon_lib_names()))
    else:
      return set(self.interpreter.non_hydrogens +
        self.interpreter.hydrogens).difference(
        set(self.mon_lib_names()))

peptide_expected_patterns = [
  "N", "h", "1h", "2h", "3h",
  "CA",
  "C", "O",
  "OXT", "hXT"]

peptide_synonym_patterns = dict_with_add({
  "OT1": "O",
  "OT2": "OXT",
  "OC":  "OXT",
  "hC":  "hXT",
  "hN": "h",
  "1hN": "1h",
  "2hN": "2h",
  "3hN": "3h",
  "1hT": "1h",
  "2hT": "2h",
  "3hT": "3h",
  "h0A": "1h",
  "h0B": "2h",
  "h0C": "3h"})

gly_interpreter = interpreter(
  peptide_expected_patterns + [
    "1hA", "2hA", "3hA"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hA", "2hA", "3hA")],
  non_hydrogens=("N","CA","C","O"),
  hydrogens=("H","HA1","HA2"))

ala_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB"],
  peptide_synonym_patterns,
  non_hydrogens=("N","CA","C","O","CB"),
  hydrogens=("H","HA","HB1","HB2","HB3"))

val_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "hB",
    "CG1", "1hG1", "2hG1", "3hG1",
    "CG2", "1hG2", "2hG2", "3hG2"],
  peptide_synonym_patterns,
  non_hydrogens=("N","CA","C","O","CB","CG1","CG2"),
  hydrogens=("H","HA","HB","HG11","HG12","HG13","HG21","HG22","HG23"))

leu_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG", "hG",
    "CD1", "1hD1", "2hD1", "3hD1",
    "CD2", "1hD2", "2hD2", "3hD2"],
  peptide_synonym_patterns + {"1hG": "hG"},
  mutually_exclusive_pairs=[
    ("1hB", "2hB", "3hB")],
  non_hydrogens=("N","CA","C","O","CB","CG","CD1","CD2"),
  hydrogens=("H","HA","HB1","HB2","HD11","HD12","HD13","HD21","HD22","HD23"))

ile_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "hB",
    "CG1", "1hG1", "2hG1", "3hG1",
    "CG2", "1hG2", "2hG2", "3hG2",
    "CD1", "1hD1", "2hD1", "3hD1"],
  peptide_synonym_patterns + {
    "CD": "CD1",
    "1hD": "1hD1",
    "2hD": "2hD1",
    "3hD": "3hD1"},
  mutually_exclusive_pairs=[
    ("1hG1", "2hG1", "3hG1")],
  non_hydrogens=("N","CA","C","O","CB","CG1","CD1","CG2"),
  hydrogens=("H","HA","HB","HG11","HG12","HD11","HD12","HD13","HG21","HG22",
    "HG23"))

met_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG", "1hG", "2hG", "3hG",
    "SD",
    "CE", "1hE", "2hE", "3hE"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "2hB", "3hB"),
    ("1hG", "2hG", "3hG")],
  non_hydrogens=("N","CA","C","O","CB","CG","SD","CE"),
  hydrogens=("H","HA","HB1","HB2","HG1","HG2","HE1","HE2","HE3"))

mse_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG", "1hG", "2hG", "3hG",
    "SE",
    "CE", "1hE", "2hE", "3hE"],
  peptide_synonym_patterns + {"SED": "SE"},
  mutually_exclusive_pairs=[
    ("1hB", "2hB", "3hB"),
    ("1hG", "2hG", "3hG")],
  non_hydrogens=("N","CA","C","O","CB","CG","SE","CE"),
  hydrogens=("H","HA","HB1","HB2","HG1","HG2","HE1","HE2","HE3"))

pro_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG", "1hG", "2hG", "3hG",
    "CD", "1hD", "2hD", "3hD"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "2hB", "3hB"),
    ("1hG", "2hG", "3hG"),
    ("1hD", "2hD", "3hD")],
  non_hydrogens=("N","CA","C","O","CB","CG","CD"),
  hydrogens=("HA","HB1","HB2","HG1","HG2","HD1","HD2"))

phe_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG",
    "CD1", "hD1",
    "CD2", "hD2",
    "CE1", "hE1",
    "CE2", "hE2",
    "CZ", "hZ"],
  peptide_synonym_patterns + {"1hZ": "hZ"},
  mutually_exclusive_pairs=[
    ("1hB", "2hB", "3hB")],
  non_hydrogens=("N","CA","C","O","CB","CG","CD1","CE1","CZ","CE2","CD2"),
  hydrogens=("H","HA","HB1","HB2","HD1","HE1","HZ","HE2","HD2"))

trp_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG",
    "CD1", "hD1",
    "CD2",
    "NE1", "hE1",
    "CE2",
    "CE3", "hE3",
    "CZ2", "hZ2",
    "CZ3", "hZ3",
    "CH2", "hH2"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "2hB", "3hB")],
  non_hydrogens=("N","CA","C","O","CB","CG","CD1","CD2","NE1","CE2","CE3","CZ2",
    "CZ3","CH2"),
  hydrogens=("H","HA","HB2","HB3","HD1","HE1","HE3","HZ2","HZ3","HH2"))

ser_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "OG", "hG"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "2hB", "3hB")],
  non_hydrogens=("N","CA","C","O","CB","OG"),
  hydrogens=("H","HA","HB1","HB2","HG"))

thr_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "hB",
    "OG1", "hG1",
    "CG2", "1hG2", "2hG2", "3hG2"],
  peptide_synonym_patterns,
  non_hydrogens=("N","CA","C","O","CB","OG1","CG2"),
  hydrogens=("H","HA","HB","HG1","HG21","HG22","HG23"))

asn_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG",
    "OD1",
    "ND2", "1hD2", "2hD2"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "2hB", "3hB")],
  non_hydrogens=("N","CA","C","O","CB","CG","OD1","ND2"),
  hydrogens=("H","HA","HB1","HB2","HD21","HD22"))

gln_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG", "1hG", "2hG", "3hG",
    "CD",
    "OE1",
    "NE2", "1hE2", "2hE2"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "2hB", "3hB"),
    ("1hG", "2hG", "3hG")],
  non_hydrogens=("N","CA","C","O","CB","CG","CD","OE1","NE2"),
  hydrogens=("H","HA","HB1","HB2","HG1","HG2","HE21","HE22"))

tyr_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB",
    "1hB", # needs to be commented for v3
    "2hB", "3hB",
    "CG",
    "CD1", "hD1",
    "CD2", "hD2",
    "CE1", "hE1",
    "CE2", "hE2",
    "CZ",
    "OH", "hH"],
  peptide_synonym_patterns,
  # needs to be commented for v3
  mutually_exclusive_pairs=[
    ("1hB", "2hB", "3hB")],
  non_hydrogens=("N","CA","C","O","CB","CG","CD1","CE1","CZ","OH","CE2","CD2"),
  hydrogens=("H","HA","HB1","HB2","HD1","HE1","HH","HE2","HD2"))

cys_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "SG", "hG"],
  peptide_synonym_patterns + {"1hG": "hG"},
  mutually_exclusive_pairs=[
    ("1hB", "2hB", "3hB")],
  non_hydrogens=("N","CA","C","O","CB","SG"),
  hydrogens=("H","HA","HB1","HB2","HG"))

lys_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG", "1hG", "2hG", "3hG",
    "CD", "1hD", "2hD", "3hD",
    "CE", "1hE", "2hE", "3hE",
    "NZ", "1hZ", "2hZ", "3hZ"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "2hB", "3hB"),
    ("1hG", "2hG", "3hG"),
    ("1hD", "2hD", "3hD"),
    ("1hE", "2hE", "3hE")],
  non_hydrogens=("N","CA","C","O","CB","CG","CD","CE","NZ"),
  hydrogens=("H","HA","HB1","HB2","HG1","HG2","HD1","HD2","HE1","HE2","HZ1",
    "HZ2","HZ3"))

arg_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG", "1hG", "2hG", "3hG",
    "CD", "1hD", "2hD", "3hD",
    "NE", "hE",
    "CZ",
    "NH1", "1hH1", "2hH1",
    "NH2", "1hH2", "2hH2",
    ],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "2hB", "3hB"),
    ("1hG", "2hG", "3hG"),
    ("1hD", "2hD", "3hD")],
  non_hydrogens=("N","CA","C","O","CB","CG","CD","NE","CZ","NH1","NH2"),
  hydrogens=("H","HA","HB1","HB2","HG1","HG2","HD1","HD2","HE","HH11","HH12",
    "HH21","HH22"))

his_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG",
    "ND1", "hD1",
    "CD2", "hD2",
    "CE1", "hE1",
    "NE2", "hE2"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "2hB", "3hB")],
  non_hydrogens=("N","CA","C","O","CB","CG","ND1","CE1","NE2","CD2"),
  hydrogens=("H","HA","HB1","HB2","HD1","HE1","HE2","HD2"))

asp_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG",
    "OD1", "hD1",
    "OD2", "hD2"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "2hB", "3hB"),
    ("hD1", "hD2", "hD2")],
  non_hydrogens=("N","CA","C","O","CB","CG","OD1","OD2"),
  hydrogens=("H","HA","HB1","HB2"))

glu_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG", "1hG", "2hG", "3hG",
    "CD",
    "OE1", "hE1",
    "OE2", "hE2"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "2hB", "3hB"),
    ("1hG", "2hG", "3hG")],
  non_hydrogens=("N","CA","C","O","CB","CG","CD","OE1","OE2"),
  hydrogens=("H","HA","HB1","HB2","HG1","HG2"))

interpreters = {
  "GLY": gly_interpreter,
  "ALA": ala_interpreter,
  "VAL": val_interpreter,
  "LEU": leu_interpreter,
  "ILE": ile_interpreter,
  "MET": met_interpreter,
  "MSE": mse_interpreter,
  "PRO": pro_interpreter,
  "PHE": phe_interpreter,
  "TRP": trp_interpreter,
  "SER": ser_interpreter,
  "THR": thr_interpreter,
  "ASN": asn_interpreter,
  "GLN": gln_interpreter,
  "TYR": tyr_interpreter,
  "CYS": cys_interpreter,
  "LYS": lys_interpreter,
  "ARG": arg_interpreter,
  "HIS": his_interpreter,
  "ASP": asp_interpreter,
  "GLU": glu_interpreter,
}
