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
        synonym_patterns,
        mutually_exclusive_pairs=[]):
    expected_patterns_dict = {}
    for expected_pattern in expected_patterns:
      expected_patterns_dict[expected_pattern] = None
    assert len(expected_patterns_dict) == len(expected_patterns)
    for synonym_pattern,expected_pattern in synonym_patterns.items():
      assert expected_pattern in expected_patterns_dict
    for pair in mutually_exclusive_pairs:
      for expected_pattern in pair:
        if (expected_pattern not in expected_patterns_dict):
          raise RuntimeError(
            "Inconsistent mutually_exclusive_pairs:\n"
            + "  pair: %s\n" % str(pair)
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
    for synonym_pattern,expected_pattern in synonym_patterns.items():
      alt = alternative_hydrogen_pattern(synonym_pattern)
      if (alt is not None): synonym_patterns[alt] = expected_pattern
    synonyms = {}
    for synonym_pattern,expected_pattern in synonym_patterns.items():
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
      expected=expected,
      unexpected=unexpected)

class matched_atom_names(object):

  def __init__(self, interpreter, expected, unexpected):
    self.interpreter = interpreter
    self.expected = expected
    self.unexpected = unexpected

  def expected_patterns_with_multiple_matches(self):
    result = {}
    for expected_pattern,names in self.expected.items():
      if (len(names) != 1):
        result[expected_pattern] = names
    return result

  def mutually_exclusive_pairs(self):
    result = []
    for pair in self.interpreter.mutually_exclusive_pairs:
      if (    pair[0] in self.expected
          and pair[1] in self.expected):
        result.append(pair)
    return result

  def show_problems(self, out=None, prefix=""):
    result = 0
    if (out is None): out = sys.stdout
    if (len(self.unexpected) != 0):
      print >> out, prefix+"unexpected atom names:", ", ".join(['"'+name+'"'
        for name in self.unexpected])
      result += 1
    for expected_pattern,names in \
          self.expected_patterns_with_multiple_matches().items():
      print >> out, prefix+"multiple matches: expected pattern=%s  names=%s" \
        % (expected_pattern, ", ".join(['"'+name+'"' for name in names]))
      result += 1
    for pair in self.mutually_exclusive_pairs():
      print >> out, prefix+"mutually exclusive: %s" % " ".join(pair)
      result += 1
    return result

peptide_expected_patterns = [
  "N", "h", "1h", "2h", "3h",
  "CA",
  "C", "O",
  "OXT", "hXT"]

peptide_synonym_patterns = dict_with_add({
  "OT1": "O",
  "OT2": "OXT",
  "1hN": "1h",
  "2hN": "2h",
  "3hN": "3h"})

gly_interpreter = interpreter(
  peptide_expected_patterns + [
    "1hA", "2hA", "3hA"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hA", "3hA")])

ala_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB"],
  peptide_synonym_patterns)

val_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "hB",
    "CG1", "1hG1", "2hG1", "3hG1",
    "CG2", "1hG2", "2hG2", "3hG2"],
  peptide_synonym_patterns)

leu_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG", "hG",
    "CD1", "1hD1", "2hD1", "3hD1",
    "CD2", "1hD2", "2hD2", "3hD2"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "3hB")])

ile_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "hB",
    "CG1", "1hG1", "2hG1", "3hG1",
    "CG2", "1hG2", "2hG2", "3hG2",
    "CD1", "1hD1", "2hD1", "3hD1"],
  peptide_synonym_patterns + {"CD": "CD1"},
  mutually_exclusive_pairs=[
    ("1hG1", "3hG1")])

met_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG", "1hG", "2hG", "3hG",
    "SD",
    "CE", "1hE", "2hE", "3hE"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "3hB"),
    ("1hG", "3hG")])

mse_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB",
    "CG", "1hG", "2hG",
    "SE",
    "CE", "1hE", "2hE", "3hE"],
  peptide_synonym_patterns + {"SED": "SE"})

pro_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG", "1hG", "2hG", "3hG",
    "CD", "1hD", "2hD", "3hD"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "3hB"),
    ("1hG", "3hG"),
    ("1hD", "3hD")])

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
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "3hB")])

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
    ("1hB", "3hB")])

ser_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "OG", "hG"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "3hB")])

thr_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "hB",
    "OG1", "hG1",
    "CG2", "1hG2", "2hG2", "3hG2"],
  peptide_synonym_patterns)

asn_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG",
    "OD1",
    "ND2", "1hD2", "2hD2"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "3hB")])

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
    ("1hB", "3hB"),
    ("1hG", "3hG")])

tyr_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG",
    "CD1", "hD1",
    "CD2", "hD2",
    "CE1", "hE1",
    "CE2", "hE2",
    "CZ",
    "OH", "hH"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "3hB")])

cys_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "SG", "hG"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "3hB")])

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
    ("1hB", "3hB"),
    ("1hG", "3hG"),
    ("1hD", "3hD"),
    ("1hE", "3hE")])

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
    ("1hB", "3hB"),
    ("1hG", "3hG"),
    ("1hD", "3hD")])

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
    ("1hB", "3hB")])

asp_interpreter = interpreter(
  peptide_expected_patterns + [
    "hA",
    "CB", "1hB", "2hB", "3hB",
    "CG",
    "OD1", "hD1",
    "OD2", "hD2"],
  peptide_synonym_patterns,
  mutually_exclusive_pairs=[
    ("1hB", "3hB"),
    ("hD1", "hD2")])

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
    ("1hB", "3hB"),
    ("1hG", "3hG")])

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
