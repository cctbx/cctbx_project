# -*- coding: utf-8; py-indent-offset: 2 -*-
from __future__ import division
from libtbx import group_args
from libtbx.utils import Sorry
from math import exp
import os

# Brown, I. D., & Altermatt, D. (1985).
# Bond-valence parameters obtained from a systematic analysis
# of the Inorganic Crystal Structure Database.
# Acta Crystallographica Section B Structural Science,
# 41(4), 244-247. doi:10.1107/S0108768185002063

def cif_param_as_list (param) :
  if (param == ".") : return None
  return param.split(",")

def cif_param_as_int (param) :
  if (param == ".") : return None
  return int(param)

def cif_param_as_float (param) :
  if (param == ".") : return None
  return float(param)

class MetalParameters (group_args) :
  def __str__ (self) :
    return "%s%+d" % (self.element.upper(), self.charge)

  def charge_as_int (self) :
    return self.charge

  def scattering_type (self) :
    charge_symbol = ""
    if (self.charge > 0) :
      charge_symbol = "+"
    elif (self.charge < 0) :
      charge_symbol = "-"
    s = "%2s%1d%s" % (self.element.strip(), abs(self.charge), charge_symbol)
    return s

class server (object) :
  def __init__ (self) :
    import iotbx.cif
    params_path = os.path.join(os.path.split(__file__)[0],
      "ion_parameters.cif")
    assert (os.path.isfile(params_path))
    cif_model = iotbx.cif.reader(file_path=params_path).model()
    self.params = cif_model["ions"]
    self._metal_params = {}
    self._charge_parms = {}
    self._resname_elem = {}
    self._default_charges = {}

  def is_supported_element (self, symbol) :
    return (symbol in self.params['_lib_valence.atom_symbol'])

  def is_supported_donor (self, symbol) :
    return (symbol in self.params['_lib_valence.donor_symbol'])

  def get_valence_params (self, atom1, atom2) :
    for i_elem, symbol in enumerate(self.params['_lib_valence.atom_symbol']) :
      if (symbol == atom1.element) :
        i_charge = int(self.params['_lib_valence.atom_charge'][i_elem])
        i_other = self.params['_lib_valence.donor_symbol'][i_elem]
        i_other_charge = int(self.params['_lib_valence.donor_charge'][i_elem])
        if ((i_charge == atom1.charge_as_int()) and
            (i_other == atom2.element) and
            (i_other_charge == atom2.charge_as_int())) :
          valence = float(self.params['_lib_valence.value'][i_elem])
          return valence, 0.37
    charge1 = atom1.charge_as_int()
    charge2 = atom2.charge_as_int()
    return None, None

  def _get_default_charge(self, element):
    if element in self._default_charges:
      return self._default_charges[element]
    p = self.params
    for i_elem, elem in enumerate(p["_lib_charge.element"]):
      if elem == element:
        charge = int(p["_lib_charge.charge"][i_elem])
        self._default_charges[element] = charge
        return charge
    return 0

  def _get_charge_params(self, resname, element = None):
    resname = resname.strip().upper()
    if element is not None: element = element.strip().upper()

    p = self.params

    if element is None:
      # Determine the element from the residue name (I.E. "HOH" => "O")
      if resname in self._resname_elem:
        element = self._resname_elem[resname]
      else:
        resn_elements = [(resn, p["_lib_charge.element"][i_resn])
                         for i_resn, resn in enumerate(p["_lib_charge.resname"])
                         if resn == resname]
        if len(resn_elements) > 1:
          raise Sorry("Ambiguous element for residue: " + resname)
        elif len(resn_elements) < 1:
          raise Sorry("Unknown element for residue: " + resname)
        element = resn_elements[0][1]
        self._resname_elem[resname] = element

    if (resname, element) in self._charge_parms:
      return self._charge_parms[(resname, element)]
    for i_resn, resn in enumerate(p["_lib_charge.resname"]):
      if resn == resname and element == p["_lib_charge.element"][i_resn]:
        elem_charge = \
          p["_lib_charge.element"][i_resn], int(p["_lib_charge.charge"][i_resn])
        break
    else:
      elem_charge = element, self._get_default_charge(element)

    self._charge_parms[(resname, element)] = elem_charge
    return elem_charge

  def get_element(self, atom):
    if isinstance(atom, str):
      resname = atom.strip().upper()
      if resname in self.params["_lib_charge.element"]:
        return resname
    else:
      if hasattr(atom, "element") and isinstance(atom.element, str):
        return atom.element.strip().upper()

      resname = atom.fetch_labels().resname.strip().upper()

    return self._get_charge_params(resname = resname)[0]

  def get_charge(self, atom):
    if isinstance(atom, str):
      atom = atom.strip().upper()
      try:
        charge = \
          self._get_charge_params(resname = atom)[1]
      except Sorry:
        charge = \
          self._get_charge_params(resname = "", element = atom)[1]
    else:
      charge = atom.charge
      if not isinstance(charge, int):
        charge = atom.charge_as_int()
      if charge != 0:
        return charge
      resname = atom.fetch_labels().resname.strip().upper()
      element = atom.element.strip().upper()
      charge = self._get_charge_params(resname = resname, element = element)[1]

    return charge

  def get_charges(self, atom):
    """
    Return all charges associated with element in our table.
    """
    element = self.get_element(atom)
    p = self.params
    charges = []
    for i_elem, elem in enumerate(p["_lib_charge.element"]):
      if elem == element:
        charges.append(int(p["_lib_charge.charge"][i_elem]))
    return charges

  def get_metal_parameters (self, element) :
    p = self.params
    for i_elem, symbol in enumerate(p['_lib_elems.element']) :
      if (symbol == element.upper()) :
        if (symbol in self._metal_params) :
          return self._metal_params[symbol]
        assert (p['_lib_ligands.element'][i_elem] == symbol)
        params = MetalParameters(
          element=symbol,
          charge=cif_param_as_int(p['_lib_elems.charge'][i_elem]),
          vec_sum_cutoff=cif_param_as_float(
            p["_lib_elems.vec_sum_cutoff"][i_elem]),
          coord_num_lower=cif_param_as_int(
            p["_lib_elems.coord_num_lower"][i_elem]),
          coord_num_upper=cif_param_as_int(
            p["_lib_elems.coord_num_upper"][i_elem]),
          min_coordinating_non_waters=cif_param_as_int(
            p["_lib_elems.min_coordinating_non_waters"][i_elem]),
          cvbs_lower=cif_param_as_float(p['_lib_elems.cvbs_lower'][i_elem]),
          cvbs_upper=cif_param_as_float(p['_lib_elems.cvbs_upper'][i_elem]),
          cvbs_expected=cif_param_as_float(
            p['_lib_elems.cvbs_expected'][i_elem]),
          allowed_coordinating_atoms=cif_param_as_list(
            p['_lib_ligands.allowed_coordinating_atoms'][i_elem]),
          allowed_coordinating_residues=cif_param_as_list(
            p['_lib_ligands.allowed_coordinating_residues'][i_elem]),
          allowed_geometries=cif_param_as_list(
            p['_lib_ligands.allowed_geometries'][i_elem]),
          allowed_backbone_atoms=cif_param_as_list(
            p['_lib_ligands.allowed_backbone_atoms'][i_elem]))
        self._metal_params[symbol] = params
        return params
    return None

  def calculate_valence (self, ion, donor, distance):
    """
    Calculates the single valence contribution of one ion donor pair,
    separated by distance. ion and donor should be AtomGuess objects.
    """

    if (not self.is_supported_donor(donor.element)) :
      return 0

    r_0, b = self.get_valence_params(ion, donor)

    if (r_0 is None) :
      # Try again, this time using the default charge for the donor
      donor = MetalParameters(
        charge = self.get_charge(donor.element),
        element = donor.element,
        )

      r_0, b = self.get_valence_params(ion, donor)

      if r_0 is None:
        return 0

    return exp((r_0 - distance) / b)

  def calculate_valences (self, ion, nearby_atoms):
    """
    Calculates all of the valence contributions between ion and each
    atom of nearby_atoms, each element of which should be a tuple of an
    atom and a vector from the ion's location.

    Returns a list of vectors, whose magnitudes are equal to the valence
    contributions from each donor atom.
    """

    vectors = []

    for contact in nearby_atoms:
      donor = MetalParameters(
        element = contact.element,
        charge = contact.charge
        )
      distance = abs(contact.vector)
      valence = self.calculate_valence(ion, donor, distance) * contact.occ

      if valence == 0:
        if ((donor.element not in ["H", "C", "AX"]) and
            (not self.is_supported_donor(donor.element))) :
          pass
      elif distance != 0:
        vectors.append(contact.vector / distance * valence)

    return vectors
