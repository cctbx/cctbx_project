from scitbx.python_utils.misc import sorted_store

class atom_info(sorted_store):

  def keys(self):
    return ("element_symbol", "ionic_state", "scattering_label")

# XXX This needs to be cleaned up at some point.

residues = {}

residues["OS4"] = {
"OS+4": atom_info("Os", 4, "Os4+"),
}
residues["ILE"] = {
"C": atom_info("C", 0, "C"),
"CD1": atom_info("C", 0, "C"),
"CB": atom_info("C", 0, "C"),
"CG1": atom_info("C", 0, "C"),
"H": atom_info("H", 0, "H'"),
"CA": atom_info("C", 0, "C"),
"CG2": atom_info("C", 0, "C"),
"O": atom_info("O", 0, "O"),
"N": atom_info("N", 0, "N"),
}
residues["NAG"] = {
"C8": atom_info("C", 0, "C"),
"H83": atom_info("H", 0, "H'"),
"H82": atom_info("H", 0, "H'"),
"H81": atom_info("H", 0, "H'"),
"C2": atom_info("C", 0, "C"),
"O7": atom_info("O", 0, "O"),
"O6": atom_info("O", 0, "O"),
"N2": atom_info("N", 0, "N"),
"O4": atom_info("O", 0, "O"),
"O3": atom_info("O", 0, "O"),
"O5": atom_info("O", 0, "O"),
"O1": atom_info("O", 0, "O"),
"H2": atom_info("H", 0, "H'"),
"H3": atom_info("H", 0, "H'"),
"H1": atom_info("H", 0, "H'"),
"H4": atom_info("H", 0, "H'"),
"H5": atom_info("H", 0, "H'"),
"HO6": atom_info("H", 0, "H'"),
"HO4": atom_info("H", 0, "H'"),
"HN2": atom_info("H", 0, "H'"),
"HO3": atom_info("H", 0, "H'"),
"HO1": atom_info("H", 0, "H'"),
"C3": atom_info("C", 0, "C"),
"H61": atom_info("H", 0, "H'"),
"C1": atom_info("C", 0, "C"),
"H62": atom_info("H", 0, "H'"),
"C7": atom_info("C", 0, "C"),
"C6": atom_info("C", 0, "C"),
"C5": atom_info("C", 0, "C"),
"C4": atom_info("C", 0, "C"),
}
residues["GLN"] = {
"C": atom_info("C", 0, "C"),
"CB": atom_info("C", 0, "C"),
"CA": atom_info("C", 0, "C"),
"CG": atom_info("C", 0, "C"),
"O": atom_info("O", 0, "O"),
"N": atom_info("N", 0, "N"),
"H": atom_info("H", 0, "H'"),
"HE22": atom_info("H", 0, "H'"),
"HE21": atom_info("H", 0, "H'"),
"CD": atom_info("C", 0, "C"),
"NE2": atom_info("N", 0, "N"),
"OE1": atom_info("O", 0, "O"),
}
residues["HG2"] = {
"HG+2": atom_info("Hg", 2, "Hg2+"),
}
residues["GLC"] = {
"C3": atom_info("C", 0, "C"),
"O6": atom_info("O", 0, "O"),
"O5": atom_info("O", 0, "O"),
"O4": atom_info("O", 0, "O"),
"O3": atom_info("O", 0, "O"),
"O2": atom_info("O", 0, "O"),
"O1": atom_info("O", 0, "O"),
"H2": atom_info("H", 0, "H'"),
"H3": atom_info("H", 0, "H'"),
"H1": atom_info("H", 0, "H'"),
"H4": atom_info("H", 0, "H'"),
"H5": atom_info("H", 0, "H'"),
"HO6": atom_info("H", 0, "H'"),
"HO4": atom_info("H", 0, "H'"),
"HO2": atom_info("H", 0, "H'"),
"HO3": atom_info("H", 0, "H'"),
"HO1": atom_info("H", 0, "H'"),
"H61": atom_info("H", 0, "H'"),
"C2": atom_info("C", 0, "C"),
"C1": atom_info("C", 0, "C"),
"H62": atom_info("H", 0, "H'"),
"C6": atom_info("C", 0, "C"),
"C5": atom_info("C", 0, "C"),
"C4": atom_info("C", 0, "C"),
}
residues["TIP"] = {
"H2": atom_info("H", 0, "H'"),
"H1": atom_info("H", 0, "H'"),
"OH2": atom_info("O", 0, "O"),
}
residues["BR1"] = {
"BR-1": atom_info("Br", -1, "Br1-"),
}
residues["GLY"] = {
"H": atom_info("H", 0, "H'"),
"CA": atom_info("C", 0, "C"),
"C": atom_info("C", 0, "C"),
"O": atom_info("O", 0, "O"),
"N": atom_info("N", 0, "N"),
}
residues["SIA"] = {
"C2": atom_info("C", 0, "C"),
"H5": atom_info("H", 0, "H'"),
"O9": atom_info("O", 0, "O"),
"C8": atom_info("C", 0, "C"),
"C9": atom_info("C", 0, "C"),
"H91": atom_info("H", 0, "H'"),
"H92": atom_info("H", 0, "H'"),
"O1B": atom_info("O", 0, "O"),
"O1A": atom_info("O", 0, "O"),
"O8": atom_info("O", 0, "O"),
"O7": atom_info("O", 0, "O"),
"O6": atom_info("O", 0, "O"),
"C11": atom_info("C", 0, "C"),
"C10": atom_info("C", 0, "C"),
"N5": atom_info("N", 0, "N"),
"H8": atom_info("H", 0, "H'"),
"O4": atom_info("O", 0, "O"),
"O10": atom_info("O", 0, "O"),
"H112": atom_info("H", 0, "H'"),
"H113": atom_info("H", 0, "H'"),
"H111": atom_info("H", 0, "H'"),
"H6": atom_info("H", 0, "H'"),
"H7": atom_info("H", 0, "H'"),
"H4": atom_info("H", 0, "H'"),
"O2": atom_info("O", 0, "O"),
"H32": atom_info("H", 0, "H'"),
"HO7": atom_info("H", 0, "H'"),
"HO4": atom_info("H", 0, "H'"),
"H31": atom_info("H", 0, "H'"),
"HO2": atom_info("H", 0, "H'"),
"HO8": atom_info("H", 0, "H'"),
"C3": atom_info("C", 0, "C"),
"HN5": atom_info("H", 0, "H'"),
"C1": atom_info("C", 0, "C"),
"HO9": atom_info("H", 0, "H'"),
"C7": atom_info("C", 0, "C"),
"C6": atom_info("C", 0, "C"),
"C5": atom_info("C", 0, "C"),
"C4": atom_info("C", 0, "C"),
}
residues["FE"] = {
"FE": atom_info("Fe", 0, "Fe"),
}
residues["BR"] = {
"BR": atom_info("Br", 0, "Br"),
}
residues["NI2"] = {
"NI+2": atom_info("Ni", 2, "Ni2+"),
}
residues["CPR"] = {
"C": atom_info("C", 0, "C"),
"CB": atom_info("C", 0, "C"),
"CA": atom_info("C", 0, "C"),
"N": atom_info("N", 0, "N"),
"CG": atom_info("C", 0, "C"),
"O": atom_info("O", 0, "O"),
"CD": atom_info("C", 0, "C"),
}
residues["HG1"] = {
"HG+1": atom_info("Hg", 1, "Hg1+"),
}
residues["GLU"] = {
"C": atom_info("C", 0, "C"),
"CB": atom_info("C", 0, "C"),
"OE2": atom_info("O", 0, "O"),
"H": atom_info("H", 0, "H'"),
"CA": atom_info("C", 0, "C"),
"N": atom_info("N", 0, "N"),
"CG": atom_info("C", 0, "C"),
"OE1": atom_info("O", 0, "O"),
"O": atom_info("O", 0, "O"),
"CD": atom_info("C", 0, "C"),
}
residues["NI"] = {
"NI": atom_info("Ni", 0, "Ni"),
}
residues["ASP"] = {
"C": atom_info("C", 0, "C"),
"O": atom_info("O", 0, "O"),
"OD1": atom_info("O", 0, "O"),
"H": atom_info("H", 0, "H'"),
"CA": atom_info("C", 0, "C"),
"CG": atom_info("C", 0, "C"),
"OD2": atom_info("O", 0, "O"),
"CB": atom_info("C", 0, "C"),
"N": atom_info("N", 0, "N"),
}
residues["CU2"] = {
"CU+2": atom_info("Cu", 2, "Cu2+"),
}
residues["CU1"] = {
"CU+1": atom_info("Cu", 1, "Cu1+"),
}
residues["LYS"] = {
"C": atom_info("C", 0, "C"),
"H": atom_info("H", 0, "H'"),
"CA": atom_info("C", 0, "C"),
"CG": atom_info("C", 0, "C"),
"CE": atom_info("C", 0, "C"),
"N": atom_info("N", 0, "N"),
"NZ": atom_info("N", 0, "N"),
"HZ1": atom_info("H", 0, "H'"),
"O": atom_info("O", 0, "O"),
"HZ3": atom_info("H", 0, "H'"),
"HZ2": atom_info("H", 0, "H'"),
"CD": atom_info("C", 0, "C"),
"CB": atom_info("C", 0, "C"),
}
residues["NA"] = {
"NA": atom_info("Na", 0, "Na"),
}
residues["PRO"] = {
"C": atom_info("C", 0, "C"),
"CB": atom_info("C", 0, "C"),
"CA": atom_info("C", 0, "C"),
"N": atom_info("N", 0, "N"),
"CG": atom_info("C", 0, "C"),
"O": atom_info("O", 0, "O"),
"CD": atom_info("C", 0, "C"),
}
residues["SR2"] = {
"SR+2": atom_info("Sr", 2, "Sr2+"),
}
residues["AU3"] = {
"AU+3": atom_info("Au", 3, "Au3+"),
}
residues["V3"] = {
"V+3": atom_info("V", 3, "V3+"),
}
residues["FE2"] = {
"FE+2": atom_info("Fe", 2, "Fe2+"),
}
residues["FE3"] = {
"FE+3": atom_info("Fe", 3, "Fe3+"),
}
residues["ASN"] = {
"C": atom_info("C", 0, "C"),
"HD22": atom_info("H", 0, "H'"),
"HD21": atom_info("H", 0, "H'"),
"CB": atom_info("C", 0, "C"),
"CA": atom_info("C", 0, "C"),
"CG": atom_info("C", 0, "C"),
"O": atom_info("O", 0, "O"),
"N": atom_info("N", 0, "N"),
"OD1": atom_info("O", 0, "O"),
"H": atom_info("H", 0, "H'"),
"ND2": atom_info("N", 0, "N"),
}
residues["PB"] = {
"PB": atom_info("Pb", 0, "Pb"),
}
residues["CR2"] = {
"CR+2": atom_info("Cr", 2, "Cr2+"),
}
residues["CO"] = {
"CO": atom_info("Co", 0, "Co"),
}
residues["THY"] = {
"O1P": atom_info("O", 0, "O"),
"C3'": atom_info("C", 0, "C"),
"C1'": atom_info("C", 0, "C"),
"O2P": atom_info("O", 0, "O"),
"C5'": atom_info("C", 0, "C"),
"P": atom_info("P", 0, "P"),
"C4'": atom_info("C", 0, "C"),
"C2'": atom_info("C", 0, "C"),
"O2'": atom_info("O", 0, "O"),
"N1": atom_info("N", 0, "N"),
"H2'": atom_info("H", 0, "H'"),
"O4": atom_info("O", 0, "O"),
"O2": atom_info("O", 0, "O"),
"N3": atom_info("N", 0, "N"),
"H3": atom_info("H", 0, "H'"),
"O5'": atom_info("O", 0, "O"),
"O4'": atom_info("O", 0, "O"),
"O3'": atom_info("O", 0, "O"),
"C5A": atom_info("C", 0, "C"),
"C2": atom_info("C", 0, "C"),
"C6": atom_info("C", 0, "C"),
"C5": atom_info("C", 0, "C"),
"C4": atom_info("C", 0, "C"),
}
residues["CL"] = {
"CL": atom_info("Cl", 0, "Cl"),
}
residues["CA"] = {
"CA": atom_info("Ca", 0, "Ca"),
}
residues["URI"] = {
"O1P": atom_info("O", 0, "O"),
"C3'": atom_info("C", 0, "C"),
"C1'": atom_info("C", 0, "C"),
"O2P": atom_info("O", 0, "O"),
"C5'": atom_info("C", 0, "C"),
"P": atom_info("P", 0, "P"),
"C4'": atom_info("C", 0, "C"),
"C2'": atom_info("C", 0, "C"),
"O2'": atom_info("O", 0, "O"),
"N1": atom_info("N", 0, "N"),
"H2'": atom_info("H", 0, "H'"),
"O4": atom_info("O", 0, "O"),
"O2": atom_info("O", 0, "O"),
"N3": atom_info("N", 0, "N"),
"H3": atom_info("H", 0, "H'"),
"O5'": atom_info("O", 0, "O"),
"O4'": atom_info("O", 0, "O"),
"C5": atom_info("C", 0, "C"),
"C2": atom_info("C", 0, "C"),
"C6": atom_info("C", 0, "C"),
"O3'": atom_info("O", 0, "O"),
"C4": atom_info("C", 0, "C"),
}
residues["XE"] = {
"XE": atom_info("Xe", 0, "Xe"),
}
residues["ADE"] = {
"O1P": atom_info("O", 0, "O"),
"C3'": atom_info("C", 0, "C"),
"C1'": atom_info("C", 0, "C"),
"O2P": atom_info("O", 0, "O"),
"C5'": atom_info("C", 0, "C"),
"P": atom_info("P", 0, "P"),
"C4'": atom_info("C", 0, "C"),
"C2'": atom_info("C", 0, "C"),
"C2": atom_info("C", 0, "C"),
"N1": atom_info("N", 0, "N"),
"H2'": atom_info("H", 0, "H'"),
"N3": atom_info("N", 0, "N"),
"N6": atom_info("N", 0, "N"),
"N7": atom_info("N", 0, "N"),
"N9": atom_info("N", 0, "N"),
"O5'": atom_info("O", 0, "O"),
"O4'": atom_info("O", 0, "O"),
"C8": atom_info("C", 0, "C"),
"O3'": atom_info("O", 0, "O"),
"H61": atom_info("H", 0, "H'"),
"O2'": atom_info("O", 0, "O"),
"H62": atom_info("H", 0, "H'"),
"C6": atom_info("C", 0, "C"),
"C5": atom_info("C", 0, "C"),
"C4": atom_info("C", 0, "C"),
}
residues["T"] = {
"O1P": atom_info("O", 0, "O"),
"C3'": atom_info("C", 0, "C"),
"C1'": atom_info("C", 0, "C"),
"O2P": atom_info("O", 0, "O"),
"C5'": atom_info("C", 0, "C"),
"P": atom_info("P", 0, "P"),
"C4'": atom_info("C", 0, "C"),
"C2'": atom_info("C", 0, "C"),
"O2'": atom_info("O", 0, "O"),
"N1": atom_info("N", 0, "N"),
"H2'": atom_info("H", 0, "H'"),
"O4": atom_info("O", 0, "O"),
"O2": atom_info("O", 0, "O"),
"N3": atom_info("N", 0, "N"),
"H3": atom_info("H", 0, "H'"),
"O5'": atom_info("O", 0, "O"),
"O4'": atom_info("O", 0, "O"),
"O3'": atom_info("O", 0, "O"),
"C5A": atom_info("C", 0, "C"),
"C2": atom_info("C", 0, "C"),
"C6": atom_info("C", 0, "C"),
"C5": atom_info("C", 0, "C"),
"C4": atom_info("C", 0, "C"),
}
residues["XYL"] = {
"H2": atom_info("H", 0, "H'"),
"H3": atom_info("H", 0, "H'"),
"H1": atom_info("H", 0, "H'"),
"H4": atom_info("H", 0, "H'"),
"HO3": atom_info("H", 0, "H'"),
"C5": atom_info("C", 0, "C"),
"HO4": atom_info("H", 0, "H'"),
"C1": atom_info("C", 0, "C"),
"HO2": atom_info("H", 0, "H'"),
"H51": atom_info("H", 0, "H'"),
"H52": atom_info("H", 0, "H'"),
"HO1": atom_info("H", 0, "H'"),
"C3": atom_info("C", 0, "C"),
"C2": atom_info("C", 0, "C"),
"O5": atom_info("O", 0, "O"),
"O4": atom_info("O", 0, "O"),
"O3": atom_info("O", 0, "O"),
"O2": atom_info("O", 0, "O"),
"O1": atom_info("O", 0, "O"),
"C4": atom_info("C", 0, "C"),
}
residues["V"] = {
"V": atom_info("V", 0, "V"),
}
residues["CS"] = {
"CS": atom_info("Cs", 0, "Cs"),
}
residues["CR"] = {
"CR": atom_info("Cr", 0, "Cr"),
}
residues["WAT"] = {
"H2": atom_info("H", 0, "H'"),
"H1": atom_info("H", 0, "H'"),
"OH2": atom_info("O", 0, "O"),
}
residues["CU"] = {
"CU": atom_info("Cu", 0, "Cu"),
}
residues["F1"] = {
"F-1": atom_info("F", -1, "F1-"),
}
residues["CO3"] = {
"CO+3": atom_info("Co", 3, "Co3+"),
}
residues["CO2"] = {
"CO+2": atom_info("Co", 2, "Co2+"),
}
residues["SR"] = {
"SR": atom_info("Sr", 0, "Sr"),
}
residues["K1"] = {
"K+1": atom_info("K", 1, "K1+"),
}
residues["PHE"] = {
"C": atom_info("C", 0, "C"),
"CD2": atom_info("C", 0, "C"),
"CB": atom_info("C", 0, "C"),
"CA": atom_info("C", 0, "C"),
"CG": atom_info("C", 0, "C"),
"O": atom_info("O", 0, "O"),
"N": atom_info("N", 0, "N"),
"CZ": atom_info("C", 0, "C"),
"CE2": atom_info("C", 0, "C"),
"CE1": atom_info("C", 0, "C"),
"H": atom_info("H", 0, "H'"),
"CD1": atom_info("C", 0, "C"),
}
residues["ALA"] = {
"C": atom_info("C", 0, "C"),
"CB": atom_info("C", 0, "C"),
"H": atom_info("H", 0, "H'"),
"CA": atom_info("C", 0, "C"),
"O": atom_info("O", 0, "O"),
"N": atom_info("N", 0, "N"),
}
residues["KR"] = {
"KR": atom_info("Kr", 0, "Kr"),
}
residues["NA1"] = {
"NA+1": atom_info("Na", 1, "Na1+"),
}
residues["CL1"] = {
"CL-1": atom_info("Cl", -1, "Cl1-"),
}
residues["LI1"] = {
"LI+1": atom_info("Li", 1, "Li1+"),
}
residues["YB2"] = {
"YB+2": atom_info("Yb", 2, "Yb2+"),
}
residues["OS"] = {
"OS": atom_info("Os", 0, "Os"),
}
residues["IR3"] = {
"IR+3": atom_info("Ir", 3, "Ir3+"),
}
residues["AG1"] = {
"AG+1": atom_info("Ag", 1, "Ag1+"),
}
residues["MAN"] = {
"C3": atom_info("C", 0, "C"),
"O6": atom_info("O", 0, "O"),
"O5": atom_info("O", 0, "O"),
"O4": atom_info("O", 0, "O"),
"O3": atom_info("O", 0, "O"),
"O2": atom_info("O", 0, "O"),
"O1": atom_info("O", 0, "O"),
"H2": atom_info("H", 0, "H'"),
"H3": atom_info("H", 0, "H'"),
"H1": atom_info("H", 0, "H'"),
"H4": atom_info("H", 0, "H'"),
"H5": atom_info("H", 0, "H'"),
"HO6": atom_info("H", 0, "H'"),
"HO4": atom_info("H", 0, "H'"),
"HO2": atom_info("H", 0, "H'"),
"HO3": atom_info("H", 0, "H'"),
"HO1": atom_info("H", 0, "H'"),
"H61": atom_info("H", 0, "H'"),
"C2": atom_info("C", 0, "C"),
"C1": atom_info("C", 0, "C"),
"H62": atom_info("H", 0, "H'"),
"C6": atom_info("C", 0, "C"),
"C5": atom_info("C", 0, "C"),
"C4": atom_info("C", 0, "C"),
}
residues["GUA"] = {
"O1P": atom_info("O", 0, "O"),
"H21": atom_info("H", 0, "H'"),
"C3'": atom_info("C", 0, "C"),
"C1'": atom_info("C", 0, "C"),
"O2P": atom_info("O", 0, "O"),
"H22": atom_info("H", 0, "H'"),
"C5'": atom_info("C", 0, "C"),
"P": atom_info("P", 0, "P"),
"C4'": atom_info("C", 0, "C"),
"C2'": atom_info("C", 0, "C"),
"N1": atom_info("N", 0, "N"),
"O6": atom_info("O", 0, "O"),
"N2": atom_info("N", 0, "N"),
"N3": atom_info("N", 0, "N"),
"H2'": atom_info("H", 0, "H'"),
"N7": atom_info("N", 0, "N"),
"N9": atom_info("N", 0, "N"),
"H1": atom_info("H", 0, "H'"),
"O5'": atom_info("O", 0, "O"),
"O4'": atom_info("O", 0, "O"),
"C8": atom_info("C", 0, "C"),
"O2'": atom_info("O", 0, "O"),
"C2": atom_info("C", 0, "C"),
"C4": atom_info("C", 0, "C"),
"C6": atom_info("C", 0, "C"),
"C5": atom_info("C", 0, "C"),
"O3'": atom_info("O", 0, "O"),
}
residues["MN"] = {
"MN": atom_info("Mn", 0, "Mn"),
}
residues["CD2"] = {
"CD+2": atom_info("Cd", 2, "Cd2+"),
}
residues["MO3"] = {
"MO+3": atom_info("Mo", 3, "Mo3+"),
}
residues["U4"] = {
"U+4": atom_info("U", 4, "U4+"),
}
residues["YB"] = {
"YB": atom_info("Yb", 0, "Yb"),
}
residues["ZN"] = {
"ZN": atom_info("Zn", 0, "Zn"),
}
residues["U3"] = {
"U+3": atom_info("U", 3, "U3+"),
}
residues["MG2"] = {
"MG+2": atom_info("Mg", 2, "Mg2+"),
}
residues["FUC"] = {
"H63": atom_info("H", 0, "H'"),
"C3": atom_info("C", 0, "C"),
"O5": atom_info("O", 0, "O"),
"O4": atom_info("O", 0, "O"),
"O3": atom_info("O", 0, "O"),
"O2": atom_info("O", 0, "O"),
"O1": atom_info("O", 0, "O"),
"H2": atom_info("H", 0, "H'"),
"H3": atom_info("H", 0, "H'"),
"H1": atom_info("H", 0, "H'"),
"H4": atom_info("H", 0, "H'"),
"H5": atom_info("H", 0, "H'"),
"HO4": atom_info("H", 0, "H'"),
"HO2": atom_info("H", 0, "H'"),
"HO3": atom_info("H", 0, "H'"),
"HO1": atom_info("H", 0, "H'"),
"H61": atom_info("H", 0, "H'"),
"C2": atom_info("C", 0, "C"),
"C1": atom_info("C", 0, "C"),
"H62": atom_info("H", 0, "H'"),
"C6": atom_info("C", 0, "C"),
"C5": atom_info("C", 0, "C"),
"C4": atom_info("C", 0, "C"),
}
residues["MSE"] = {
"C": atom_info("C", 0, "C"),
"CE": atom_info("C", 0, "C"),
"O": atom_info("O", 0, "O"),
"CB": atom_info("C", 0, "C"),
"H": atom_info("H", 0, "H'"),
"CA": atom_info("C", 0, "C"),
"CG": atom_info("C", 0, "C"),
"SE": atom_info("Se", 0, "Se"),
"N": atom_info("N", 0, "N"),
}
residues["HG"] = {
"HG": atom_info("Hg", 0, "Hg"),
}
residues["CMP"] = {
"O1P": atom_info("O", 0, "O"),
"C3'": atom_info("C", 0, "C"),
"C1'": atom_info("C", 0, "C"),
"O2P": atom_info("O", 0, "O"),
"C5'": atom_info("C", 0, "C"),
"P": atom_info("P", 0, "P"),
"C4'": atom_info("C", 0, "C"),
"C2'": atom_info("C", 0, "C"),
"C2": atom_info("C", 0, "C"),
"N1": atom_info("N", 0, "N"),
"H2'": atom_info("H", 0, "H'"),
"N3": atom_info("N", 0, "N"),
"N6": atom_info("N", 0, "N"),
"N7": atom_info("N", 0, "N"),
"N9": atom_info("N", 0, "N"),
"O5'": atom_info("O", 0, "O"),
"O4'": atom_info("O", 0, "O"),
"C8": atom_info("C", 0, "C"),
"O3'": atom_info("O", 0, "O"),
"H61": atom_info("H", 0, "H'"),
"O2'": atom_info("O", 0, "O"),
"H62": atom_info("H", 0, "H'"),
"C6": atom_info("C", 0, "C"),
"C5": atom_info("C", 0, "C"),
"C4": atom_info("C", 0, "C"),
}
residues["CYS"] = {
"C": atom_info("C", 0, "C"),
"CB": atom_info("C", 0, "C"),
"H": atom_info("H", 0, "H'"),
"CA": atom_info("C", 0, "C"),
"SG": atom_info("S", 0, "S"),
"O": atom_info("O", 0, "O"),
"N": atom_info("N", 0, "N"),
}
residues["HIS"] = {
"C": atom_info("C", 0, "C"),
"HE2": atom_info("H", 0, "H'"),
"CD2": atom_info("C", 0, "C"),
"HD1": atom_info("H", 0, "H'"),
"CB": atom_info("C", 0, "C"),
"CA": atom_info("C", 0, "C"),
"CG": atom_info("C", 0, "C"),
"O": atom_info("O", 0, "O"),
"N": atom_info("N", 0, "N"),
"CE1": atom_info("C", 0, "C"),
"H": atom_info("H", 0, "H'"),
"ND1": atom_info("N", 0, "N"),
"NE2": atom_info("N", 0, "N"),
}
residues["VAL"] = {
"C": atom_info("C", 0, "C"),
"CB": atom_info("C", 0, "C"),
"CG1": atom_info("C", 0, "C"),
"H": atom_info("H", 0, "H'"),
"CA": atom_info("C", 0, "C"),
"CG2": atom_info("C", 0, "C"),
"O": atom_info("O", 0, "O"),
"N": atom_info("N", 0, "N"),
}
residues["PT"] = {
"PT": atom_info("Pt", 0, "Pt"),
}
residues["CYT"] = {
"O1P": atom_info("O", 0, "O"),
"C3'": atom_info("C", 0, "C"),
"C1'": atom_info("C", 0, "C"),
"O2P": atom_info("O", 0, "O"),
"C5'": atom_info("C", 0, "C"),
"P": atom_info("P", 0, "P"),
"C4'": atom_info("C", 0, "C"),
"C2'": atom_info("C", 0, "C"),
"O2'": atom_info("O", 0, "O"),
"N1": atom_info("N", 0, "N"),
"H2'": atom_info("H", 0, "H'"),
"N3": atom_info("N", 0, "N"),
"N4": atom_info("N", 0, "N"),
"O2": atom_info("O", 0, "O"),
"O5'": atom_info("O", 0, "O"),
"O4'": atom_info("O", 0, "O"),
"O3'": atom_info("O", 0, "O"),
"H41": atom_info("H", 0, "H'"),
"C2": atom_info("C", 0, "C"),
"H42": atom_info("H", 0, "H'"),
"C6": atom_info("C", 0, "C"),
"C5": atom_info("C", 0, "C"),
"C4": atom_info("C", 0, "C"),
}
residues["I1"] = {
"I-1": atom_info("I", -1, "I1-"),
}
residues["WO4"] = {
"O4": atom_info("O", 0, "O"),
"O3": atom_info("O", 0, "O"),
"O1": atom_info("O", 0, "O"),
"W": atom_info("W", 0, "W"),
"O2": atom_info("O", 0, "O"),
}
residues["F"] = {
"F": atom_info("F", 0, "F"),
}
residues["HISH"] = {
"C": atom_info("C", 0, "C"),
"HE2": atom_info("H", 0, "H'"),
"CD2": atom_info("C", 0, "C"),
"HD1": atom_info("H", 0, "H'"),
"CB": atom_info("C", 0, "C"),
"CA": atom_info("C", 0, "C"),
"CG": atom_info("C", 0, "C"),
"O": atom_info("O", 0, "O"),
"N": atom_info("N", 0, "N"),
"CE1": atom_info("C", 0, "C"),
"H": atom_info("H", 0, "H'"),
"ND1": atom_info("N", 0, "N"),
"NE2": atom_info("N", 0, "N"),
}
residues["AL3"] = {
"AL+3": atom_info("Al", 3, "Al3+"),
}
residues["HO"] = {
"HO": atom_info("Ho", 0, "Ho"),
}
residues["V2"] = {
"V+2": atom_info("V", 2, "V2+"),
}
residues["GAL"] = {
"C3": atom_info("C", 0, "C"),
"O6": atom_info("O", 0, "O"),
"O5": atom_info("O", 0, "O"),
"O4": atom_info("O", 0, "O"),
"O3": atom_info("O", 0, "O"),
"O2": atom_info("O", 0, "O"),
"O1": atom_info("O", 0, "O"),
"H2": atom_info("H", 0, "H'"),
"H3": atom_info("H", 0, "H'"),
"H1": atom_info("H", 0, "H'"),
"H4": atom_info("H", 0, "H'"),
"H5": atom_info("H", 0, "H'"),
"HO6": atom_info("H", 0, "H'"),
"HO4": atom_info("H", 0, "H'"),
"HO2": atom_info("H", 0, "H'"),
"HO3": atom_info("H", 0, "H'"),
"HO1": atom_info("H", 0, "H'"),
"H61": atom_info("H", 0, "H'"),
"C2": atom_info("C", 0, "C"),
"C1": atom_info("C", 0, "C"),
"H62": atom_info("H", 0, "H'"),
"C6": atom_info("C", 0, "C"),
"C5": atom_info("C", 0, "C"),
"C4": atom_info("C", 0, "C"),
}
residues["HO3"] = {
"HO+3": atom_info("Ho", 3, "Ho3+"),
}
residues["CS1"] = {
"CS+1": atom_info("Cs", 1, "Cs1+"),
}
residues["PO4"] = {
"P": atom_info("P", 0, "P"),
"O4": atom_info("O", 0, "O"),
"O3": atom_info("O", 0, "O"),
"O2": atom_info("O", 0, "O"),
"O1": atom_info("O", 0, "O"),
}
residues["LEU"] = {
"C": atom_info("C", 0, "C"),
"CD1": atom_info("C", 0, "C"),
"CD2": atom_info("C", 0, "C"),
"CB": atom_info("C", 0, "C"),
"H": atom_info("H", 0, "H'"),
"CA": atom_info("C", 0, "C"),
"CG": atom_info("C", 0, "C"),
"O": atom_info("O", 0, "O"),
"N": atom_info("N", 0, "N"),
}
residues["A"] = {
"O1P": atom_info("O", 0, "O"),
"C3'": atom_info("C", 0, "C"),
"C1'": atom_info("C", 0, "C"),
"O2P": atom_info("O", 0, "O"),
"C5'": atom_info("C", 0, "C"),
"P": atom_info("P", 0, "P"),
"C4'": atom_info("C", 0, "C"),
"C2'": atom_info("C", 0, "C"),
"C2": atom_info("C", 0, "C"),
"N1": atom_info("N", 0, "N"),
"H2'": atom_info("H", 0, "H'"),
"N3": atom_info("N", 0, "N"),
"N6": atom_info("N", 0, "N"),
"N7": atom_info("N", 0, "N"),
"N9": atom_info("N", 0, "N"),
"O5'": atom_info("O", 0, "O"),
"O4'": atom_info("O", 0, "O"),
"C8": atom_info("C", 0, "C"),
"O3'": atom_info("O", 0, "O"),
"H61": atom_info("H", 0, "H'"),
"O2'": atom_info("O", 0, "O"),
"H62": atom_info("H", 0, "H'"),
"C6": atom_info("C", 0, "C"),
"C5": atom_info("C", 0, "C"),
"C4": atom_info("C", 0, "C"),
}
residues["C"] = {
"O1P": atom_info("O", 0, "O"),
"C3'": atom_info("C", 0, "C"),
"C1'": atom_info("C", 0, "C"),
"O2P": atom_info("O", 0, "O"),
"C5'": atom_info("C", 0, "C"),
"P": atom_info("P", 0, "P"),
"C4'": atom_info("C", 0, "C"),
"C2'": atom_info("C", 0, "C"),
"O2'": atom_info("O", 0, "O"),
"N1": atom_info("N", 0, "N"),
"H2'": atom_info("H", 0, "H'"),
"N3": atom_info("N", 0, "N"),
"N4": atom_info("N", 0, "N"),
"O2": atom_info("O", 0, "O"),
"O5'": atom_info("O", 0, "O"),
"O4'": atom_info("O", 0, "O"),
"O3'": atom_info("O", 0, "O"),
"H41": atom_info("H", 0, "H'"),
"C2": atom_info("C", 0, "C"),
"H42": atom_info("H", 0, "H'"),
"C6": atom_info("C", 0, "C"),
"C5": atom_info("C", 0, "C"),
"C4": atom_info("C", 0, "C"),
}
residues["MN3"] = {
"MN+3": atom_info("Mn", 3, "Mn3+"),
}
residues["G"] = {
"O1P": atom_info("O", 0, "O"),
"H21": atom_info("H", 0, "H'"),
"C3'": atom_info("C", 0, "C"),
"C1'": atom_info("C", 0, "C"),
"O2P": atom_info("O", 0, "O"),
"H22": atom_info("H", 0, "H'"),
"C5'": atom_info("C", 0, "C"),
"P": atom_info("P", 0, "P"),
"C4'": atom_info("C", 0, "C"),
"C2'": atom_info("C", 0, "C"),
"N1": atom_info("N", 0, "N"),
"O6": atom_info("O", 0, "O"),
"N2": atom_info("N", 0, "N"),
"N3": atom_info("N", 0, "N"),
"H2'": atom_info("H", 0, "H'"),
"N7": atom_info("N", 0, "N"),
"N9": atom_info("N", 0, "N"),
"H1": atom_info("H", 0, "H'"),
"O5'": atom_info("O", 0, "O"),
"O4'": atom_info("O", 0, "O"),
"C8": atom_info("C", 0, "C"),
"O2'": atom_info("O", 0, "O"),
"C2": atom_info("C", 0, "C"),
"C4": atom_info("C", 0, "C"),
"C6": atom_info("C", 0, "C"),
"C5": atom_info("C", 0, "C"),
"O3'": atom_info("O", 0, "O"),
}
residues["I"] = {
"I": atom_info("I", 0, "I"),
}
residues["K"] = {
"K": atom_info("K", 0, "K"),
}
residues["IR"] = {
"IR": atom_info("Ir", 0, "Ir"),
}
residues["PB2"] = {
"PB+2": atom_info("Pb", 2, "Pb2+"),
}
residues["HOH"] = {
"H2": atom_info("H", 0, "H'"),
"H1": atom_info("H", 0, "H'"),
"O": atom_info("O", 0, "O"),
}
residues["THR"] = {
"C": atom_info("C", 0, "C"),
"CB": atom_info("C", 0, "C"),
"N": atom_info("N", 0, "N"),
"H": atom_info("H", 0, "H'"),
"CA": atom_info("C", 0, "C"),
"CG2": atom_info("C", 0, "C"),
"OG1": atom_info("O", 0, "O"),
"O": atom_info("O", 0, "O"),
"HG1": atom_info("H", 0, "H'"),
}
residues["AS"] = {
"AS": atom_info("As", 0, "As"),
}
residues["AR"] = {
"AR": atom_info("Ar", 0, "Ar"),
}
residues["U"] = {
"U": atom_info("U", 0, "U"),
}
residues["SO4"] = {
"S": atom_info("S", 0, "S"),
"O4": atom_info("O", 0, "O"),
"O3": atom_info("O", 0, "O"),
"O2": atom_info("O", 0, "O"),
"O1": atom_info("O", 0, "O"),
}
residues["ZN2"] = {
"ZN+2": atom_info("Zn", 2, "Zn2+"),
}
residues["TRP"] = {
"C": atom_info("C", 0, "C"),
"HE1": atom_info("H", 0, "H'"),
"CZ2": atom_info("C", 0, "C"),
"CZ3": atom_info("C", 0, "C"),
"CA": atom_info("C", 0, "C"),
"CG": atom_info("C", 0, "C"),
"O": atom_info("O", 0, "O"),
"N": atom_info("N", 0, "N"),
"CH2": atom_info("C", 0, "C"),
"CE3": atom_info("C", 0, "C"),
"CD1": atom_info("C", 0, "C"),
"CD2": atom_info("C", 0, "C"),
"H": atom_info("H", 0, "H'"),
"CE2": atom_info("C", 0, "C"),
"CB": atom_info("C", 0, "C"),
"NE1": atom_info("N", 0, "N"),
}
residues["TIP3"] = {
"H2": atom_info("H", 0, "H'"),
"H1": atom_info("H", 0, "H'"),
"OH2": atom_info("O", 0, "O"),
}
residues["CD"] = {
"CD": atom_info("Cd", 0, "Cd"),
}
residues["PT2"] = {
"PT+2": atom_info("Pt", 2, "Pt2+"),
}
residues["SER"] = {
"C": atom_info("C", 0, "C"),
"CB": atom_info("C", 0, "C"),
"OG": atom_info("O", 0, "O"),
"H": atom_info("H", 0, "H'"),
"CA": atom_info("C", 0, "C"),
"HG": atom_info("H", 0, "H'"),
"O": atom_info("O", 0, "O"),
"N": atom_info("N", 0, "N"),
}
residues["AL"] = {
"AL": atom_info("Al", 0, "Al"),
}
residues["CR3"] = {
"CR+3": atom_info("Cr", 3, "Cr3+"),
}
residues["MG"] = {
"MG": atom_info("Mg", 0, "Mg"),
}
residues["MN2"] = {
"MN+2": atom_info("Mn", 2, "Mn2+"),
}
residues["CA2"] = {
"CA+2": atom_info("Ca", 2, "Ca2+"),
}
residues["MO"] = {
"MO": atom_info("Mo", 0, "Mo"),
}
residues["YB3"] = {
"YB+3": atom_info("Yb", 3, "Yb3+"),
}
residues["AU1"] = {
"AU+1": atom_info("Au", 1, "Au1+"),
}
residues["M3L"] = {
"C": atom_info("C", 0, "C"),
"CB": atom_info("C", 0, "C"),
"CA": atom_info("C", 0, "C"),
"CG": atom_info("C", 0, "C"),
"CE": atom_info("C", 0, "C"),
"N": atom_info("N", 0, "N"),
"NZ": atom_info("N", 0, "N"),
"O": atom_info("O", 0, "O"),
"H": atom_info("H", 0, "H'"),
"CD": atom_info("C", 0, "C"),
"CM3": atom_info("C", 0, "C"),
"CM2": atom_info("C", 0, "C"),
"CM1": atom_info("C", 0, "C"),
}
residues["TYR"] = {
"C": atom_info("C", 0, "C"),
"CE1": atom_info("C", 0, "C"),
"OH": atom_info("O", 0, "O"),
"CB": atom_info("C", 0, "C"),
"CA": atom_info("C", 0, "C"),
"CG": atom_info("C", 0, "C"),
"O": atom_info("O", 0, "O"),
"N": atom_info("N", 0, "N"),
"CZ": atom_info("C", 0, "C"),
"HH": atom_info("H", 0, "H'"),
"CE2": atom_info("C", 0, "C"),
"CD2": atom_info("C", 0, "C"),
"H": atom_info("H", 0, "H'"),
"CD1": atom_info("C", 0, "C"),
}
residues["ARG"] = {
"C": atom_info("C", 0, "C"),
"H": atom_info("H", 0, "H'"),
"CB": atom_info("C", 0, "C"),
"CA": atom_info("C", 0, "C"),
"CG": atom_info("C", 0, "C"),
"NE": atom_info("N", 0, "N"),
"O": atom_info("O", 0, "O"),
"N": atom_info("N", 0, "N"),
"HH22": atom_info("H", 0, "H'"),
"CZ": atom_info("C", 0, "C"),
"HH21": atom_info("H", 0, "H'"),
"NH1": atom_info("N", 0, "N"),
"HH12": atom_info("H", 0, "H'"),
"HH11": atom_info("H", 0, "H'"),
"NH2": atom_info("N", 0, "N"),
"CD": atom_info("C", 0, "C"),
"HE": atom_info("H", 0, "H'"),
}
residues["AG"] = {
"AG": atom_info("Ag", 0, "Ag"),
}
residues["AU"] = {
"AU": atom_info("Au", 0, "Au"),
}
residues["MET"] = {
"C": atom_info("C", 0, "C"),
"O": atom_info("O", 0, "O"),
"CB": atom_info("C", 0, "C"),
"H": atom_info("H", 0, "H'"),
"CA": atom_info("C", 0, "C"),
"N": atom_info("N", 0, "N"),
"CG": atom_info("C", 0, "C"),
"CE": atom_info("C", 0, "C"),
"SD": atom_info("S", 0, "S"),
}

fallback_dictionary = {
"OXT": atom_info("O", 0, "O"),
"OT1": atom_info("O", 0, "O"),
"OT2": atom_info("O", 0, "O"),
"OT": atom_info("O", 0, "O"),
"SE": atom_info("Se", 0, "Se"),
}

def get(residue_name, atom_name):
  try: return residues[residue_name][atom_name.replace("*", "'")]
  except KeyError: pass
  try: return fallback_dictionary[atom_name]
  except KeyError: pass
  raise KeyError, "Unknown residue name (%s) or atom name (%s)." % (
    residue_name, atom_name)

if (__name__ == "__main__"):
  get("GLN", "CA").show()
