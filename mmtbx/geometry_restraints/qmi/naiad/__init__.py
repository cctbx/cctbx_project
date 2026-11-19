"""naiad: H-bond-network-aware, cation-aware placement of water hydrogens.

Map-free placement of the two H on every bare HOH/DOD water: H point at real
H-bond acceptors, avoid clashes, and stay off metal cations, purely from
geometry (no experimental map, no residue tables, no monomer library). Used by
the qmi metal-restraint pipeline to protonate coordinating waters that
mmtbx.reduce2 leaves bare. Public entry point:
``water_protonation.place_water_hydrogens``.
"""
