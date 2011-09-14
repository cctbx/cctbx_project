from libtbx import metric_prefixes

boltzmann_constant_j_per_k = 1.3806504e-23 # 1.3806504(24) 2006 CODATA value
avogadro_constant_per_mol = 6.022045e23 # 6.02214179(30) 2006 CODATA value
kcal_j = 4184.0 # exact value
amu_kg = 1.660538782e-27 # 1.660538782(83) 2006 CODATA value
angstroms_m = 1e-10 # exact value

akma_system_of_units = """\
AKMA = Angstroms, Kilocalories/Mol, Atomic mass units (as used e.g. by
CHARMM). Distances are measured in Angstroms, energies in kcal/mol,
mass in atomic mass units."""

boltzmann_constant_akma = avogadro_constant_per_mol \
                        * boltzmann_constant_j_per_k / kcal_j
akma_time_as_pico_seconds = ( \
  angstroms_m**2 * amu_kg * avogadro_constant_per_mol / kcal_j)**0.5 \
    / metric_prefixes.pico
