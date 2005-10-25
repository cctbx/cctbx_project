# Scattering factors for Dummy Bond Electrons (dbe)
#...................................................
# Created: Pavel Afonine
#...................................................
# References:
#...................................................
# dbe currently available for: H3C-NH-CO-CH3 and Phe ring
#
#    Bond      dbe name
# [H3]C-N[H]     DNC
# [H]N-C[O]      DNCO
# C=O            DO2
# [O]C-C[H3]     DCC
# [H2]C-H        DCH3
# N-H            DNH
#

class one_gaussian_dbe(object):
  source="Afonine P.V et al. (2004). Acta Cryst., D60"
  source_short = "Afonine P.V et al. (2004). Acta Cryst., D60"
  peptide_table = {
    "DNC":  gaussian([0.07975,1.21804]),
    "DNCO": gaussian([0.15475,1.48297]),
    "DO2":  gaussian([0.13738,1.26719]),
    "DCC":  gaussian([0.20249,1.75293]),
    "DCH3": gaussian([0.25038,1.54089]),
    "DNH":  gaussian([0.17120,1.24466])
                  }
  phe_ring_table = {
    "DCGB": gaussian([0.17071,1.63858]),
    "DCDG": gaussian([0.30051,2.06493]),
    "DCED": gaussian([0.27138,1.89454]),
    "DCZE": gaussian([0.26426,1.84996]),
    "DCH":  gaussian([0.26459,1.55955])
                   }
