from eltbx.caasf_wk1995 import CAASF_WK1995
import sftbx
sf = CAASF_WK1995("Si4+")
s = sftbx.XrayScatterer("Si", sf, (0.1, 0.2, 0.3), 1., 0.03, 0j)
print s
