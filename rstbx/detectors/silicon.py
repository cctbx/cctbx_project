from __future__ import absolute_import, division, print_function
from six.moves import range
crystalline_density = 2.3290 / 1000.# g/(mm^3) #wikipedia

# Custer JS et al (1994) Density of amorphous Si. Applied Phys. Lett. 64: 437-439.
amorphous_density = crystalline_density / 1.018

class Si_mass_attenuation:
  """Parameterized fit provided by Miroslav Kobas, Dectris Ltd.
  valid range, 2000 - 50000 eV
  Source data for parameter fit were obtained at http://www-cxro.lbl.gov/optical_constants/pert_form.html"""

  # XXX problem with this fit?  I don't get the same fit-to-tabulated that Miro got. 0.1% error
  # XXX Later; try to refit the tabulated data
  def from_energy_eV(self,eV):
    a = 197.625
    b = -1.66146E7
    c = 5.64362E11
    d = -7.64354E15
    e = 1.36449E20
    f = -1.15528E24
    g = 6.54438E27
    h = -2.53016E31
    i = 6.54905E34
    j = -1.08207E38
    k = 1.02849E41
    l = -4.26636E43
    """Returns photoabsorption cross section in (mm^2)/g"""
    if not 2000 <= eV <= 50000: raise Exception("Input energy %7.1feV out of range"%eV)
    xx = 1./eV
    return a + xx*(b + xx*(c + xx*(d + xx*(e + xx*(f + xx*(g + xx*(h + xx*(i + xx*(j + xx*(k + xx*l))))))))))
    #return a + b*xx + c*xx**2 + d*xx**3 + e*xx**4 + f*xx**5 + g*xx**6 + h*xx**7 + i*xx**8 + j*xx**9 + k*xx**10 + l*xx**11

  def from_wavelength_Angstrom(self,Ang):
    return self.from_energy_eV(eV_from_Ang(Ang))

def eV_from_Ang(Ang):
  return 12398.424/Ang

if __name__=="__main__":
  for eV in range(2000, 30500, 500):
    print(eV, "%.2f"%Si_mass_attenuation().from_energy_eV(eV))
  print(Si_mass_attenuation().from_wavelength_Angstrom(1.3))
  print(Si_mass_attenuation().from_wavelength_Angstrom(1.297461)*amorphous_density)
