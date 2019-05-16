from __future__ import absolute_import, division, print_function
from dials.array_family import flex
import math
from six.moves import zip

def residual_map_special_deltapsi_add_on( reflections,experiments,matches,hkllist, predicted,plot,eta_deg,deff ):

        detector = experiments[0].detector
        crystal = experiments[0].crystal
        unit_cell = crystal.get_unit_cell()
        pxlsz = detector[0].get_pixel_size()
        model_millers = reflections["miller_index"]
        dpsi = flex.double()
        for match in matches:

          obs_miller = hkllist[match["pred"]]
          model_index= model_millers.first_index(obs_miller)

          raw_delta_psi = reflections["delpsical.rad"][model_index]
          deltapsi_envelope = (unit_cell.d(obs_miller)/deff) + math.pi*eta_deg/180.
          normalized_delta_psi = raw_delta_psi/deltapsi_envelope

          dpsi.append( normalized_delta_psi )

        from matplotlib.colors import Normalize
        dnorm = Normalize()
        dnorm.autoscale(dpsi.as_numpy_array())

        CMAP = plot.get_cmap("bwr")
        for match,dcolor in zip(matches,dpsi):

          #print dcolor, dnorm(dcolor),  CMAP(dnorm(dcolor))
          #blue represents negative delta psi:  outside Ewald sphere; red, positive, inside Ewald sphere
          plot.plot([predicted[match["pred"]][1]/pxlsz[1]],[-predicted[match["pred"]][0]/pxlsz[0]],color=CMAP(dnorm(dcolor)),
          marker=".", markersize=5)
