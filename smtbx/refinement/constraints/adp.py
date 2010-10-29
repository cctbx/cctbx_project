import iotbx.constraints.adp as _input
import smtbx.refinement.constraints as _

class u_iso_proportional_to_pivot_u_eq(_input.u_iso_proportional_to_pivot_u_eq):

  def add_to(self, reparametrisation):
    param = reparametrisation.add(
      _.u_iso_proportional_to_pivot_u_eq,
      pivot_u=reparametrisation.add_new_thermal_displacement_parameter(
        self.u_eq_scatterer_idx),
      scatterer = reparametrisation.structure.scatterers()[
        self.u_iso_scatterer_idx],
      multiplier=self.multiplier)
    reparametrisation.asu_scatterer_parameters[
      self.u_iso_scatterer_idx].u = param
