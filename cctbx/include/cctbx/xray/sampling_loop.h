#if !defined(CCTBX_XRAY_SAMPLED_MODEL_DENSITY_H) \
 && !defined(CCTBX_XRAY_FAST_GRADIENTS_H)
#error "Do not include this file directly."
#endif
#ifndef DOXYGEN_SHOULD_SKIP_THIS
      // highly hand-optimized loop over points in shell
      grid_point_type g_min = pivot - shell.radii;
      grid_point_type g_max = pivot + shell.radii;
      grid_point_type gp;
      grid_point_element_type g01, g0112;
      FloatType f0, f1, f2;
      FloatType c00, c01, c11;
      for(gp[0] = g_min[0]; gp[0] <= g_max[0]; gp[0]++) {
        g01 = math::mod_positive(gp[0],grid_f[0]) * grid_a[1];
        f0 = FloatType(gp[0]) / grid_f[0] - coor_frac[0];
        c00 = orth_mx[0] * f0;
      for(gp[1] = g_min[1]; gp[1] <= g_max[1]; gp[1]++) {
        g0112 = (g01+math::mod_positive(gp[1],grid_f[1])) * grid_a[2];
        f1 = FloatType(gp[1]) / grid_f[1] - coor_frac[1];
        c01 = orth_mx[1] * f1 + c00;
        c11 = orth_mx[4] * f1;
      for(gp[2] = g_min[2]; gp[2] <= g_max[2]; gp[2]++) {
        f2 = FloatType(gp[2]) / grid_f[2] - coor_frac[2];
        scitbx::vec3<FloatType> d(
          orth_mx[2] * f2 + c01,
          orth_mx[5] * f2 + c11,
          orth_mx[8] * f2);
        FloatType d_sq = d.length_sq();
        if (d_sq > shell.max_d_sq) continue;
        std::size_t i_map = g0112 + math::mod_positive(gp[2],grid_f[2]);
#endif // DOXYGEN_SHOULD_SKIP_THIS
