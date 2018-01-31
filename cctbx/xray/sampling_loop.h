#if !defined(CCTBX_XRAY_SAMPLED_MODEL_DENSITY_H) \
 && !defined(CCTBX_XRAY_FAST_GRADIENTS_H)
#error "Do not include this file directly."
#endif
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define CCTBX_XRAY_SAMPLING_LOOP_END }}}
      // highly hand-optimized loop over points in sampling box
      scitbx::vec3<int> g_min = sampling_box.box_min;
      scitbx::vec3<int> g_max = sampling_box.box_max;
      scitbx::vec3<int> gp;
      gp1g.clear();
      o1f1_.clear();
      o4f1_.clear();
      gp2g.clear();
      o2f2_.clear();
      o5f2_.clear();
      o8f2_.clear();
      int n_gp1 = sampling_box.box_edges[1];
      int n_gp2 = sampling_box.box_edges[2];
      gp1g.reserve(n_gp1);
      o1f1_.reserve(n_gp1);
      o4f1_.reserve(n_gp1);
      gp2g.reserve(n_gp2);
      o2f2_.reserve(n_gp2);
      o5f2_.reserve(n_gp2);
      o8f2_.reserve(n_gp2);
      for(gp[1] = g_min[1]; gp[1] <= g_max[1]; gp[1]++) {
        gp1g.push_back(scitbx::math::mod_positive(gp[1],grid_f[1]));
      }
      for(gp[1] = g_min[1]; gp[1] <= g_max[1]; gp[1]++) {
        FloatType f1 = FloatType(gp[1]) / grid_f[1] - coor_frac[1];
        o1f1_.push_back(orth_mx[1]*f1);
        o4f1_.push_back(orth_mx[4]*f1);
      }
      for(gp[2] = g_min[2]; gp[2] <= g_max[2]; gp[2]++) {
        gp2g.push_back(scitbx::math::mod_positive(gp[2],grid_f[2]));
      }
      for(gp[2] = g_min[2]; gp[2] <= g_max[2]; gp[2]++) {
        FloatType f2 = FloatType(gp[2]) / grid_f[2] - coor_frac[2];
        o2f2_.push_back(orth_mx[2]*f2);
        o5f2_.push_back(orth_mx[5]*f2);
        o8f2_.push_back(orth_mx[8]*f2);
      }
      int g_min0 = g_min[0];
      int g_max0 = g_max[0];
      for(int gp0=g_min0;gp0<=g_max0;gp0++) {
        int g01 = scitbx::math::mod_positive(gp0, grid_f[0]) * grid_a[1];
        FloatType f0 = FloatType(gp0) / grid_f[0] - coor_frac[0];
        FloatType c00 = orth_mx[0] * f0;
      for(int i_gp1=0;i_gp1<n_gp1;i_gp1++) {
        /* It is possible that (g01+gp1g[i_gp1]) and grid_a[2] are large enough
           so that their product is greater than INT_MAX, resulting in g0112
           negative. Later asigning it to std::size_t this leads to a disaster.
        int g0112 = (g01+gp1g[i_gp1]) * grid_a[2];
        */
        std::size_t g0112 = (std::size_t(g01)+std::size_t(gp1g[i_gp1])) *
                                                         std::size_t(grid_a[2]);
        FloatType c01 = o1f1_[i_gp1] + c00;
        FloatType c11 = o4f1_[i_gp1];
      for(int i_gp2=0;i_gp2<n_gp2;i_gp2++) {
        scitbx::vec3<FloatType> d(
          o2f2_[i_gp2] + c01,
          o5f2_[i_gp2] + c11,
          o8f2_[i_gp2]);
        FloatType d_sq = d.length_sq();
        if (d_sq > sampling_box.max_d_sq) continue;
        /* See comment above: g0112 can become negative.
        std::size_t i_map = g0112 + gp2g[i_gp2];
        */
        std::size_t i_map = g0112 + std::size_t(gp2g[i_gp2]);
#endif // DOXYGEN_SHOULD_SKIP_THIS
