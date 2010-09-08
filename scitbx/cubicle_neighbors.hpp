#ifndef SCITBX_CUBICLE_NEIGHBORS_HPP
#define SCITBX_CUBICLE_NEIGHBORS_HPP

#include <scitbx/cubicles.h>
#include <vector>

namespace scitbx {

  template <typename FloatType=double>
  struct cubicle_neighbors
  {
    struct cubicle_site
    {
      std::size_t i_seq;
      vec3<FloatType> site_cart;

      cubicle_site(
        std::size_t const& i_seq_,
        vec3<FloatType> const& site_cart_)
      :
        i_seq(i_seq_),
        site_cart(site_cart_)
      {}
    };

    protected:
      typedef std::vector<cubicle_site> cubicle_content_t;
      cubicles<cubicle_content_t, FloatType> cubicles_;

      public:

    cubicle_neighbors() {}

    cubicle_neighbors(
      af::const_ref<vec3<FloatType> > const& main_sites_cart,
      FloatType const& cubicle_edge,
      FloatType const& epsilon)
    :
      cubicles_(main_sites_cart, cubicle_edge, epsilon)
    {
      for(std::size_t i_seq=0;i_seq<main_sites_cart.size();i_seq++) {
        std::size_t i1d_cub = cubicles_.ref.accessor()(
          cubicles_.i_cubicle(main_sites_cart[i_seq]));
        cubicles_.ref[i1d_cub].push_back(
          cubicle_site(i_seq, main_sites_cart[i_seq]));
      }
    }

    std::map<int, std::vector<unsigned> >
    neighbors_of(
      af::const_ref<vec3<FloatType> > const& other_sites_cart,
      FloatType const& distance_cutoff_sq)
    {
      SCITBX_ASSERT(distance_cutoff_sq < fn::pow2(cubicles_.cubicle_edge));
      std::map<int, std::vector<unsigned> > result;
      af::c_grid<3, unsigned> const& n_cubicles = cubicles_.ref.accessor();
      for(std::size_t j_seq=0;j_seq<other_sites_cart.size();j_seq++) {
        vec3<int> j_cub = cubicles_.j_cubicle(other_sites_cart[j_seq]);
        vec3<int> i_cub;
        for(i_cub[0]=j_cub[0]-1;i_cub[0]<=j_cub[0]+1;i_cub[0]++) {
          if (i_cub[0] < 0 || i_cub[0] >= n_cubicles[0]) continue;
        for(i_cub[1]=j_cub[1]-1;i_cub[1]<=j_cub[1]+1;i_cub[1]++) {
          if (i_cub[1] < 0 || i_cub[1] >= n_cubicles[1]) continue;
        for(i_cub[2]=j_cub[2]-1;i_cub[2]<=j_cub[2]+1;i_cub[2]++) {
          if (i_cub[2] < 0 || i_cub[2] >= n_cubicles[2]) continue;
          cubicle_content_t cub = cubicles_.ref(i_cub);
          for(unsigned i=0;i<cub.size();i++) {
            FloatType distance_sq = (
              cub[i].site_cart - other_sites_cart[j_seq]).length_sq();
            if (distance_sq <= distance_cutoff_sq) {
              result[j_seq].push_back(cub[i].i_seq);
            }
          }
        }}}
      }
      return result;
    }
  };

} // namespace scitbx

#endif // GUARD
