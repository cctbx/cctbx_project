#ifndef CCTBX_SGTBX_DIRECT_SPACE_ASU_H
#define CCTBX_SGTBX_DIRECT_SPACE_ASU_H

#include <cctbx/sgtbx/sym_equiv_sites.h>
#include <scitbx/math/minimum_covering_sphere.h>

namespace cctbx { namespace sgtbx {

//! Direct-space asymmetric units.
namespace direct_space_asu {

  //! Floating-point parameterization of a cut plane.
  template <typename FloatType=double>
  class float_cut_plane
  {
    public:
      //! Default constructor. Some data members are not initialized!
      float_cut_plane() {}

      //! Initialization with normal vector n and constant c.
      float_cut_plane(
        fractional<FloatType> const& n_,
        FloatType const& c_)
      :
        n(n_),
        c(c_)
      {}

      //! Normal vector.
      fractional<FloatType> n;

      //! Constant.
      FloatType c;

      //! Returns n * point + c.
      FloatType
      evaluate(fractional<FloatType> const& point) const
      {
        return n * point + c;
      }

      //! Equivalent to evaluate(point) >= -epsilon.
      bool
      is_inside(
        fractional<FloatType> const& point,
        FloatType const& epsilon=0) const
      {
        if (evaluate(point) < -epsilon) return false;
        return true;
      }

      //! Returns a point in the cut plane.
      fractional<FloatType>
      get_point_in_plane() const
      {
        fractional<FloatType> result(0,0,0);
        for(std::size_t i=0;i<3;i++) {
          if (n[i] != 0) {
            result[i] = -c / n[i];
            return result;
          }
        }
        throw error("float_cut_plane normal vector is the null vector.");
      }

      //! Shifts the plane by the distance specified as thickness.
      float_cut_plane
      add_buffer(
        uctbx::unit_cell const& unit_cell,
        FloatType const& thickness) const
      {
        typedef fractional<FloatType> f_t;
        typedef cartesian<FloatType> c_t;
        f_t x_frac = get_point_in_plane();
        c_t x_cart = unit_cell.orthogonalize(x_frac);
        c_t n_cart = (n * unit_cell.fractionalization_matrix()).normalize();
        c_t y_cart = x_cart - n_cart * thickness;
        f_t y_frac = unit_cell.fractionalize(y_cart);
        return float_cut_plane(n, -n * y_frac);
      }
  };

  //! Floating-point parameterization of an asymmetric unit.
  template <typename FloatType=double>
  class float_asu
  {
    public:
      //! Array type for facets.
      typedef af::small<float_cut_plane<FloatType>, 9> facets_t;

      //! Default constructor. Some data members are not initialized!
      float_asu() {}

      //! Initialization with unit cell and list of facets.
      float_asu(
        uctbx::unit_cell const& unit_cell,
        facets_t const& facets,
        FloatType const& is_inside_epsilon=1.e-6)
      :
        unit_cell_(unit_cell),
        facets_(facets),
        is_inside_epsilon_(is_inside_epsilon),
        have_box_(false)
      {}

      //! Unit cell as passed to the constructor.
      uctbx::unit_cell const&
      unit_cell() const { return unit_cell_; }

      //! Facets as passed to the constructor.
      facets_t const&
      facets() const { return facets_; }

      //! Epsilon for is_inside() as passed to the constructor.
      FloatType
      is_inside_epsilon() const { return is_inside_epsilon_; }

      //! True if is_inside() is true for all facets().
      /*! is_inside_epsilon() is used in the test.
       */
      bool
      is_inside(fractional<FloatType> const& point) const
      {
        for(std::size_t i=0;i<facets_.size();i++) {
          if (!facets_[i].is_inside(point, is_inside_epsilon_)) return false;
        }
        return true;
      }

      /*! \brief New asymmetric unit with all facets shifted by the distance
          specified as thickness.
       */
      float_asu
      add_buffer(FloatType const& thickness) const
      {
        facets_t buffer_facets;
        for(std::size_t i=0;i<facets_.size();i++) {
          buffer_facets.push_back(
            facets_[i].add_buffer(unit_cell_, thickness));
        }
        return float_asu(unit_cell_, buffer_facets, is_inside_epsilon_);
      }

      //! List of vertices. Duplicates are possible.
      /*! epsilon is used to detect singular matrices and must be
          strictly greater than zero.
       */
      af::shared<scitbx::vec3<FloatType> >
      volume_vertices(
        bool cartesian=false,
        FloatType const& epsilon=1.e-6) const
      {
        CCTBX_ASSERT(epsilon > 0);
        af::shared<scitbx::vec3<FloatType> > result;
        scitbx::mat3<FloatType> m;
        scitbx::vec3<FloatType> b;
        for(std::size_t i0=0   ;i0<facets_.size()-2;i0++) {
          m.set_row(0, facets_[i0].n);
          b[0] = -facets_[i0].c;
        for(std::size_t i1=i0+1;i1<facets_.size()-1;i1++) {
          m.set_row(1, facets_[i1].n);
          b[1] = -facets_[i1].c;
        for(std::size_t i2=i1+1;i2<facets_.size()  ;i2++) {
          m.set_row(2, facets_[i2].n);
          b[2] = -facets_[i2].c;
          FloatType d = m.determinant();
          scitbx::mat3<FloatType> c = m.co_factor_matrix_transposed();
          if (scitbx::fn::absolute(d) > c.max_abs() * epsilon) {
            c /= d;
            fractional<FloatType> vertex = c * b;
            if (is_inside(vertex)) {
              if (cartesian) {
                result.push_back(unit_cell_.orthogonalize(vertex));
              }
              else {
                result.push_back(vertex);
              }
            }
          }
        }}}
        return result;
      }

      /*! \brief Fractional coordinates of lower-left corner of
          box covering the asymmetric unit.
       */
      fractional<FloatType> const&
      box_min() const
      {
        if (!have_box_) compute_box();
        return box_min_;
      }

      /*! \brief Fractional coordinates of upper-right corner of
          box covering the asymmetric unit.
       */
      fractional<FloatType> const&
      box_max() const
      {
        if (!have_box_) compute_box();
        return box_max_;
      }

    protected:
      uctbx::unit_cell unit_cell_;
      facets_t facets_;
      FloatType is_inside_epsilon_;
      mutable bool have_box_;
      mutable fractional<FloatType> box_min_;
      mutable fractional<FloatType> box_max_;

      void
      compute_box() const
      {
        af::shared<scitbx::vec3<FloatType> > vertices_ = volume_vertices();
        af::const_ref<scitbx::vec3<FloatType> > vertices = vertices_.ref();
        CCTBX_ASSERT(vertices.size() >= 4);
        box_min_ = box_max_ = vertices[0];
        for(std::size_t i=1;i<vertices.size();i++) {
          scitbx::vec3<FloatType> const& vertex = vertices[i];
          for(std::size_t j=0;j<3;j++) {
            scitbx::math::update_min(box_min_[j], vertex[j]);
            scitbx::math::update_max(box_max_[j], vertex[j]);
          }
        }
        have_box_ = true;
      }
  };

  //! Mapping of a site to an asymmetric unit.
  template <typename FloatType=double, typename IntShiftType=int>
  class asu_mapping
  {
    public:
      //! Default constructor. Some data members are not initialized!
      asu_mapping() {}

      //! Grouping of parameters that define the mapping.
      /*! Not available in Python.
       */
      asu_mapping(
        std::size_t i_sym_op,
        scitbx::vec3<IntShiftType> const& unit_shifts,
        fractional<FloatType> const& mapped_site)
      :
        i_sym_op_(i_sym_op),
        unit_shifts_(unit_shifts),
        mapped_site_(mapped_site)
      {}

      //! Index of symmetry operation.
      std::size_t
      i_sym_op() const { return i_sym_op_; }

      //! Additional unit shifts.
      scitbx::vec3<IntShiftType> const&
      unit_shifts() const { return unit_shifts_; }

      //! Fractional coordinates of the mapped site.
      /*! mapped_site = space_group(i_sym_op()) * orginal_site + unit_shifts()
       */
      fractional<FloatType> const&
      mapped_site() const { return mapped_site_; }

    private:
      std::size_t i_sym_op_;
      scitbx::vec3<IntShiftType> unit_shifts_;
      fractional<FloatType> mapped_site_;
  };

  //! Container for mapping of sites to an asymmetric unit.
  /*! The array of mappings is kept together with the parameters
      used to compute them to maximize consistency.
   */
  template <typename FloatType=double, typename IntShiftType=int>
  class asu_mappings
  {
    public:
      //! Type of array of mappings for one site.
      typedef std::vector<asu_mapping<FloatType, IntShiftType> >
        array_of_mappings_for_one_site;

      //! Type of array of array of mappings for one site.
      typedef af::shared<array_of_mappings_for_one_site>
        array_of_array_of_mappings_for_one_site;

      //! Default constructor. Some data members are not initialized!
      asu_mappings();

      //! Grouping of the data needed to compute the mappings.
      asu_mappings(
        sgtbx::space_group const& space_group,
        float_asu<FloatType> const& asu,
        FloatType const& buffer_thickness,
        FloatType const& sym_equiv_epsilon=1.e-6)
      :
        space_group_(space_group),
        asu_(asu),
        buffer_thickness_(buffer_thickness),
        sym_equiv_epsilon_(sym_equiv_epsilon),
        asu_buffer_(asu.add_buffer(buffer_thickness)),
        sym_equiv_tolerance_(
          std::pow(asu.unit_cell().volume(), 1/3.) * sym_equiv_epsilon),
        sym_equiv_minimum_distance_(sym_equiv_tolerance_ * 10),
        buffer_covering_sphere_(
          scitbx::math::minimum_covering_sphere_3d<FloatType>(
            asu.volume_vertices(true).const_ref())
          .expand(buffer_thickness+sym_equiv_tolerance_))
      {}

      //! Pre-allocates memory for mappings(); for efficiency.
      void
      reserve(std::size_t n_sites_final) { mappings_.reserve(n_sites_final); }

      //! Space group as passed to the constructor.
      sgtbx::space_group const&
      space_group() const { return space_group_; }

      //! Asymmetric unit as passed to the constructor.
      float_asu<FloatType> const&
      asu() const { return asu_; }

      //! Equivalent to asu().unit_cell().
      uctbx::unit_cell const&
      unit_cell() const { return asu_.unit_cell(); }

      //! Buffer thickness as passed to the constructor.
      FloatType
      buffer_thickness() const { return buffer_thickness_; }

      /*! \brief Reference to internalized result of
          asu().add_buffer(buffer_thickness()).
       */
      float_asu<FloatType> const&
      asu_buffer() const { return asu_buffer_; }

      //! Value of sym_equiv_epsilon as passed to the constructor.
      FloatType
      sym_equiv_epsilon() const { return sym_equiv_epsilon_; }

      //! Sphere covering the asymmetric unit + buffer_thickness().
      /*! The sphere is computed as the minimum covering sphere
          for the vertices of asu(), followed by adding buffer_thickness()
          to the radius. In general some of the vertices of buffer_asu()
          may be outside buffer_covering_sphere().
       */
      scitbx::math::sphere_3d<FloatType> const&
      buffer_covering_sphere() const { return buffer_covering_sphere_; }

      //! Processes one site and appends the results to mappings().
      fractional<FloatType>
      process(fractional<FloatType> const& original_site)
      {
        mappings_.push_back(array_of_mappings_for_one_site());
        array_of_mappings_for_one_site& site_mappings = mappings_.back();
        sgtbx::sym_equiv_sites<FloatType> equiv_sites(
          asu_.unit_cell(),
          space_group_,
          original_site,
          sym_equiv_minimum_distance_,
          sym_equiv_tolerance_);
        af::const_ref<typename sgtbx::sym_equiv_sites<FloatType>::coor_t>
          coordinates = equiv_sites.coordinates().const_ref();
        af::const_ref<std::size_t>
          sym_op_indices = equiv_sites.sym_op_indices().const_ref();
        bool have_site_in_asu = false;
        for(std::size_t i_sym_eq=0;i_sym_eq<coordinates.size();i_sym_eq++) {
          scitbx::vec3<FloatType> const& site = coordinates[i_sym_eq];
          scitbx::vec3<IntShiftType> unit_shifts_min;
          scitbx::vec3<IntShiftType> unit_shifts_max;
          for(std::size_t i=0;i<3;i++) {
            unit_shifts_min[i] = scitbx::math::iceil(
              asu_buffer_.box_min()[i] - site[i]);
            unit_shifts_max[i] = scitbx::math::ifloor(
              asu_buffer_.box_max()[i] - site[i]);
          }
          scitbx::vec3<IntShiftType> u;
          fractional<FloatType> mapped_site;
          for(u[0]=unit_shifts_min[0];u[0]<=unit_shifts_max[0];u[0]++) {
            mapped_site[0] = site[0] + u[0];
          for(u[1]=unit_shifts_min[1];u[1]<=unit_shifts_max[1];u[1]++) {
            mapped_site[1] = site[1] + u[1];
          for(u[2]=unit_shifts_min[2];u[2]<=unit_shifts_max[2];u[2]++) {
            mapped_site[2] = site[2] + u[2];
            if (   asu_buffer_.is_inside(mapped_site)
                && buffer_covering_sphere_.is_inside(
                     asu_.unit_cell().orthogonalize(mapped_site))) {
              asu_mapping<FloatType, IntShiftType>
                mapping(sym_op_indices[i_sym_eq], u, mapped_site);
              if (!have_site_in_asu && asu_.is_inside(mapped_site)) {
                site_mappings.insert(site_mappings.begin(), mapping);
                have_site_in_asu = true;
              }
              else {
                site_mappings.push_back(mapping);
              }
            }
          }}}
        }
        CCTBX_ASSERT(have_site_in_asu);
        return site_mappings[0].mapped_site();
      }

      //! Accumulated mappings due to repeated calls of process().
      array_of_array_of_mappings_for_one_site const&
      mappings() const { return mappings_; }

    protected:
      sgtbx::space_group space_group_;
      float_asu<FloatType> asu_;
      FloatType buffer_thickness_;
      FloatType sym_equiv_epsilon_;
      float_asu<FloatType> asu_buffer_;
      FloatType sym_equiv_tolerance_;
      FloatType sym_equiv_minimum_distance_;
      scitbx::math::sphere_3d<FloatType> buffer_covering_sphere_;
      array_of_array_of_mappings_for_one_site mappings_;
  };

}}} // namespace cctbx::sgtbx::direct_space_asu

#endif // CCTBX_SGTBX_DIRECT_SPACE_ASU_H
