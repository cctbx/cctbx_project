#ifndef CCTBX_CRYSTAL_DIRECT_SPACE_ASU_H
#define CCTBX_CRYSTAL_DIRECT_SPACE_ASU_H

#include <cctbx/sgtbx/sym_equiv_sites.h>
#include <cctbx/sgtbx/site_symmetry_table.h>
#include <scitbx/math/minimum_covering_sphere.h>

namespace cctbx { namespace crystal {

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
        c_t n_cart = unit_cell.v_times_fractionalization_matrix_transpose(
          /* v */ n).normalize();
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
      //! Array type for cuts.
      typedef af::small<float_cut_plane<FloatType>, 12> cuts_t;

      //! Default constructor. Some data members are not initialized!
      float_asu() {}

      //! Initialization with unit cell and list of cuts.
      float_asu(
        uctbx::unit_cell const& unit_cell,
        cuts_t const& cuts,
        FloatType const& is_inside_epsilon=1e-6)
      :
        unit_cell_(unit_cell),
        cuts_(cuts),
        is_inside_epsilon_(is_inside_epsilon),
        have_box_(false)
      {}

      //! Unit cell as passed to the constructor.
      uctbx::unit_cell const&
      unit_cell() const { return unit_cell_; }

      //! cuts as passed to the constructor.
      cuts_t const&
      cuts() const { return cuts_; }

      //! Epsilon for is_inside() as passed to the constructor.
      FloatType
      is_inside_epsilon() const { return is_inside_epsilon_; }

      //! True if is_inside() is true for all cuts().
      /*! is_inside_epsilon() is used in the test.
       */
      bool
      is_inside(fractional<FloatType> const& point) const
      {
        for(std::size_t i=0;i<cuts_.size();i++) {
          if (!cuts_[i].is_inside(point, is_inside_epsilon_)) return false;
        }
        return true;
      }

      //! Array version of is_inside(), given fractional coordinates.
      af::shared<bool>
      is_inside_frac(
        af::const_ref<scitbx::vec3<FloatType> > const& sites_frac)
      {
        af::shared<bool>
          result(sites_frac.size(), af::init_functor_null<bool>());
        bool* res = result.begin();
        for(std::size_t i=0;i<sites_frac.size();i++) {
          *res++ = is_inside(sites_frac[i]);
        }
        return result;
      }

      //! Array version of is_inside(), given cartesian coordinates.
      af::shared<bool>
      is_inside_cart(
        af::const_ref<scitbx::vec3<FloatType> > const& sites_cart)
      {
        af::shared<bool>
          result(sites_cart.size(), af::init_functor_null<bool>());
        uctbx::unit_cell const& uc = unit_cell();
        bool* res = result.begin();
        for(std::size_t i=0;i<sites_cart.size();i++) {
          *res++ = is_inside(uc.fractionalize(sites_cart[i]));
        }
        return result;
      }

      /*! \brief New asymmetric unit with all cuts shifted by the distance
          specified as thickness.
       */
      float_asu
      add_buffer(FloatType const& thickness) const
      {
        cuts_t buffer_cuts;
        for(std::size_t i=0;i<cuts_.size();i++) {
          buffer_cuts.push_back(
            cuts_[i].add_buffer(unit_cell_, thickness));
        }
        return float_asu(unit_cell_, buffer_cuts, is_inside_epsilon_);
      }

      //! List of vertices. Duplicates are possible.
      /*! epsilon is used to detect singular matrices and must be
          strictly greater than zero.
       */
      af::shared<scitbx::vec3<FloatType> >
      shape_vertices(
        bool cartesian=false,
        FloatType const& epsilon=1e-6) const
      {
        CCTBX_ASSERT(epsilon > 0);
        af::shared<scitbx::vec3<FloatType> > result;
        if (cuts_.size() < 3) return result;
        scitbx::mat3<FloatType> m;
        scitbx::vec3<FloatType> b;
        for(std::size_t i0=0   ;i0<cuts_.size()-2;i0++) {
          m.set_row(0, cuts_[i0].n);
          b[0] = -cuts_[i0].c;
        for(std::size_t i1=i0+1;i1<cuts_.size()-1;i1++) {
          m.set_row(1, cuts_[i1].n);
          b[1] = -cuts_[i1].c;
        for(std::size_t i2=i1+1;i2<cuts_.size()  ;i2++) {
          m.set_row(2, cuts_[i2].n);
          b[2] = -cuts_[i2].c;
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

      /*! \brief Coordinates of lower-left corner of
          box covering the asymmetric unit.
       */
      scitbx::vec3<FloatType> const&
      box_min(bool cartesian=false) const
      {
        if (!have_box_) compute_box();
        if (cartesian) return box_min_cart_;
        return box_min_frac_;
      }

      /*! \brief Coordinates of upper-right corner of
          box covering the asymmetric unit.
       */
      scitbx::vec3<FloatType> const&
      box_max(bool cartesian=false) const
      {
        if (!have_box_) compute_box();
        if (cartesian) return box_max_cart_;
        return box_max_frac_;
      }

      /*! box_max() - box_min(). Not available in Python.
       */
      scitbx::vec3<FloatType>
      box_span(bool cartesian=false) const
      {
        if (!have_box_) compute_box();
        if (cartesian) return box_max_cart_ - box_min_cart_;
        return box_max_frac_ - box_min_frac_;
      }

    protected:
      uctbx::unit_cell unit_cell_;
      cuts_t cuts_;
      FloatType is_inside_epsilon_;
      mutable bool have_box_;
      mutable scitbx::vec3<FloatType> box_min_frac_;
      mutable scitbx::vec3<FloatType> box_max_frac_;
      mutable scitbx::vec3<FloatType> box_min_cart_;
      mutable scitbx::vec3<FloatType> box_max_cart_;

      void
      compute_box() const
      {
        af::shared<scitbx::vec3<FloatType> > vertices_ = shape_vertices();
        af::const_ref<scitbx::vec3<FloatType> > vertices = vertices_.ref();
        CCTBX_ASSERT(vertices.size() >= 4);
        box_min_frac_ = box_max_frac_ = vertices[0];
        box_min_cart_ = box_max_cart_ = unit_cell_.orthogonalize(vertices[0]);
        for(std::size_t i=1;i<vertices.size();i++) {
          scitbx::vec3<FloatType> vertex = vertices[i];
          for(std::size_t j=0;j<3;j++) {
            scitbx::math::update_min(box_min_frac_[j], vertex[j]);
            scitbx::math::update_max(box_max_frac_[j], vertex[j]);
          }
          vertex = unit_cell_.orthogonalize(vertex);
          for(std::size_t j=0;j<3;j++) {
            scitbx::math::update_min(box_min_cart_[j], vertex[j]);
            scitbx::math::update_max(box_max_cart_[j], vertex[j]);
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
        unsigned i_sym_op,
        scitbx::vec3<IntShiftType> const& unit_shifts,
        cartesian<FloatType> const& mapped_site)
      :
        i_sym_op_(i_sym_op),
        unit_shifts_(unit_shifts),
        mapped_site_(mapped_site)
      {}

      //! Index of symmetry operation.
      unsigned
      i_sym_op() const { return i_sym_op_; }

      //! Additional unit shifts.
      scitbx::vec3<IntShiftType> const&
      unit_shifts() const { return unit_shifts_; }

      //! Cartesian coordinates of the mapped site.
      /*! mapped_site = space_group(i_sym_op()) * orginal_site + unit_shifts()
       */
      cartesian<FloatType> const&
      mapped_site() const { return mapped_site_; }

    private:
      unsigned i_sym_op_;
      scitbx::vec3<IntShiftType> unit_shifts_;
      cartesian<FloatType> mapped_site_;
  };

  //! Grouping of indices for site in asu_mappings container.
  /*! Not available in Python.
   */
  struct asu_mapping_index
  {
    asu_mapping_index() {}

    asu_mapping_index(unsigned i_seq_, unsigned i_sym_)
    :
      i_seq(i_seq_),
      i_sym(i_sym_)
    {}

    bool
    operator<(asu_mapping_index const& other) const
    {
      if (i_seq < other.i_seq) return true;
      if (i_seq > other.i_seq) return false;
      if (i_sym < other.i_sym) return true;
      return false;
    }

    unsigned i_seq;
    unsigned i_sym;
  };

  //! Grouping of indices for pair of sites in asu_mappings container.
  struct asu_mapping_index_pair
  {
    //! Main index of first site.
    unsigned i_seq;
    //! Main index of second site.
    unsigned j_seq;
    //! Symmetry index of second site.
    unsigned j_sym;

    /*! If minimal is false:
          Inside-inside upper triangle (i_seq < j_seq) or
          inside-outside (j_sym != 0)
        If minimal is true:
          i_seq <= j_seq
     */
    bool
    is_active(bool minimal=false) const
    {
      return i_seq <= j_seq || (!minimal && j_sym != 0);
    }
  };

  //! asu_mapping_index_pair plus difference vector and distance squared.
  template <typename FloatType=double>
  struct asu_mapping_index_pair_and_diff : asu_mapping_index_pair
  {
    //! Difference vector.
    cartesian<FloatType> diff_vec;
    //! Distance squared.
    FloatType dist_sq;
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
        FloatType const& buffer_thickness)
      :
        space_group_(space_group),
        asu_(asu),
        buffer_thickness_(buffer_thickness),
        space_group_ops_(space_group.all_ops()),
        space_group_ops_const_ref_(space_group_ops_.const_ref()),
        asu_buffer_(asu.add_buffer(buffer_thickness)),
        buffer_covering_sphere_(
          scitbx::math::minimum_covering_sphere_3d<FloatType>(
            asu.shape_vertices(true).const_ref())
          .expand(buffer_thickness)
          .expand_relative(2*asu.is_inside_epsilon())),
        mappings_const_ref_(mappings_.const_ref()),
        n_sites_in_asu_and_buffer_(0),
        mapped_sites_min_(0,0,0),
        mapped_sites_max_(0,0,0)
      {}

      //! Pre-allocates memory for mappings(); for efficiency.
      void
      reserve(std::size_t n_sites_final)
      {
        site_symmetry_table_.reserve(n_sites_final);
        mappings_.reserve(n_sites_final);
        mappings_const_ref_ = mappings_.const_ref();
      }

      //! Support for site_cluster_analysis.
      /*! Not available in Python.
       */
      void
      discard_last()
      {
        site_symmetry_table_.discard_last();
        mappings_.pop_back();
        mappings_const_ref_ = mappings_.const_ref();
      }

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

      //! Sphere covering the asymmetric unit + buffer_thickness().
      /*! The sphere is computed as the minimum covering sphere
          for the vertices of asu(), followed by adding buffer_thickness()
          to the radius. In general some of the vertices of buffer_asu()
          may be outside buffer_covering_sphere().
       */
      scitbx::math::sphere_3d<FloatType> const&
      buffer_covering_sphere() const { return buffer_covering_sphere_; }

      //! Processes one site and appends the results to mappings().
      asu_mappings&
      process(
        fractional<FloatType> const& original_site,
        FloatType const& min_distance_sym_equiv=0.5)
      {
        return process(
          original_site,
          sgtbx::site_symmetry(
            asu_.unit_cell(),
            space_group_,
            original_site,
            min_distance_sym_equiv,
            /*assert_min_distance_sym_equiv*/ true));
      }

      //! Processes one site and appends the results to mappings().
      asu_mappings&
      process(
        fractional<FloatType> const& original_site,
        sgtbx::site_symmetry_ops const& site_symmetry_ops)
      {
        CCTBX_ASSERT(mappings_.begin()
                  == mappings_const_ref_.begin());
        mappings_.push_back(array_of_mappings_for_one_site());
        mappings_const_ref_ = mappings_.const_ref();
        array_of_mappings_for_one_site& site_mappings = mappings_.back();
        site_symmetry_table_.process(site_symmetry_ops);
        sgtbx::sym_equiv_sites<FloatType> equiv_sites(
          asu_.unit_cell(),
          space_group_,
          original_site,
          site_symmetry_ops);
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
              asu_buffer_.box_min()[i] - site[i] - 2*asu_.is_inside_epsilon());
            unit_shifts_max[i] = scitbx::math::ifloor(
              asu_buffer_.box_max()[i] - site[i] + 2*asu_.is_inside_epsilon());
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
              asu_mapping<FloatType, IntShiftType> mapping(
                sym_op_indices[i_sym_eq],
                u,
                asu_.unit_cell().orthogonalize(mapped_site));
              if (!have_site_in_asu && asu_.is_inside(mapped_site)) {
                site_mappings.insert(site_mappings.begin(), mapping);
                have_site_in_asu = true;
              }
              else {
                site_mappings.push_back(mapping);
              }
              n_sites_in_asu_and_buffer_++;
              if (   site_mappings.size() == 1
                  && mappings_const_ref_.size() == 1) {
                mapped_sites_min_ = mapping.mapped_site();
                mapped_sites_max_ = mapping.mapped_site();
              }
              else {
                for(std::size_t i=0;i<3;i++) {
                  FloatType const& e = mapping.mapped_site()[i];
                  scitbx::math::update_min(mapped_sites_min_[i], e);
                  scitbx::math::update_max(mapped_sites_max_[i], e);
                }
              }
            }
          }}}
        }
        CCTBX_ASSERT(have_site_in_asu);
        return *this;
      }

      //! Calls process() for each original site.
      asu_mappings&
      process_sites_frac(
        af::const_ref<scitbx::vec3<FloatType> > const& original_sites,
        FloatType const& min_distance_sym_equiv=0.5)
      {
        for(std::size_t i=0;i<original_sites.size();i++) {
          process(original_sites[i], min_distance_sym_equiv);
        }
        return *this;
      }

      //! Calls process() for each original site.
      asu_mappings&
      process_sites_frac(
        af::const_ref<scitbx::vec3<FloatType> > const& original_sites,
        sgtbx::site_symmetry_table const& site_symmetry_table)
      {
        CCTBX_ASSERT(site_symmetry_table.indices_const_ref().size()
                  == original_sites.size());
        for(std::size_t i=0;i<original_sites.size();i++) {
          process(original_sites[i], site_symmetry_table.get(i));
        }
        return *this;
      }

      //! Calls process() for each original site.
      asu_mappings&
      process_sites_cart(
        af::const_ref<scitbx::vec3<FloatType> > const& original_sites,
        FloatType const& min_distance_sym_equiv=0.5)
      {
        uctbx::unit_cell const& uc = unit_cell();
        for(std::size_t i=0;i<original_sites.size();i++) {
          process(
            uc.fractionalize(original_sites[i]), min_distance_sym_equiv);
        }
        return *this;
      }

      //! Calls process() for each original site.
      asu_mappings&
      process_sites_cart(
        af::const_ref<scitbx::vec3<FloatType> > const& original_sites,
        sgtbx::site_symmetry_table const& site_symmetry_table)
      {
        CCTBX_ASSERT(site_symmetry_table.indices_const_ref().size()
                  == original_sites.size());
        uctbx::unit_cell const& uc = unit_cell();
        for(std::size_t i=0;i<original_sites.size();i++) {
          process(
            uc.fractionalize(original_sites[i]), site_symmetry_table.get(i));
        }
        return *this;
      }

      /*! \brief Total number of sites in the asymmetric unit and the
          surrounding buffer.
       */
      std::size_t
      n_sites_in_asu_and_buffer() const { return n_sites_in_asu_and_buffer_; }

      //! Accumulated mappings due to repeated calls of process().
      array_of_array_of_mappings_for_one_site const&
      mappings() const { return mappings_; }

      //! Use for maximum performance.
      /*! Not available in Python.
       */
      af::const_ref<array_of_mappings_for_one_site> const&
      mappings_const_ref() const { return mappings_const_ref_; }

      //! Minimum coordinates of all mapped sites.
      cartesian<FloatType> const&
      mapped_sites_min() const { return mapped_sites_min_; }

      //! Maximum coordinates of all mapped sites.
      cartesian<FloatType> const&
      mapped_sites_max() const { return mapped_sites_max_; }

      //! mapped_sites_max() - mapped_sites_min().
      cartesian<FloatType>
      mapped_sites_span() const
      {
        return mapped_sites_max_ - mapped_sites_min_;
      }

      //! Special position operator for the i_seq'th processed site.
      sgtbx::rt_mx const&
      special_op(std::size_t i_seq) const
      {
        return site_symmetry_table_.get(i_seq).special_op();
      }

      //! Table of unique site symmetries.
      sgtbx::site_symmetry_table const&
      site_symmetry_table() const { return site_symmetry_table_; }

      //! mappings()[i_seq][i_sym] with range checking of i_seq and i_sym.
      /*! Not available in Python.
       */
      asu_mapping<FloatType, IntShiftType> const&
      get_asu_mapping(std::size_t i_seq, std::size_t i_sym) const
      {
        CCTBX_ASSERT(mappings_const_ref_.begin() == mappings_.begin());
        CCTBX_ASSERT(i_seq < mappings_const_ref_.size());
        CCTBX_ASSERT(i_sym < mappings_const_ref_[i_seq].size());
        return mappings_const_ref_[i_seq][i_sym];
      }

      //! Symmetry operation original_site -> site in asu.
      sgtbx::rt_mx
      get_rt_mx(asu_mapping<FloatType, IntShiftType> const& mapping) const
      {
        sgtbx::rt_mx const&
          rt = space_group_ops_const_ref_[mapping.i_sym_op()];
        int t_den = rt.t().den();
        return rt + sgtbx::tr_vec(mapping.unit_shifts()*t_den, t_den);
      }

      //! Symmetry operation original_site -> site in asu.
      /*! Shorthand for: get_rt_mx(get_asu_mapping(i_seq, i_sym))
       */
      sgtbx::rt_mx
      get_rt_mx(std::size_t i_seq, std::size_t i_sym) const
      {
        return get_rt_mx(get_asu_mapping(i_seq, i_sym));
      }

      //! Symmetry operation original_site -> site in asu.
      /*! Shorthand for: get_rt_mx(pair.i_seq, 0)
       */
      sgtbx::rt_mx
      get_rt_mx_i(asu_mapping_index_pair const& pair) const
      {
        return get_rt_mx(pair.i_seq, 0);
      }

      //! Symmetry operation original_site -> site in asu.
      /*! Shorthand for: get_rt_mx(pair.j_seq, pair.j_sym)
       */
      sgtbx::rt_mx
      get_rt_mx_j(asu_mapping_index_pair const& pair) const
      {
        return get_rt_mx(pair.j_seq, pair.j_sym);
      }

      //! Difference vector for the given pair.
      /*! result = site(j_seq,j_sym) - site(i_seq,0).
       */
      cartesian<FloatType>
      diff_vec(asu_mapping_index_pair const& pair) const
      {
        return get_asu_mapping(pair.j_seq, pair.j_sym).mapped_site()
             - get_asu_mapping(pair.i_seq, 0).mapped_site();
      }

      //! Maps a moved site (e.g. during refinement) to the asymmetric unit.
      cartesian<FloatType>
      map_moved_site_to_asu(
        cartesian<FloatType> const& moved_original_site,
        std::size_t i_seq,
        std::size_t i_sym) const
      {
        if (r_cart_.size() == 0) {
          scitbx::mat3<FloatType> o(unit_cell().orthogonalization_matrix());
          scitbx::mat3<FloatType> f(unit_cell().fractionalization_matrix());
          r_cart_.reserve(space_group_.order_z());
          t_cart_.reserve(space_group_.order_z());
          for(std::size_t i_sym_op=0;
              i_sym_op<space_group_.order_z();
              i_sym_op++) {
            sgtbx::rt_mx const& s = space_group_ops_const_ref_[i_sym_op];
            typedef scitbx::type_holder<FloatType> t_h;
            r_cart_.push_back(o*s.r().as_floating_point(t_h())*f);
            t_cart_.push_back(o*s.t().as_floating_point(t_h()));
          }
        }
        asu_mapping<FloatType, IntShiftType> const&
          am = get_asu_mapping(i_seq, i_sym);
        return r_cart_[am.i_sym_op()] * moved_original_site
             + t_cart_[am.i_sym_op()]
             + scitbx::vec3<FloatType>(
                 unit_cell().orthogonalization_matrix() * am.unit_shifts());
      }

      /*! \brief Rotation part of
          space_group(mappings()[i_seq][i_sym].i_sym_op()).inverse()
          in the cartesian system.
       */
      /*! Useful for mapping difference vectors in cartesian space
          from the asymmetric unit to the original sites.

          The rotation matrices are cached for maximum performance.
       */
      scitbx::mat3<FloatType>
      r_inv_cart(std::size_t i_seq, std::size_t i_sym) const
      {
        if (r_inv_cart_.size() == 0) {
          scitbx::mat3<FloatType> o(unit_cell().orthogonalization_matrix());
          scitbx::mat3<FloatType> f(unit_cell().fractionalization_matrix());
          r_inv_cart_.reserve(space_group_.order_z());
          for(std::size_t i_sym_op=0;
              i_sym_op<space_group_.order_z();
              i_sym_op++) {
            sgtbx::rt_mx const& s = space_group_ops_const_ref_[i_sym_op];
            typedef scitbx::type_holder<FloatType> t_h;
            r_inv_cart_
              .push_back(o*s.r().inverse().as_floating_point(t_h())*f);
          }
        }
        return r_inv_cart_[get_asu_mapping(i_seq, i_sym).i_sym_op()];
      }

      /*! \brief True if the interaction is between sites in general
          positions as passed to process().
       */
      /*! For sites on special positions the return value is always false.
       */
      bool
      is_simple_interaction(asu_mapping_index_pair const& pair) const
      {
        CCTBX_ASSERT(
             pair.i_seq < mappings_const_ref_.size()
          && pair.j_seq < mappings_const_ref_.size()
          && pair.j_sym < mappings_const_ref_[pair.j_seq].size());
        if (   site_symmetry_table_.indices_const_ref()[pair.i_seq]
            || site_symmetry_table_.indices_const_ref()[pair.j_seq]) {
          return false;
        }
        asu_mapping<FloatType, IntShiftType> const&
          am_i = mappings_const_ref_[pair.i_seq][0];
        asu_mapping<FloatType, IntShiftType> const&
          am_j = mappings_const_ref_[pair.j_seq][pair.j_sym];
        sgtbx::rt_mx const& rt_i = space_group_ops_const_ref_[am_i.i_sym_op()];
        sgtbx::rt_mx const& rt_j = space_group_ops_const_ref_[am_j.i_sym_op()];
        CCTBX_ASSERT(rt_i.r().den() == rt_j.r().den()
                  && rt_i.t().den() == rt_j.t().den());
        if (rt_i.r().num() != rt_j.r().num()) return false;
        scitbx::vec3<IntShiftType> const& u_i = am_i.unit_shifts();
        scitbx::vec3<IntShiftType> const& u_j = am_j.unit_shifts();
        int t_den = rt_i.t().den();
        for(unsigned i=0;i<3;i++) {
          if (   rt_i.t().num()[i] + u_i[i] * t_den
              != rt_j.t().num()[i] + u_j[i] * t_den) return false;
        }
        return true;
      }

      //! Returns a new pair after checking the indices.
      asu_mapping_index_pair
      make_trial_pair(unsigned i_seq, unsigned j_seq, unsigned j_sym) const
      {
        CCTBX_ASSERT(mappings_const_ref_.begin() == mappings_.begin());
        CCTBX_ASSERT(i_seq < mappings_const_ref_.size());
        CCTBX_ASSERT(j_seq < mappings_const_ref_.size());
        CCTBX_ASSERT(j_sym < mappings_const_ref_[j_seq].size());
        asu_mapping_index_pair new_pair;
        new_pair.i_seq = i_seq;
        new_pair.j_seq = j_seq;
        new_pair.j_sym = j_sym;
        return new_pair;
      }

      /*! Returns a new pair after checking the indices and asserting
          pair.is_active().
       */
      asu_mapping_index_pair
      make_pair(unsigned i_seq, unsigned j_seq, unsigned j_sym) const
      {
        asu_mapping_index_pair new_pair = make_trial_pair(i_seq, j_seq, j_sym);
        CCTBX_ASSERT(new_pair.is_active());
        return new_pair;
      }

      //! Determination of the index i_sym corresponding to the given rt_mx.
      /*! The result value i_sym satisfies the relation:

               mappings()[i_seq][i_sym] * special_op(i_seq)
            == rt_mx * special_op(i_seq)

          The result value is -1 if the site corresponding to rt_mx is not
          in the asymmetric unit or the surrounding buffer region.
       */
      int
      find_i_sym(unsigned i_seq, sgtbx::rt_mx const& rt_mx) const
      {
        CCTBX_ASSERT(i_seq < mappings_const_ref_.size());
        std::size_t i_sst = site_symmetry_table_.indices_const_ref()[i_seq];
        if (i_sst == 0) {
          sgtbx::rt_mx rt_mx_sp = rt_mx.cancel();
          for(int i_sym=0; i_sym<mappings_const_ref_[i_seq].size(); i_sym++) {
            if (get_rt_mx(i_seq, i_sym).cancel() == rt_mx_sp) {
              return i_sym;
            }
          }
        }
        else {
          sgtbx::rt_mx const& special_op
            = site_symmetry_table_.table_const_ref()[i_sst].special_op();
          sgtbx::rt_mx rt_mx_sp = rt_mx.multiply(special_op);
          for(int i_sym=0; i_sym<mappings_const_ref_[i_seq].size(); i_sym++) {
            if (get_rt_mx(i_seq, i_sym).multiply(special_op) == rt_mx_sp) {
              return i_sym;
            }
          }
        }
        return -1;
      }

    protected:
      sgtbx::space_group space_group_;
      float_asu<FloatType> asu_;
      FloatType buffer_thickness_;
      af::shared<sgtbx::rt_mx> space_group_ops_;
      af::const_ref<sgtbx::rt_mx> space_group_ops_const_ref_;
      float_asu<FloatType> asu_buffer_;
      scitbx::math::sphere_3d<FloatType> buffer_covering_sphere_;
      sgtbx::site_symmetry_table site_symmetry_table_;
      array_of_array_of_mappings_for_one_site mappings_;
      af::const_ref<array_of_mappings_for_one_site> mappings_const_ref_;
      std::size_t n_sites_in_asu_and_buffer_;
      cartesian<FloatType> mapped_sites_min_;
      cartesian<FloatType> mapped_sites_max_;
      mutable std::vector<scitbx::mat3<FloatType> > r_cart_;
      mutable std::vector<scitbx::vec3<FloatType> > t_cart_;
      mutable std::vector<scitbx::mat3<FloatType> > r_inv_cart_;
  };

}}} // namespace cctbx::crystal::direct_space_asu

#endif // CCTBX_CRYSTAL_DIRECT_SPACE_ASU_H
