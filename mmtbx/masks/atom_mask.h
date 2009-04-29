#ifndef MMTBX_MASKS_ATOM_MASK_H
#define MMTBX_MASKS_ATOM_MASK_H

#include <mmtbx/error.h>
#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/space_group_type.h>
#include <cctbx/sgtbx/direct_space_asu/proto/direct_space_asu.h>
#include <cctbx/sgtbx/direct_space_asu/proto/small_vec_math.h>
#include <cctbx/miller.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_interval_grid.h>

namespace mmtbx {

  //! Masks for bulk solvent modelling
  /*!
  Mask calculations are based on the following paper.
  Jiang, J.-S. & Brünger, A. T. (1994). J. Mol. Biol. 243, 100-115.
  "Protein hydration observed by X-ray diffraction. Solvation properties
  of penicillopepsin and neuraminidase crystal structures."
   */
  namespace masks {

  namespace af = scitbx::af;

  using cctbx::sgtbx::asu::direct_space_asu;

  typedef af::shared< scitbx::double3 > coord_array_t;
  typedef af::shared< double > double_array_t;
  typedef af::c_interval_grid<3> asu_grid_t;
  typedef af::c_interval_grid<3> grid_t;
  //typedef signed char data_type; this is insufficient
  typedef short data_type;
  namespace {
    BOOST_STATIC_ASSERT( std::numeric_limits<data_type>::is_integer );
  }
  typedef af::versa<data_type, grid_t > mask_array_t;
  typedef af::versa<data_type, asu_grid_t > mask_asu_array_t;
  const data_type  mark = std::min(
      std::abs(std::numeric_limits<data_type>::max()-1),
      std::abs(std::numeric_limits<data_type>::min()+1)
      );

  //! Flat solvent mask and structure factors.
  /*! \class atom_mask atom_mask.h mmtbx/masks/atom_mask.h
   This mask will have valid values only inside asymmetric unit. Outside
   they will be zero. Inside the asu 0 - means macromolecule,
   non-zero - solvent.
   */
  class atom_mask
  {
    public:
      //! Allocates memory for mask calculation.
      /*! gridding_n_real must be compatible with the space_group.
       * At least the following must be true: every grid point has
       * to be transformed into grid point by every symmetry operator
       * of the group.
       */
      atom_mask(
        const cctbx::uctbx::unit_cell & unit_cell,
        const cctbx::sgtbx::space_group &group_,
        const grid_t::index_type & gridding_n_real,
        double solvent_radius_,
        double shrink_truncation_radius_)
      :
        solvent_radius(solvent_radius_),
        shrink_truncation_radius(shrink_truncation_radius_),
        accessible_surface_fraction(-1.0),
        contact_surface_fraction(-1.0),
        asu(group_.type()),
        cell(unit_cell),
        group(group_)
      {
        MMTBX_ASSERT( mark > group.order_z() && -mark < -group.order_z() );
        MMTBX_ASSERT(solvent_radius >= 0.0);
        MMTBX_ASSERT(shrink_truncation_radius >= 0.0);
        MMTBX_ASSERT(gridding_n_real.const_ref().all_gt(0));
        this->full_cell_grid_size = scitbx::int3(gridding_n_real);
        this->determine_boundaries(); // also allocates memory
        // TODO: ? put more computation here, eg.: mask_asu
      }

      //! Allocates memory for mask calculation.
      /*! This constructor will calculate grid size appropriate for the
       * spacegroup based on resolution and grid_step_factor.
       */
      atom_mask(
        const cctbx::uctbx::unit_cell & unit_cell,
        const cctbx::sgtbx::space_group &group_,
        double resolution,
        double grid_step_factor = 4.0,
        double solvent_radius_ = 1.11,
        double shrink_truncation_radius_ = 0.9)
      :
        solvent_radius(solvent_radius_),
        shrink_truncation_radius(shrink_truncation_radius_),
        accessible_surface_fraction(-1.0),
        contact_surface_fraction(-1.0),
        asu(group_.type()),
        cell(unit_cell),
        group(group_),
        asu_atoms()
      {
        MMTBX_ASSERT( mark > group.order_z() && -mark < -group.order_z() );
        MMTBX_ASSERT(solvent_radius >= 0.0);
        MMTBX_ASSERT(shrink_truncation_radius >= 0.0);
        this->determine_gridding(this->full_cell_grid_size, resolution,
          grid_step_factor);
        this->determine_boundaries(); // also allocates memory
        // TODO: ? put more computation here, eg.: mask_asu
      }

      //! Clears current, and calculates new mask based on atomic data.
      void compute(
        const coord_array_t & sites_frac,
        const double_array_t & atom_radii
      );

      const mask_array_t & get_mask() const { return data; }

      //! Computes mask structure factors.
      scitbx::af::shared< std::complex<double> > structure_factors(
        const scitbx::af::const_ref< cctbx::miller::index<> > &indices );

      scitbx::int3 grid_size() const { return full_cell_grid_size; }
      size_t grid_size_1d() const
      {
        MMTBX_ASSERT( scitbx::ge_all(this->grid_size(), scitbx::int3(0,0,0)) );
        return static_cast<size_t>(this->grid_size()[0])
          * static_cast<size_t>(this->grid_size()[1])
          * static_cast<size_t>(this->grid_size()[2]);
      }

      size_t n_asu_atoms() const
      {
        return asu_atoms.size();
      }

      const double solvent_radius;
      const double shrink_truncation_radius;

      //! Solvent volume, if atom radius = vdw_radius + solvent_radius
      double accessible_surface_fraction;

      //! Solvent volume
      double contact_surface_fraction;

      // execution time in millisecs, DO NOT USE
      long debug_mask_asu_time, debug_atoms_to_asu_time,
           debug_accessible_time, debug_contact_time, debug_fft_time;

      bool debug_has_enclosed_box;

    private:

      typedef std::pair< scitbx::double3, double > atom_t;
      typedef std::vector< atom_t > atom_array_t;
      void compute_contact_surface();
      void compute_accessible_surface(const atom_array_t &atoms);
      void mask_asu();
      void atoms_to_asu(
        const coord_array_t & sites_frac,
        const double_array_t & atom_radii);

      void determine_gridding(cctbx::sg_vec3 &grid, double resolution,
        double factor = 4.0) const;
      void get_asu_boundaries(scitbx::int3 &low, scitbx::int3 &high) const;
      void get_expanded_asu_boundaries(scitbx::double3 &low,
        scitbx::double3 &high) const;
      void get_expanded_asu_boundaries(scitbx::int3 &low,
        scitbx::int3 &high) const;
      void determine_boundaries();

      // these 3 should be some kind of safe reference
      const direct_space_asu asu;
      const cctbx::uctbx::unit_cell cell;
      const cctbx::sgtbx::space_group group;

      scitbx::int3 full_cell_grid_size, asu_low, asu_high;
      scitbx::double3 expanded_box[2];
      atom_array_t asu_atoms;
      mask_array_t data;
  }; // class atom_mask

}} // namespace mmtbx::masks
#endif
