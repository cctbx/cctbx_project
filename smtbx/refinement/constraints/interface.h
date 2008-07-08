/** Warning! The inclusion of this header will result in compilation error.
It is only aimed at documenting once and for all the interface all constraints
abide to.
*/

namespace smtbx { namespace refinement { namespace constraints {

/// A set of constraints on some parameters of several scatterers
/** One important class abiding to this interface is
    \code smtbx::refinement::constraints::constraint_array \endcode.
    Concerning the latter class, see also
    \code few_scatterer_constraints \endcode .
*/
template<typename FloatType, class XrayScattererType,
         template<class> class SharedArray1D=af::shared>
class many_scatterer_constraints
{
  public:
    typedef FloatType float_type;
    typedef XrayScattererType xray_scatterer_type;
    typedef parameter_map<xray_scatterer_type> parameter_map_type;

    /// Construct a constraint
    /** The constraint is on the structure defined by the given unit cell,
        site symmetry table and array of scatterers. The client code must also
        have constructed a parameter map for that structure and pass it to
        this constructor. Although the information is the latter is already
        in the array of scatterers, the parameter map is best constructed only
        once for efficiency reason.

        This object keeps a reference to all those objects, which explains
        that the other member functions do not need to be passed the information
        about the structure.

        The argument \code constraint_flags \endcode is used as follow:
        \code constraint_flags[i] \endcode must be a flag summarizing the
        constraints on the i-th scatterer at the time of construction of this object. E.g. \code constraint_flags[i].grad_site() \endcode being true
        means that the site of the i-th scatterer is still free
        whereas its being false means it has already been constrained
    */
    many_scatterer_constraints(
      uctbx::unit_cell const &unit_cell_,
      sgtbx::site_symmetry_table const &site_symmetry_table_,
      af::shared<xray_scatterer_type> scatterers_,
      parameter_map_type const &crystallographic_parameter_map_,
      af::ref<xray::scatterer_flags> constraint_flags);

    /// The constraints already active when this object was constructed.
    /** This map associates to each scatterer index the corresponding
        constraint flags are they were passed to the constructor.
    */
    std::map<int, xray::scatterer_flags> already_constrained();

    /// Compute the gradients using the chain rule in the reverse mode
    /** \arg crystallographic_gradients the array of gradients
        wrt crystallographic parameters, as obtained from
        \code cctbx::xray::structure_factors::gradients_direct<>::packed()
        \endcode for example. The application of the chain rule may result
        in the elements of that arguments being modified.

        \arg reparameterization_gradients if the constraint requires the
        introduction of extra parameters, then the derivatives wrt those will
        be appended to that array.

    */
    void compute_gradients(
      af::ref<float_type> const &crystallographic_gradients,
      SharedArray1D<float_type> reparameterization_gradients);

    /// Apply the shifts of crystallographic and extra parameters to the
    /// scatterers in the structure.
    /** The shifts of the extra parameters must be in the exact same order
        as the derivates were filed in by \code compute_gradients \endcode.
        This is essential so that this object knows where to look for those
        shifts it is concerned with (as opposed to those shifts which other
        constraints object are concerned with).
    */
    void apply_shifts(
      af::const_ref<float_type> const &crystallographic_shifts,
      af::const_ref<float_type> const &reparametrization_shifts);

    /// Place the scatterers where they should be according to the constraints
    void place_constrained_scatterers();
};


/// Constraints on some parameters of a few scatterers
/** Instances of classes abiding to this interface are to be used as
    elements of \code smtbx::refinement::constraints::constraint_array \endcode
    Examples are the geometrical hydrogen's constraints.

    Its member function
    \code compute_gradients \endcode,
    \code apply_shifts \endcode,
    \code place_constrained_scatters \endcode
    perform the same task as the member functions with the same names in
    \code many_scatterer_constraints \endcode. The only difference is that
    this class follows a lightweight pattern and therefore those member
    functions are passed the necessary context, that the
    \code constraint_array \endcode
    an instance belongs to passes to it. That pattern also requires an
    after-construction initialisation with the necessary context, which is done
    by \code initialise_in_context \endcode.

*/
template<typename FloatType, class XrayScattererType,
         template<class> class SharedArray1D=af::shared>

class few_scatterer_constraints
{
  public:
    typedef XrayScattererType xray_scatterer_type;
    typedef FloatType float_type;
    typedef parameter_map<xray_scatterer_type> parameter_map_type;
    typedef cartesian<float_type> cart_t;
    typedef fractional<float_type> frac_t;

    /// Initialise the object in the given context
    /** This method is called when this object is added to a
        \code constraint_array \endcode. The meaning of the arguments
        is the same as for the constructor of \code many_scatterer_constraints
        \endcode. The only difference is the presence of
        \code already_constrained \endcode where this object should record
        which and how scatterers may have already been constrained.
    */
    void initialise_in_context(
      uctbx::unit_cell const &unit_cell,
      af::const_ref<xray_scatterer_type> const &scatterers,
      af::ref<xray::scatterer_flags> const &constraint_flags,
      std::map<int, xray::scatterer_flags> &already_constrained);

    /// Called when the \code constraint_array \endcode it belong to
    /// runs its member function \code compute_gradients \endcode.
    void compute_gradients(
      uctbx::unit_cell const &unit_cell,
      sgtbx::site_symmetry_table const &site_symmetry_table,
      af::const_ref<xray_scatterer_type> const &scatterers,
      parameter_map_type const &crystallographic_parameter_map,
      af::ref<float_type> const &crystallographic_gradients,
      SharedArray1D<float_type> reparametrization_gradients);

    /// Called when the \code constraint_array \endcode it belong to
    /// runs its member function \code apply_shifts \endcode.
    void apply_shifts(
      uctbx::unit_cell const &unit_cell,
      sgtbx::site_symmetry_table const &site_symmetry_table,
      af::ref<xray_scatterer_type> const &scatterers,
      parameter_map_type const &crystallographic_parameter_map,
      af::const_ref<float_type> const &crystallographic_shifts,
      af::const_ref<float_type> const &reparametrization_shifts);

    /// Called when the \code constraint_array \endcode it belong to
    /// runs its member function \code place_constrained_scatterers \endcode.
    void place_constrained_scatterers(
      uctbx::unit_cell const &unit_cell,
      sgtbx::site_symmetry_table const &site_symmetry_table,
      af::ref<xray_scatterer_type> const &scatterers);
};


}}}
