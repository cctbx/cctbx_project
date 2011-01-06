#ifndef SMTBX_REFINEMENT_CONSTRAINTS_REPARAMETRISATION_H
#define SMTBX_REFINEMENT_CONSTRAINTS_REPARAMETRISATION_H

#include <scitbx/array_family/ref.h>
#include <scitbx/sparse/matrix.h>
#include <smtbx/import_scitbx_af.h>
#include <cctbx/uctbx.h>
#include <cctbx/xray/scatterer.h>
#include <smtbx/import_cctbx.h>
#include <smtbx/error.h>

#include <boost/shared_array.hpp>
#include <boost/foreach.hpp>
#include <boost/range/iterator_range.hpp>
#include <vector>
#include <iterator>
#include <queue>
#include <iosfwd>

/// Constraints handled as reparametrisation
/** Constraints are in the first place restrictions placed on quantities
    depending on the refined crystallographic parameters. For example,
    a tetrahedral C-CH3 is modelled by stating that bond lengths and bond
    angles, which are functions of the sites of the Hydrogen atoms
    and of the Carbon they are bonded to, shall keep some fixed values.

    Traditionally, those constraints are solved analytically by expressing
    some of the parameters involved in the constraints as functions
    of the others, and sometimes of extra parameter, which are usually referred
    to as "free variables" following the ShelX nomenclature.
    In the above example, the Hydrogen sites are a function of the Carbon site
    and of one rotation angle around the C-C bond.

    The former presentation is what is known as constraints in numerical
    optimisation, whereas the latter is how that word is understood in
    crystallography.
 */
namespace smtbx { namespace refinement { namespace constraints {


/// Graph traversal classic colouring
//@{
unsigned const white=1, grey=2, black=3;
//@}

/// Type of sparse matrix used throughout the module
typedef scitbx::sparse::matrix<double> sparse_matrix_type;

/// Cartesian and fractional coordinates used throughout the module
//@{
typedef cartesian<double> cart_t;
typedef fractional<double> frac_t;
//@}

/// Anisotropic displacement tensor used throughout the module
typedef scitbx::sym_mat3<double> tensor_rank_2_t;


/// A range of index [fist, first+size)
/** The member valid tells whether the members first and size have any meaning
    at all. An invalid instance is mapped to an empty tuple in Python,
    and to a tuple (first, size) otherwise.
 */
class index_range
{
public:
  bool is_valid() const { return valid; }
  std::size_t first() const { return first_; }
  std::size_t last() const { return first_ + size_; }
  std::size_t size() const { return size_; }

  index_range()
  : valid(false)
  {}

  index_range(std::size_t first, std::size_t size)
  : valid(true), first_(first), size_(size)
  {}

private:
  bool valid;
  std::size_t first_, size_;
};


/// A parameter which may be vector-valued
/** Its components may be independent scalar parameters eventually passed
    to the minimisation engine, or it may be function of other scalar- or
    vector-valued parameters, which we will refer to as its arguments.

    Implementation note: this class will be inherited virtually to allow
    multiple inheritance with a correct handling of the dreaded diamond.
    As a result, the most derived class will need to *explicitly*
    call the constructor of this class parameter. Thus the usual pattern
    of calling this class constructor in the constructor of direct heirs
    is redundant and we won't do it. But then the compiler will try to
    invoke a constructor for class parameter without argument. Thus we
    provide a trivial one. The pattern repeats itself down the inheritance
    hierarchy: either a class defines no constructor (and gets a compiler-
    generated default constructor without arguments) or it defines a non-trivial
    constructor, in which case a trivial one needs be provided too. In any case,
    any non-trivial constructor shall only initialise the class member and
    not attempt to initialise the virtual bases higher in the inheritance
    hierarchy since, again, that will have to be done by the most derived
    classes anyway. (phew!)
 */
class parameter
{
public:
  template <class>
  friend class dfs_visitor;

  friend class variability_visitor;

  template <class>
  friend class topologist;

  template <class>
  friend class parameter_tagger;

  friend class reparametrisation;

  /// Trivial constructor made necessary by virtual inheritance pattern we use.
  /** See documentation header for this class */
  parameter() {}

  /// Construct a parameter with the given number of arguments
  parameter(std::size_t n_arguments)
    : colour_(white), variable(true), root(false), n_args(n_arguments),
      arg(new parameter *[n_arguments])
  {}

  virtual ~parameter();

  /// This vector parameter has its components indexed in the range
  /// [ index(), index() + size() ).
  /** This is used to label rows and columns of the Jacobian matrix
   */
  std::size_t index() const { return index_; }

  /// Number of arguments
  std::size_t n_arguments() const { return n_args; }

  /// Whether this parameter is independent in the refinement
  bool is_independent() const { return !n_arguments(); }

  /// i-th argument
  parameter *argument(std::size_t i) const { return arg[i]; }

  //@{
  /// Colour accessor
  unsigned colour() const { return colour_; }

  void set_colour(unsigned c) { colour_ = c; }
  //@}

  /// Whether this parameter is a root in the dependency tree it belongs to
  bool is_root() const { return root; }

  /// Whether this parameter will be refined by the minimisation engine
  virtual bool is_variable() const;

  /// Dimension of this vector parameter
  std::size_t size() const {
    return components().size();
  }

  /// A reference to the vector of component values
  virtual af::ref<double> components() = 0;

  /// A reference to the vector of component values
  af::const_ref<double> components() const {
    return const_cast<parameter *>(this)->components();
  }

  /// Compute the value of this parameter but not its Jacobian
  void evaluate(uctbx::unit_cell const &unit_cell) {
    return linearise(unit_cell, 0);
  }

  /// Evaluate value and Jacobian
  /** It shall use the forward method,
      assuming its arguments have already been linearised.
      Moreover if jacobian_transpose == 0, then only the value shall
      be computed.
   */
  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose) = 0;

protected:
  /// This is to be used by heirs to set the arguments at construction
  //@{
  void set_argument(std::size_t i, parameter *p) {
    SMTBX_ASSERT(p);
    arg[i] = p;
  }

  void set_arguments(parameter *p0) {
    set_argument(0, p0);
  }
  void set_arguments(parameter *p0, parameter *p1) {
    set_argument(0, p0);
    set_argument(1, p1);
  }
  void set_arguments(parameter *p0, parameter *p1, parameter *p2) {
    set_argument(0, p0);
    set_argument(1, p1);
    set_argument(2, p2);
  }
  void set_arguments(parameter *p0, parameter *p1, parameter *p2, parameter *p3)
  {
    set_argument(0, p0);
    set_argument(1, p1);
    set_argument(2, p2);
    set_argument(3, p3);
  }
  void set_arguments(parameter *p0, parameter *p1, parameter *p2, parameter *p3,
                     parameter *p4)
  {
    set_argument(0, p0);
    set_argument(1, p1);
    set_argument(2, p2);
    set_argument(3, p3);
    set_argument(4, p4);
  }
  //@}

  /// This is to be used by heirs and friends to set the variability status
  virtual void set_variable(bool f);

  /// This is to be used by heirs and friends to set the root status
  void set_root(bool f) { root = f; }

  /// Set the index of this parameter
  void set_index(std::size_t idx) { index_ = idx; }

private:
  /* Don't use bit fields to reduce memory footprint any more
     because it triggered a nasty bug when index_ >= 256 at least
     on MacOS X with gcc 4.2 and clang.
   */
  unsigned char colour_ ;
  bool          variable;
  bool          root    ;
  unsigned char n_args  ;
  std::size_t   index_  ;
  parameter     **arg   ;
};


/// Exception carrying the parameter at fault
class error : public smtbx::error
{
public:
  error(std::string const &msg, parameter *p)
    : smtbx::error(msg), p(p)
  {
    std::ostringstream o;
    o << " parameter at address " << std::hex << p << ".";
    msg_ += o.str();
  }

private:
  parameter *p;
};


/// scalar parameter
class scalar_parameter : public virtual parameter
{
public:
  virtual af::ref<double> components();

  double value;
};


/// Independent scalar parameter
class independent_scalar_parameter : public scalar_parameter
{
public:
  independent_scalar_parameter(double value, bool variable=true)
  : parameter(0)
  {
    this->value = value;
    set_variable(variable);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
};


template <int N>
class small_vector_parameter : public virtual parameter
{
public:
  virtual af::ref<double> components() { return value.ref(); }

  af::small<double, N> value;
};


template <int N>
class independent_small_vector_parameter
  : public small_vector_parameter<N>
{
public:
  independent_small_vector_parameter(af::small<double, N> const &value,
                                     bool variable=true)
  : parameter(0)
  {
    this->value = value;
    this->set_variable(variable);
  }

  /// Construct with an intial value equal to the zero vector of dimension n
  independent_small_vector_parameter(int n, bool variable=true)
  : parameter(0)
  {
    this->value.resize(n, 0);
    this->set_variable(variable);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose)
  {}

};


/// Site, isotropic or anisotropic displacement, etc., or combination of those.
/** They belong to one or more scatterers in the asymmetric unit.

    Heirs of this class implementing behaviour for specific scatterer parameters
    (e.g. site, occupancy, etc) shall read their values from their associated
    scatterer objects only once, at construction time. They shall hold the
    values internally and then write those back to the scatterer when
    member function store() is called but all reparametrisation computation
    are performed with the internal values. Consequently changing the values
    of the scatterer parameters directly (e.g. sc.site = ...) in user code
    will have no effect on the reparametrisation computation.
 */
class asu_parameter : public virtual parameter
{
public:
  typedef xray::scatterer<> scatterer_type;
  typedef af::const_ref<scatterer_type *> scatterer_sequence_type;

  /// The scatterers in the asu this models some parameters of.
  virtual scatterer_sequence_type scatterers() const = 0;

  /// Indices of those components associated with the given asu scatterer.
  /** The resulting range may be invalid if the scatterer is not one of those
      this parameter refers to.
   */
  virtual index_range
  component_indices_for(scatterer_type const *scatterer) const = 0;

  /// Write annotations for each component associated with the given scatterer
  /** It shall be written as a comma separated list which shall be empty
      if the scatterer is not one of those this parameter refers to.
   */
  virtual void
  write_component_annotations_for(scatterer_type const *scatterer,
                                  std::ostream &output) const = 0;

  /// Store its components into the corresponding scatterers
  virtual void store(uctbx::unit_cell const &unit_cell) const = 0;
};


/// Parameter of a single scatterer in the asu
class single_asu_scatterer_parameter : public virtual asu_parameter
{
public:
  single_asu_scatterer_parameter() {}

  single_asu_scatterer_parameter(scatterer_type *scatterer)
  : scatterer(scatterer)
  {}

  virtual scatterer_sequence_type scatterers() const;

  virtual index_range
  component_indices_for(scatterer_type const *scatterer) const;

protected:
  /// The scatterer this parameter belongs to
  scatterer_type *scatterer;
};


/// Scatterer site.
/** A parameter whose components are the fractional coordinates.
 */
class site_parameter : public virtual parameter
{
public:
  virtual af::ref<double> components();

  /// The site value in Cartesian coordinates
  fractional<double> value;
};



/// asu scatterer site.
/** A parameter whose components are the fractional coordinates of a
    scatterer in the asu.
 */
class asu_site_parameter : public site_parameter,
                           public virtual single_asu_scatterer_parameter
{
public:
  /// Variability property, directly linked to the scatterer grad_site flag
  //@{
  virtual void set_variable(bool f);

  virtual bool is_variable() const;
  //@}

  virtual void
  write_component_annotations_for(scatterer_type const *scatterer,
                                  std::ostream &output) const;

  virtual void store(uctbx::unit_cell const &unit_cell) const;

};


/// Scatterer site that is an independent parameter
class independent_site_parameter : public asu_site_parameter
{
public:
  independent_site_parameter(scatterer_type *scatterer)
  : parameter(0),
    single_asu_scatterer_parameter(scatterer)
  {
    value = scatterer->site;
  }

  /// Does nothing in this class
  /** This optimisation relies on class reparametrisation implementation
   */
  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
};


/// Anisotropic displacement parameters of a site
/** A parameter whose components are the coefficients of the tensor
    in fractional coordinates.
 */
class u_star_parameter : public virtual parameter
{
public:
  virtual af::ref<double> components();

  /// The site value in Cartesian coordinates
  tensor_rank_2_t value;
};


/// Anisotropic displacement parameters of a site in the asu
class asu_u_star_parameter : public u_star_parameter,
                             public virtual single_asu_scatterer_parameter
{
public:
  /// Variability property, directly linked to the scatterer grad_site flag
  //@{
  virtual void set_variable(bool f);

  virtual bool is_variable() const;
  //@}

  virtual void
  write_component_annotations_for(scatterer_type const *scatterer,
                                  std::ostream &output) const;

  virtual void store(uctbx::unit_cell const &unit_cell) const;

};


class independent_u_star_parameter : public asu_u_star_parameter
{
public:
  independent_u_star_parameter(scatterer_type *scatterer)
  : parameter(0),
    single_asu_scatterer_parameter(scatterer)
  {
    value = scatterer->u_star;
  }

  /// Does nothing in this class
  /** This optimisation relies on class reparametrisation implementation
   */
  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
};


/// Occupancy of a scatterer
class asu_occupancy_parameter : public scalar_parameter,
                                public virtual single_asu_scatterer_parameter
{
public:
  /// Variability property, directly linked to the scatterer grad_occupancy flag
  //@{
  virtual void set_variable(bool f);

  virtual bool is_variable() const;
  //@}

  virtual void
  write_component_annotations_for(scatterer_type const *scatterer,
                                  std::ostream &output) const;

  virtual void store(uctbx::unit_cell const &unit_cell) const;
};


class independent_occupancy_parameter : public asu_occupancy_parameter
{
public:
  independent_occupancy_parameter(scatterer_type *scatterer)
  : parameter(0),
    single_asu_scatterer_parameter(scatterer)
  {
    value = scatterer->occupancy;
  }

  /// Does nothing in this class
  /** This optimisation relies on class reparametrisation implementation
   */
  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
};


/// Isotropic thermal displacement parameter of a scatterer
class asu_u_iso_parameter : public scalar_parameter,
                            public virtual single_asu_scatterer_parameter
{
public:
  /// Variability property, directly linked to the scatterer grad_u_iso flag
  //@{
  virtual void set_variable(bool f);

  virtual bool is_variable() const;
  //@}

  virtual void
  write_component_annotations_for(scatterer_type const *scatterer,
                                  std::ostream &output) const;

  virtual void store(uctbx::unit_cell const &unit_cell) const;
};


class independent_u_iso_parameter : public asu_u_iso_parameter
{
public:
  independent_u_iso_parameter(scatterer_type *scatterer)
  : parameter(0),
    single_asu_scatterer_parameter(scatterer)
  {
    value = scatterer->u_iso;
  }

  /// Does nothing in this class
  /** This optimisation relies on class reparametrisation implementation
   */
  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
};


class computing_graph_has_cycle_error : public error
{
public:
  computing_graph_has_cycle_error(parameter *p)
    : error("Cycle detected in constraints computing graph at", p)
  {}
};

/// A DFS visit following from each parameter to its arguments.
/** This uses a CRTP pattern.
 */
template <class DerivedType>
class dfs_visitor
{
public:
  void visit(parameter *p) {
    if (p->colour() == white) {
      heir()->start(p);
      if (heir()->shall_start_visit_from(p)) visit_from(p);
    }
  }

private:
  DerivedType *heir() { return static_cast<DerivedType *>(this); }

  void visit_from(parameter *p) {
    heir()->discover(p);
    p->set_colour(grey);
    for (int i=0; i<p->n_arguments(); ++i) {
      parameter *q = p->argument(i);
      heir()->examine_edge(p, q);
      if (!shall_cross(p, q)) continue;
      if (q->colour() == white) visit_from(q);
      else if(q->colour() == grey) throw computing_graph_has_cycle_error(q);
    }
    heir()->finish(p);
    p->set_colour(black);
  }

  /// Callbacks overridable by heirs
  //@{

  /// Called with the parameter the visit starts from
  /** The default is to do nothing */
  void start(parameter *p) { }

  /// Whether the visit shall proceed by starting a DFS from p
  /** The default is to do so */
  bool shall_start_visit_from(parameter *p) { return true; }

  /// Called just after p is first seen
  /** The default is to do nothing */
  void discover(parameter *p) { }

  /// Call just after the edge from p to its argument q is discovered
  /** the default is to do nothing */
  void examine_edge(parameter *p, parameter *q) { }

  /// Whether the visit shall proceed by crossing from p to its argument q.
  /** The default is to do so. */
  bool shall_cross(parameter *p, parameter *q) { return true; }

  /// Called just before p is seen for the last time.
  /** This is after the entire subtree hanging from p has been explored.
      The default is to do nothing
   */
  void finish(parameter *p) { }

  //@}
};


class variability_visitor : public dfs_visitor<variability_visitor>
{
public:
  variability_visitor(uctbx::unit_cell const &unit_cell)
    : unit_cell(unit_cell)
  {}

  void finish(parameter *p) {
    if (p->n_arguments()) {
      p->set_variable(false);
      for (int i=0; i<p->n_arguments(); ++i) {
        parameter *q = p->argument(i);
        if (q->is_variable()) p->set_variable(true);
      }
    }
    if (!p->is_variable()) p->evaluate(unit_cell);
  }

private:
  uctbx::unit_cell const &unit_cell;
};


class evaluator : public dfs_visitor<evaluator>
{
public:
  evaluator(uctbx::unit_cell const &unit_cell,
            sparse_matrix_type *jacobian_transpose)
    : unit_cell(unit_cell), jacobian_transpose(jacobian_transpose)
  {}

  bool shall_start_visit_from(parameter *p) { return p->is_variable(); }

  bool shall_cross(parameter *p, parameter *q) { return q->is_variable(); }

  void finish(parameter *p) { p->linearise(unit_cell, jacobian_transpose); }

private:
  uctbx::unit_cell const &unit_cell;
  sparse_matrix_type *jacobian_transpose;
};


template <class OutputPointerType>
class topologist : public dfs_visitor<topologist<OutputPointerType> >
{
public:
  topologist(OutputPointerType all)
    : all(all)
  {}

  void start(parameter *p) { p->set_root(true); }

  void examine_edge(parameter *p, parameter *q) { q->set_root(false); }

  void finish(parameter *p) { *all++ = p; }

private:
  OutputPointerType all;
};



class cant_modify_parameter_value_error : public error
{
public:
  cant_modify_parameter_value_error(parameter *p)
    : error("No access to value for independent variable",
            p)
  {}
};

/// Crystal structure reparametrisation
/** This holds parameters dependent on each others, most of them being
    crystallographic parameters, some of which being roots of the computational
    graph.
 */
class reparametrisation
{
private:
  typedef std::vector<parameter *> parameter_array_t;

public:
  typedef asu_parameter::scatterer_type scatterer_type;

  typedef parameter_array_t::iterator iterator;

  typedef boost::iterator_range<iterator> range;

  /// Construction in one go from parameters from the range [first, last)
  /// and all their direct or indirect arguments.
  /** This object becomes the owner of all those parameters and
      they will therefore be deallocated when this object is deleted.
   */
  template <class ForwardIteratorType>
  reparametrisation(uctbx::unit_cell const &unit_cell,
                    boost::iterator_range<ForwardIteratorType> const &params)
    : unit_cell(unit_cell)
  {
    // Classify parameters
    typedef std::back_insert_iterator<std::vector<parameter *> >
            all_param_inserter_t;
    topologist<all_param_inserter_t> t(std::back_inserter(all));
    BOOST_FOREACH(parameter *p, params) t.visit(p);
    whiten(); // only time we need to call that explicitly

    analyse_variability();
  }

  /// Incremental construction
  //@{

  /// Construct an object without any reparametrisation
  reparametrisation(uctbx::unit_cell const &unit_cell)
  : unit_cell(unit_cell)
  {}

  /// Add a new parameter and all its direct or indirect arguments
  /** One shall call finalise() after the desired parameters have been added.

      This object takes ownership of those parameters, which will therefore
      be deallocated when this object is destroyed.
   */
  void add(parameter *p);

  /// Ready this for linearise(), etc.
  void finalise();
  //@}

  /// Walks the computational graph to find constant branches
  /** This member function is to be called every time the variability of
      a parameter changes after this has been constructed
   */
  void analyse_variability();

  /// Destroy all parameters in the computational graph.
  ~reparametrisation();

  /// The range of all parameters held by this object
  range parameters() { return boost::make_iterator_range(all); }

  /// Total number of parameter components
  std::size_t n_components() {
    return n_independents() + n_intermediates() + n_non_trivial_roots();
  }

  /// Number of independent parameter components
  std::size_t n_independents() { return n_independents_; }

  /// Number of components of parameters which are neither independents or roots
  std::size_t n_intermediates() { return n_intermediates_; }

  /// Number of components of roots which are not independent parameters
  std::size_t n_non_trivial_roots() { return n_non_trivial_roots_; }

  /// Call parameter::linearise() on each computational graph vertex.
  /** In the right order. */
  void linearise();

  /// Apply the given shifts to the independent parameters
  /** Indexing is, as for the Jacobian, enforced by parameter::index()
   */
  void apply_shifts(af::const_ref<double> const &shifts);

  /// Norm of the vector of independent parameters
  double norm_of_independent_parameter_vector();

  /// Store all crystallographic parameter values into their respective
  /// scatterers.
  void store();

  /// Let the given visitor visit all parameters.
  /** All nodes are whitened before returning, allowing another visit to proceed
      immediately.
   */
  template <class Visitor>
  void accept(Visitor &v) {
    BOOST_FOREACH(parameter *p, all) if (p->is_root()) v.visit(p);
    whiten();
  }

private:
  void whiten();

public:
  /// The transpose of the Jacobian of the function transforming independent
  /// parameters into root parameters.
  /** Mathematically, \f$[ \partial x_j/\partial x_i ]_ij\f$
   where \f$x_i\f$ is the parameter p s.t. p.index() == i.
   With p, q, r = n_independents(), n_intermediates(), n_non_trivial_roots(),
   the indexing is done as follow:
   0   .. p-1     : independent parameters
   p   .. p+q-1   : intermediate parameters
   p+q .. p+q+r-1 : those roots that are not also independent
   */
  sparse_matrix_type jacobian_transpose;

private:
  uctbx::unit_cell unit_cell;
  parameter_array_t all;
  std::size_t n_independents_, n_intermediates_, n_non_trivial_roots_;
};


}}}


#endif // GUARD
