#ifndef SMTBX_REFINEMENT_CONSTRAINTS_REPARAMETRISATION_H
#define SMTBX_REFINEMENT_CONSTRAINTS_REPARAMETRISATION_H

#include <scitbx/array_family/ref.h>
#include <scitbx/sparse/matrix.h>
#include <smtbx/import_scitbx_af.h>
#include <cctbx/uctbx.h>
#include <cctbx/xray/scatterer.h>
#include <smtbx/import_cctbx.h>
#include <smtbx/error.h>

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range.hpp>
#include <boost/functional/hash.hpp>
#include <vector>
#include <iterator>
#include <queue>

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
    In the above example, the Hydrogen sites are function of the Carbon site
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


/// A parameter which may be vector-valued
/** Its components may be independent scalar parameters eventually passed
    to the minimisation engine, or it may be function of other scalar- or
    vector-valued parameters, which we will refer to as its arguments.
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
  virtual std::size_t size() const = 0;

  /// A pointer to the vector of component values
  /** Zero in this class: it shall be overriden by heirs allowing
      the value to be accessed, merely independent parameters.
   */
  virtual double *components();

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
  void set_argument(std::size_t i, parameter *p) { arg[i] = p; }

  void set_arguments(parameter *p0) { arg[0] = p0; }
  void set_arguments(parameter *p0, parameter *p1) { arg[0] = p0; arg[1] = p1; }
  void set_arguments(parameter *p0, parameter *p1, parameter *p2) {
    arg[0] = p0; arg[1] = p1; arg[2] = p2;
  }
  void set_arguments(parameter *p0, parameter *p1, parameter *p2, parameter *p3)
  {
    arg[0] = p0; arg[1] = p1; arg[2] = p2; arg[3] = p3;
  }
  void set_arguments(parameter *p0, parameter *p1, parameter *p2, parameter *p3,
                     parameter *p4)
  {
    arg[0] = p0; arg[1] = p1; arg[2] = p2; arg[3] = p3; arg[4] = p4;
  }
  //@}

  /// This is to be used by heirs and friends to set the variability status
  virtual void set_variable(bool f);

  /// This is to be used by heirs and friends to set the root status
  void set_root(bool f) { root = f; }

  /// Set the index of this parameter
  void set_index(std::size_t idx) { index_ = idx; }

private:
  /* Use bit fields to reduce memory footprint */
  #if defined(__GNUC__) && defined(__x86_64__)
    #define SMTBX_REFINEMENT_CONSTRAINTS_INDEX_FIELD_WIDTH : sizeof(std::size_t)
  #else
    // size specs trigger crash on Win32 with VS 2008 and 32-bit RHEL 3.9
    #define SMTBX_REFINEMENT_CONSTRAINTS_INDEX_FIELD_WIDTH
  #endif
  unsigned colour_            : 2                   ;
  bool               variable : 1                   ;
  bool               root     : 1                   ;
  unsigned           n_args   : 4                   ;
  std::size_t        index_   SMTBX_REFINEMENT_CONSTRAINTS_INDEX_FIELD_WIDTH;
  parameter        **arg                            ;
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


/// An independent scalar parameter which is not a crystallographic parameter.
/** It means it is not one of those attributes features by xray::scatterer<>
    or its equivalent.
 */
class independent_scalar_parameter : public parameter
{
public:
  independent_scalar_parameter(double value, bool variable=true)
    : parameter(0), value(value)
  {
    set_variable(variable);
  }

  virtual std::size_t size() const;

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

  virtual double *components();

  double value;
};


/// Site, isotropic or anisotropic displacement, etc., or combination of those.
/** They may belong to one or more scatterers
 */
class crystallographic_parameter : public parameter
{
public:
  typedef xray::scatterer<> scatterer_type;
  typedef boost::shared_ptr<scatterer_type> scatterer_pointer;

  crystallographic_parameter(std::size_t n_arguments)
    : parameter(n_arguments)
  {}

  /// Store its components into the corresponding scatterers
  virtual void store(uctbx::unit_cell const &unit_cell) const = 0;
};


/// Scatterer site.
/** A parameter whose components are the fractional coordinates.
 */
class site_parameter : public crystallographic_parameter
{
public:
  virtual std::size_t size() const;

  site_parameter(scatterer_pointer &scatterer, std::size_t n_arguments)
    : crystallographic_parameter(n_arguments),
      scatterer(scatterer)
  {}

  virtual void store(uctbx::unit_cell const &unit_cell) const;

  /// The site value in Cartesian coordinates
  fractional<double> value;

protected:
  /// The scatterer this parameter belongs to
  scatterer_pointer scatterer;
};


/// Scatterer site that is an independent parameter
class independent_site_parameter : public site_parameter
{
public:
  independent_site_parameter(scatterer_pointer &scatterer)
    : site_parameter(scatterer, 0)
  {}

  /// Variability property, directly linked to the scatterer grad_site flag
  //@{
  virtual void set_variable(bool f);

  virtual bool is_variable() const;
  //@}

  /// Read the site value from the referenced scatterer
  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

  virtual double *components();
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
      visit_from(p);
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

  /// Called just after p is first seen
  /** The default is to do nothing */
  void discover(parameter *p) { }

  /// Call just after the edge from p to its argument q is discovered
  /** the default it to do nothing */
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
  typedef crystallographic_parameter::scatterer_type scatterer_type;

  typedef parameter_array_t::iterator iterator;

  typedef boost::iterator_range<iterator> range;

  /// Construct from parameters from the range [first, last)
  /** This object becomes the owner of all the parameter in that range and
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
    whiten(); // only time we need to call that explicitely

    analyse_variability();
  }

  /// Walks the computational graph to find constant branches
  /** This member function is to be called every time the variability of
      a parameter changes after this has been constructed
   */
  void analyse_variability() {
    /* Assign variability to each parameter.
     It also evaluates constant parameters once and for all.
     */
    variability_visitor var(unit_cell);
    accept(var);

    // Assign indices to parameters
    n_independents_ = n_intermediates_ = n_non_trivial_roots_ = 0;
    BOOST_FOREACH(parameter *p, all) {
      std::size_t s = p->size();
      if      (!p->is_variable()) n_intermediates_ += s;
      else if (!p->n_arguments()) n_independents_ += s;
      else if (p->is_root())      n_non_trivial_roots_ += s;
      else                        n_intermediates_ += s;
    }
    std::size_t i_independent = 0,
    i_intermediate = n_independents(),
    i_non_trivial_root = n_independents() + n_intermediates();
    BOOST_FOREACH(parameter *p, all) {
      std::size_t s = p->size();
      if      (!p->is_variable()) {
        p->set_index(i_intermediate);
        i_intermediate += s;
      }
      else if (!p->n_arguments()) {
        p->set_index(i_independent);
        i_independent += s;
      }
      else if (p->is_root()) {
        p->set_index(i_non_trivial_root);
        i_non_trivial_root += s;
      }
      else {
        p->set_index(i_intermediate);
        i_intermediate += s;
      }
    }

    // Initialise Jacobian transpose: [ dx_j/dx_i ]_ij
    /* The block of independent parameters is initialised to the identity matrix.
     Logically, it should be done in independent_xxxx_parameter::linearise,
     but it is more efficient to do it once and for all here.
     */
    sparse_matrix_type jt(n_independents(), n_parameters());
    for (std::size_t j=0; j<n_independents(); ++j) jt(j, j) = 1.;
    jacobian_transpose = jt;
  }

  /// Destroy all parameters in the computational graph.
  ~reparametrisation() {
    BOOST_FOREACH(parameter *p, all) delete p;
  }

  /// The range of all parameters held by this object
  range parameters() { return boost::make_iterator_range(all); }

  /// Total number of parameter components
  std::size_t n_parameters() {
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
  void linearise() {
    // Initialise to zero Jacobian columns of intermediate and non trivial roots
    for (std::size_t j=n_independents(); j<n_parameters(); ++j) {
      jacobian_transpose.col(j).zero();
    }
    evaluator eval(unit_cell, &jacobian_transpose);
    accept(eval);
  }

  /// Apply the given shifts to the independent parameters
  /** Indexing is, as for the Jacobian, enforced by parameter::index()
   */
  void apply_shifts(af::const_ref<double> const &shifts) {
    SMTBX_ASSERT(shifts.size() == n_independents());
    BOOST_FOREACH(parameter *p, all) {
      if (!p->n_arguments() && p->is_variable()) {
        double *x = p->components();
        if (!x) throw cant_modify_parameter_value_error(p);
        std::size_t n = p->size();
        double const *s = &shifts[p->index()];
        for (std::size_t i=0; i<n; ++i) x[i] += s[i];
      }
    }
  }

  /// Store all crystallographic parameter values into their respective
  /// scatterers.
  void store() {
    BOOST_FOREACH(parameter *p, all) {
      crystallographic_parameter
      *cp = dynamic_cast<crystallographic_parameter *> (p);
      if (cp) cp->store(unit_cell);
    }
  }

  /// Let the given visitor visits all parameters.
  /** All nodes are whitened before returning, allowing another visit to proceed
      immediately.
   */
  template <class Visitor>
  void accept(Visitor &v) {
    BOOST_FOREACH(parameter *p, all) if (p->is_root()) v.visit(p);
    whiten();
  }

private:
  void whiten() {
    BOOST_FOREACH(parameter *p, all) p->set_colour(white);
  }

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
  uctbx::unit_cell const &unit_cell;
  parameter_array_t all;
  std::size_t n_independents_, n_intermediates_, n_non_trivial_roots_;
};


}}}


#endif // GUARD
