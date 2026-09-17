#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_H

/// Crystallographic least-squares

#include <scitbx/lstbx/normal_equations.h>
#include <scitbx/sparse/matrix.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/matrix/tensors.h>
#include <scitbx/matrix/matrix_vector_operations.h>

#include <cctbx/xray/fc_correction.h>
#include <cctbx/xray/observations.h>

#include <smtbx/error.h>
#include <smtbx/structure_factors/direct/standard_xray.h>
#include <smtbx/refinement/least_squares_twinning.h>
#include <smtbx/refinement/ml_target.h>
#include <smtbx/refinement/weighting_schemes.h>

#include <algorithm>
#include <vector>
#if defined(_OPENMP)
  #include <omp.h>
#endif
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/smart_ptr/scoped_ptr.hpp>
#include <boost/thread.hpp>

#ifdef HAVE_GCCVISIBILITYPATCH
#define DllExport __attribute__ ((visibility("default")))
#else
#ifdef _MSC_VER
#define DllExport   __declspec( dllexport )
#endif
#ifdef __BORLANDC__
#define DllExport __export
#endif
#ifdef __GNUC__
#define DllExport
#endif
#endif

// returns false to interrupt, true to continue
typedef bool (*ProgressListener)(size_t max, size_t pos);

namespace smtbx { namespace refinement { namespace least_squares {
  ProgressListener& GetRefinementProgressListener();

  namespace lstbx = scitbx::lstbx;

  /** builder_base hands the design matrix out as af::versa<FloatType>, which a
      narrowed store cannot satisfy. It does not need to: the products the
      conjugate gradients want are done inside the builder (see
      build_design_matrix::times), and the reason the matrix stays in C++ at all
      is to avoid a second copy of the whole matrix. So the accessor works when the
      store is the working type and refuses otherwise, rather than quietly
      widening the matrix back out and allocating what was just saved.
   */
  template <typename FloatType, typename StoreType>
  struct design_matrix_of_float_type {
    static af::versa<FloatType, af::c_grid<2> > const &
    get(af::versa<StoreType, af::c_grid<2> > const &) {
      SMTBX_NOT_IMPLEMENTED();
      throw 1;
    }
  };

  template <typename FloatType>
  struct design_matrix_of_float_type<FloatType, FloatType> {
    static af::versa<FloatType, af::c_grid<2> > const &
    get(af::versa<FloatType, af::c_grid<2> > const &m) { return m; }
  };

  template <typename FloatType>
  class builder_base {
  public:
    builder_base()
      : interrupted(false)
    {}
    virtual ~builder_base() {}
    virtual af::versa<FloatType, af::c_grid<2> > const& design_matrix() const = 0;
    virtual bool has_design_matrix() const = 0;
    virtual cctbx::xray::observations<FloatType> const& reflections() const = 0;
    virtual af::shared<std::complex<FloatType> > const& f_calc() const = 0;
    virtual af::shared<FloatType> const& observables() const = 0;
    virtual af::shared<FloatType> const& weights() const = 0;

    static int get_available_threads() {
      int& available = available_threads_var();
      if (available == -1) {
        available = std::max(1,
          static_cast<int>(boost::thread::physical_concurrency()));
      }
      return available;
    }

    static void set_available_threads(int thread_count) {
      // limit to the logical cores count
      available_threads_var() =
        std::max(1, std::min(
          static_cast<int>(boost::thread::hardware_concurrency()),
          thread_count));
    }

    static bool has_openmp() {
#if defined(_OPENMP)
      return true;
#endif
      return false;
    }
    bool OnProgress(size_t max, size_t pos) const {
      if (interrupted) {
        return false;
      }
      ProgressListener l = GetRefinementProgressListener();
      if (l != 0) {
        interrupted = !(*l)(max, pos);
      }
      return !interrupted;
    }
    void interrupt() {
      interrupted = true;
    }

    /** @brief How long the main thread waits on a worker before looking up.

    The workers cannot simply be joined and forgotten: a refinement has to stay
    interruptible, so the progress listener must keep being called while they
    run, and it is that listener's answer which cancels the build.

    That used to be done by polling a flag with a fixed 100 ms sleep between
    looks, which put a 100 ms floor under every threaded build. The
    accumulation itself is a few milliseconds on anything short of a protein,
    so on a small structure that wait *was* the cost of threading -- it is why a
    threaded build measured a flat tenth of a second whatever the structure, and
    why parallelising a small one came out several times slower than leaving it
    serial.

    A bounded join replaces the poll: the thread is joined the instant it
    finishes, and the wait gives up every 100 ms only so the listener can be
    called. A sleep-poll cannot do that at any interval -- asking for a shorter
    sleep does not help much, because the granularity of a Windows sleep is
    about a millisecond however small the argument, and asking for a longer one
    is the original problem.

    Joining is also what makes reading a worker's results safe: `running` is a
    plain bool, so observing it false is not on its own a guarantee that
    everything written before it is visible here.
    */
    static boost::chrono::milliseconds progress_interval() {
      return boost::chrono::milliseconds(100);
    }

    /** @brief Bounds on how many workers accumulate a normal matrix each.

    A ceiling and a floor, not a target. Every worker keeps a private copy of
    the normal matrix, so both the memory a build costs and the O(n^2) it
    spends allocating, zeroing and merging those copies grow with the thread
    count, while the arithmetic there is to share does not. Returns therefore
    fall away, and past some point another worker costs more than it brings.

    The ceiling guards a machine with a great many cores against putting a copy
    of the matrix on every one of them; on anything with fewer hardware threads
    than this it does nothing at all. The floor guards the other end, where a
    structure big enough that one matrix overruns the memory budget would
    otherwise be left with a single worker and no sharing at all -- some is
    worth having even when the budget cannot pay for it in full.
     */
    //@{
    static int max_accumulator_threads() { return 16; }
    static int min_accumulator_threads() { return 3; }
    //@}

    /** @brief A per-thread normal equations object, told its memory budget
        where the type in question has one to be told.

    Not every normal equations type the builder is instantiated on takes a
    buffer size -- only the ones whose accumulator buffers rows. The two
    overloads pick themselves: the first exists only when the three-argument
    constructor does, and the int/long argument makes it the better match when
    both are viable. Plain overload resolution, so no C++11 or later required.
     */
    //@{
    template <typename NE>
    static NE *new_chunk_equations(int n, std::size_t buffer_bytes, int,
                                   decltype(NE(0, true, std::size_t(0)))* = 0)
    {
      return new NE(n, true, buffer_bytes);
    }

    template <typename NE>
    static NE *new_chunk_equations(int n, std::size_t, long) {
      return new NE(n);
    }
    //@}

    /** @brief How much work one worker needs before it pays for itself.

    Model a threaded build as

        time ~ F/T + n^2 n_refl/(T r) + T n^2 m

    -- the structure factor pass and the arithmetic both shared over T threads,
    and against them the O(n^2) each worker pays for a private normal matrix,
    which grows with T. Minimising over T,

        T* = sqrt(n_refl / k),   k = r m

    and **n cancels**: how many workers are worth having follows the number of
    reflections, not the number of parameters. That is a property of the cost
    structure and holds anywhere; only k, a ratio of merge cost to arithmetic
    rate, belongs to the hardware.

    Two constants because the accumulators differ greatly in m. The packed one
    keeps half a matrix and folds each row in as it arrives, so it merges less
    and wants more workers for the same work; the buffered one carries a full
    matrix, a row buffer and a syrk besides. They are told apart as in
    new_chunk_equations above -- the buffered accumulator is the one that takes
    a buffer size.

    Both values are deliberately conservative. Overshooting T* costs more than
    undershooting it, every extra worker being another copy of the matrix, and
    the optimum is flat around its minimum, so erring low gives up very little.
     */
    //@{
    template <typename NE>
    static double accumulator_work_constant(
      int, decltype(NE(0, true, std::size_t(0)))* = 0) { return 195.; }

    template <typename NE>
    static double accumulator_work_constant(long) { return 30.; }

    static int threads_for_work(double work_constant, std::size_t n_reflections)
    {
      double const t = std::sqrt(double(n_reflections)/work_constant);
      return std::max(1, static_cast<int>(t + 0.5));
    }
    //@}
  private:
    mutable bool interrupted;
    static int& available_threads_var() {
      static int available = -1;
      return available;
    }
  };


  /** \brief Build normal equations for the given data, model, weighting
  and constraints. Optionally builds the design matrix.

  The constraints is performed with a reparametrisation whose Jacobian
  transpose is passed as an argument.

  What is being built, in the usual notation: a crystallographic refinement
  minimises

      L(x) = sum_h w_h ( K y_c(h; x) - y_o(h) )^2

  over the refined parameters x, with y the observable (normally |Fc|^2 against
  Fo^2), w the weights and K an overall scale. That is a non-linear problem
  solved by repeated linearisation -- each cycle linearises y_c about the
  current x, giving the design matrix J with one row per reflection,

      J_{h,k} = dy_c(h) / dx_k

  and then solves the *normal equations* of that linear problem for the shift,

      (J^T W J) dx = J^T W r,      r_h = y_o(h) - K y_c(h)

  W being the diagonal of weights. The whole cost of a refinement cycle is
  filling in J^T W J: the loop below visits each reflection once, computes its
  y_c and its row of J -- see least_squares_fc.h for the two changes of variable
  that row is dragged through -- and hands the row to an accumulator which adds
  its rank-1 contribution w_h j_h j_h^T. Summing rank-1 updates is exactly what
  makes this a BLAS-3 operation if the rows are batched (syrk over a block of
  rows) rather than a BLAS-2 one, which is what
  scitbx/matrix/symmetric_rank_1_update.h is about.

  The scale factor K is not among the x. It is eliminated analytically at each
  cycle, the problem being linear in it -- separable, or variable-projection,
  least squares; scitbx/lstbx/normal_equations.h carries that derivation.

  The accumulator is fed whether or not the design matrix is being stored, so
  build_design_matrix really does build both, as the name says. That matters
  for the stored-J conjugate gradients: they want J and, alongside it, the
  scale factor, the right hand side and the preconditioner blocks, and those
  cost a few flops per reflection while the gradients are still in cache
  against a second pass over the whole matrix if they are recovered from J
  afterwards.
  What is accumulated is the caller's choice, the accumulator being a template
  parameter -- see smtbx/refinement/least_squares_matrix_free.h for one which
  gathers exactly that summary and never forms a normal matrix.
  */
  /** StoreType is what the design matrix is *held* in, which need not be what
      the refinement computes in. Holding it in float halves a matrix that is
      large problem, and the conjugate gradients still accumulate their
      products in FloatType -- see matrix_vector_mixed in
      scitbx/matrix/matrix_vector_operations.h for what that is worth and what
      it costs. Everything else, the gradients included, stays FloatType: only
      the stored row is narrowed, and only as it is written.
   */
  template <typename FloatType,
    bool build_design_matrix,
    typename StoreType = FloatType>
    struct build_design_matrix_and_normal_equations : public builder_base<FloatType> {

    typedef builder_base<FloatType> parent_t;
    typedef StoreType store_type;

    typedef f_calc_function_base<FloatType>
      f_calc_function_base_t;

    typedef boost::shared_ptr<f_calc_function_base_t>
      one_miller_index_fcalc_ptr_t;
    typedef boost::shared_ptr<fc_correction<FloatType> >
      fc_correction_ptr_t;

    build_design_matrix_and_normal_equations(
      cctbx::xray::observations<FloatType> const& reflections,
      MaskData<FloatType> const& f_mask_data,
      boost::optional<FloatType> scale_factor,
      f_calc_function_base_t& f_calc_function,
      scitbx::sparse::matrix<FloatType> const&
        jacobian_transpose_matching_grad_fc,
      cctbx::xray::fc_correction<FloatType> const& fc_cr,
      bool objective_only = false,
      bool may_parallelise = false,
      bool use_openmp = false,
      int max_memory = 300,
      ml_data<FloatType> const& ml = ml_data<FloatType>())
      :
      reflections_(reflections),
      f_mask_data(f_mask_data),
      scale_factor(scale_factor),
      f_calc_function(f_calc_function),
      jacobian_transpose_matching_grad_fc(jacobian_transpose_matching_grad_fc),
      fc_cr(fc_cr),
      objective_only(objective_only),
      may_parallelise(may_parallelise),
      built(false),
      use_openmp(use_openmp),
      max_memory(max_memory),
      ml(ml),
      f_calc_(reflections.size()),
      observables_(reflections.size()),
      weights_(reflections.size()),
      /* Left uninitialised: the pass writes every row of every column, so
         value-initialising first is a second pass over the whole matrix for
         nothing -- and that is also a full set of first-touch page faults
         taken twice over.
       */
      design_matrix_(af::c_grid<2>(build_design_matrix ? reflections.size() : 0,
        build_design_matrix ? jacobian_transpose_matching_grad_fc.n_rows() : 0),
        af::init_functor_null<StoreType>())
    {
      /* A design matrix of gradients is exactly what an objective-only build
         declines to compute, so asking for both is a caller error rather than
         a case to support: the gradients are never gathered and the matrix
         would come out zero while has_design_matrix() reported it built.
         Both are Python arguments -- see smtbx/refinement/boost_python/
         least_squares.cpp -- so the combination has to be refused here rather
         than left to the one caller which happens not to ask for it.
       */
      SMTBX_ASSERT(!(build_design_matrix && objective_only));
      check_dispersion_correction();
      check_maximum_likelihood();
    }

    /** @brief Refuse a radial f'/f'' correction on twinned data.

    The correction accumulates its gradients during one computation of Fc, and
    twinning_processor computes Fc once per twin component, each call starting
    the accumulation over. Only the last component's contribution would survive.
    Summing them with their component scales is what it would take, in
    least_squares_twinning.h; until then this says so rather than returning a
    wrong gradient. Batch scaling (HKLF 2) is a different thing and is fine: it
    scales one already-complete gradient vector.
    */
    void check_dispersion_correction() const {
      cctbx::xray::dispersion_radial_correction<FloatType> const *dc =
        f_calc_function.get_dispersion_correction();
      SMTBX_ASSERT(!(dc && dc->grad && reflections_.is_twinned()));
    }

    /** @brief What maximum likelihood cannot be combined with, refused up front.

    The likelihood is a statement about one amplitude and its distribution, so
    the things it cannot be mixed with are the ones that change what the
    computed observable means:

    - twinning, because the components are summed as intensities and the
      amplitudes of different components do not add;
    - a non-trivial Fc correction (extinction, SWAT), because it returns an
      intensity multiplier which has no place in the amplitude distribution;
    - anything that is not one alpha, beta, epsilon and centric flag per
      reflection.

    Each of these would otherwise produce numbers rather than an error, which is
    the worse outcome. OpenMP is refused separately, in build().
    */
    void check_maximum_likelihood() const {
      if (!ml.active()) {
        return;
      }
      SMTBX_ASSERT(ml.is_consistent_with(reflections_.size()));
      SMTBX_ASSERT(!reflections_.is_twinned());
      SMTBX_ASSERT(fc_cr.is_trivial());
    }

    template<class NormalEquations>
    build_design_matrix_and_normal_equations(
      NormalEquations& normal_equations,
      cctbx::xray::observations<FloatType> const& reflections,
      MaskData<FloatType> const& f_mask_data,
      IWeightingScheme<FloatType> const& weighting_scheme,
      boost::optional<FloatType> scale_factor,
      f_calc_function_base_t &f_calc_function,
      scitbx::sparse::matrix<FloatType> const &
        jacobian_transpose_matching_grad_fc,
      cctbx::xray::fc_correction<FloatType> const& fc_cr,
      bool objective_only = false,
      bool may_parallelise = false,
      bool use_openmp = false,
      int max_memory = 300,
      ml_data<FloatType> const& ml = ml_data<FloatType>())
      :
      reflections_(reflections),
      f_mask_data(f_mask_data),
      scale_factor(scale_factor),
      f_calc_function(f_calc_function),
      jacobian_transpose_matching_grad_fc(jacobian_transpose_matching_grad_fc),
      fc_cr(fc_cr),
      objective_only(objective_only),
      may_parallelise(may_parallelise),
      built(false),
      use_openmp(use_openmp),
      max_memory(max_memory),
      ml(ml),
      f_calc_(reflections.size()),
      observables_(reflections.size()),
      weights_(reflections.size()),
      /* Left uninitialised: the pass writes every row of every column, so
         value-initialising first is a second pass over the whole matrix for
         nothing -- and that is also a full set of first-touch page faults
         taken twice over.
       */
      design_matrix_(af::c_grid<2>(build_design_matrix ? reflections.size() : 0,
        build_design_matrix ? jacobian_transpose_matching_grad_fc.n_rows() : 0),
        af::init_functor_null<StoreType>())
    {
      // as above: refused before anything is built
      SMTBX_ASSERT(!(build_design_matrix && objective_only));
      check_dispersion_correction();
      check_maximum_likelihood();
      build(normal_equations, weighting_scheme);
    }

    template<class NormalEquations>
    void build(NormalEquations& normal_equations,
      IWeightingScheme<FloatType> const& weighting_scheme)
    {
      typedef boost::shared_ptr<NormalEquations>
              normal_equations_ptr_t;
      typedef accumulate_reflection_chunk<NormalEquations>
        accumulate_reflection_chunk_t;
      typedef boost::shared_ptr<accumulate_reflection_chunk_t>
              accumulate_reflection_chunk_ptr_t;
      if (built) {
        return;
      }
      // Accumulate equations Fo(h) ~ Fc(h)
      reflections_.update_prime_fraction();
      twinning_processor<FloatType> twp(reflections_, f_mask_data, !objective_only,
        jacobian_transpose_matching_grad_fc);
      if (may_parallelise) {
        //!!
        scitbx::matrix::tensors::initialise<FloatType>();
#if defined(_OPENMP)
        /* The OpenMP accumulation calls add_equations_omp, which is specialised
           per accumulator in normal_equations_omp.h and is not provided for the
           fixed-scale accumulator maximum likelihood uses. Refusing here beats
           discovering it as a throw from inside a parallel region, and beats
           far more a specialisation that silently accumulates the wrong matrix.
         */
        SMTBX_ASSERT(!(use_openmp && ml.active()));
        if (use_openmp) {
          typedef accumulate_reflection_chunk_omp<NormalEquations>
            accumulate_reflection_chunk_omp_t;
          /**
           * @brief A pointer to the normal equations object for local refinement.
           */
          normal_equations_ptr_t local_NE(new NormalEquations(normal_equations.n_parameters()));
          accumulate_reflection_chunk_omp_t job(
            *this,
            local_NE,
            reflections_, f_mask_data, twp, weighting_scheme, scale_factor,
            one_miller_index_fcalc_ptr_t(&f_calc_function, null_deleter()),
            jacobian_transpose_matching_grad_fc,
            fc_cr, objective_only,
            f_calc_.ref(), observables_.ref(), weights_.ref(),
            design_matrix_, max_memory);
          boost::thread th(boost::ref(job));
          //job();
          while (!th.try_join_for(parent_t::progress_interval())) {
            if (!this->OnProgress(~0, ~0)) {
              this->interrupt();
              th.join();
              throw SMTBX_ERROR("external_interrupt");
            }
          }
          if (job.exception_) {
            throw* job.exception_.get();
          }
          /* Unlike the two paths below, this one does not feed the accumulator
             while building the design matrix: add_equations_omp works from a
             batch of gradient rows which that build never gathers, and it is
             specific to the packed normal matrix in any case. So there is
             nothing to merge here, and callers wanting both at once have to
             stay off use_openmp -- which smtbx/refinement/cgls.py does.
           */
          if (!build_design_matrix) {
            normal_equations = *local_NE;
            normal_equations.finalise(objective_only);
          }
          built = true;
          return;
        }
#endif
        /* Every worker keeps a normal matrix of its own, so what a threaded
           build costs in memory grows with the thread count, and nothing in
           `available_threads` knows that. On a large structure this is the
           dominant allocation and can run to many times the budget asked for.

           Bounding the workers by what the budget holds is not the trade it
           appears to be: past a handful of them the extra copies buy overhead
           and duplicate storage rather than arithmetic, the heavy kernel being
           threaded by the BLAS already. A build bounded this way is usually
           quicker as well as smaller.
         */
        int thread_count = parent_t::get_available_threads();
        thread_count = std::min(thread_count, parent_t::max_accumulator_threads());
        /* The work law is calibrated on structure factors, whose cost per
           reflection is roughly the same for all of them. It must not be
           applied to a functor whose reflections are far dearer than that: the
           dynamical electron diffraction one computes a beam group per
           reflection, orders of magnitude more work, so counting reflections
           understates the work by about as much and would leave a build that
           parallelises perfectly well with a handful of threads.

           Told apart by raw_gradients(), as the Jacobian product above is: a
           functor handing back gradients already in the basis of the refined
           parameters is the dynamical one, and it is the same distinction
           rather than a new flag invented for this.
         */
        if (f_calc_function.raw_gradients()) {
          thread_count = std::min(thread_count, parent_t::threads_for_work(
            parent_t::template accumulator_work_constant<NormalEquations>(0),
            reflections_.size()));
        }
        if (max_memory > 0 && normal_equations.n_parameters() > 0) {
          std::size_t const n = normal_equations.n_parameters();
          std::size_t const result_bytes = n*n*sizeof(FloatType);
          /* A quarter of the budget, not all of it. The private matrices are
             not what the budget was written for -- the row buffers are -- and
             giving them the whole of it lets a *generous* budget produce a
             slower and hungrier build than a modest one, by making room for
             copies that do not pay for themselves.

             Floored, and on a very large structure the floor matters more than
             the cap: once one matrix is itself larger than the share, the
             division alone would say a single worker and give up threading
             altogether, costing far more than the memory it saves.
           */
          std::size_t const budget = (std::size_t(max_memory) << 20)/4;
          int const by_memory = static_cast<int>(
            std::max<std::size_t>(parent_t::min_accumulator_threads(),
              budget/std::max<std::size_t>(1, result_bytes)));
          thread_count = std::min(thread_count, by_memory);
        }
        /* Held here rather than in a boost::thread_group so that each thread
           can be joined individually with a bounded wait; a thread_group offers
           join_all and nothing timed, and the progress listener has to keep
           being called while the wait is in progress.
         */
        std::vector<boost::shared_ptr<boost::thread> > pool;
        std::vector<accumulate_reflection_chunk_ptr_t> accumulators;
        Scheduler scheduler(reflections_.size());
        /* Each worker buffers rows of its own and, left to itself, takes the
           accumulator's default -- generous, and multiplied by the thread count
           far past any budget the caller set. Only the main accumulator ever
           honoured max_memory; these did not.

           The budget is shared out instead. Fewer rows to a chunk means more
           chunks and one more pass over the result apiece, which is the trade;
           the floor keeps that from becoming silly on a machine with many
           threads and a small budget.
         */
        std::size_t const thread_buffer_bytes = max_memory > 0
          ? std::max<std::size_t>(4u << 20,
              ((std::size_t(max_memory) << 20)/4)/thread_count)
          : 0;
        for(int thread_idx=0; thread_idx<thread_count; thread_idx++) {
          normal_equations_ptr_t chunk_normal_equations(
            parent_t::template new_chunk_equations<NormalEquations>(
              normal_equations.n_parameters(), thread_buffer_bytes, 0));
          accumulate_reflection_chunk_ptr_t accumulator(
            new accumulate_reflection_chunk_t(
              *this,
              scheduler,
              chunk_normal_equations,
              reflections_, f_mask_data, twp, weighting_scheme, scale_factor,
              one_miller_index_fcalc_ptr_t(f_calc_function.fork()),
              jacobian_transpose_matching_grad_fc,
              fc_correction_ptr_t(fc_cr.fork()),
              objective_only,
              f_calc_.ref(), observables_.ref(), weights_.ref(),
              design_matrix_, ml));
          accumulators.push_back(accumulator);
          pool.push_back(boost::shared_ptr<boost::thread>(
            new boost::thread(boost::ref(*accumulator))));
        }
        /* Joined in turn, each with a bounded wait so that the progress
           listener still runs while they work. Waiting on them one after
           another rather than all at once loses nothing: the listener is called
           at the same cadence whichever of them is outstanding, and the last
           join returns when the last thread does either way.
         */
        for (int thread_idx = 0; thread_idx < thread_count; thread_idx++) {
          while (!pool[thread_idx]->try_join_for(
                   parent_t::progress_interval()))
          {
            if (!this->OnProgress(~0, ~0)) {
              this->interrupt();
              for (int j = thread_idx; j < thread_count; j++) {
                pool[j]->join();
              }
              throw SMTBX_ERROR("external_interrupt");
            }
          }
        }

        for (int thread_idx = 0; thread_idx < thread_count; thread_idx++) {
          // note this check has to run whether or not there is anything to
          // merge: it is the only place a worker thread's exception surfaces
          if (accumulators[thread_idx]->exception_) {
            throw* accumulators[thread_idx]->exception_.get();
          }
          normal_equations += accumulators[thread_idx]->normal_equations;
        }
        normal_equations.finalise(objective_only);
      }
      else {
        Scheduler scheduler(reflections_.size());
        accumulate_reflection_chunk_t job(
          *this,
          scheduler,
          normal_equations_ptr_t(&normal_equations, null_deleter()),
          reflections_, f_mask_data, twp, weighting_scheme, scale_factor,
          one_miller_index_fcalc_ptr_t(f_calc_function.fork()),
          jacobian_transpose_matching_grad_fc,
          fc_correction_ptr_t(fc_cr.fork()),
          objective_only,
          f_calc_.ref(), observables_.ref(), weights_.ref(),
          design_matrix_, ml);
        job();
        if (job.exception_) {
          throw *job.exception_.get();
        }
        normal_equations.finalise(objective_only);
      }
      built = true;
    }

    virtual cctbx::xray::observations<FloatType> const& reflections() const {
      return reflections_;
    }

    virtual af::shared<std::complex<FloatType> > const& f_calc() const { return f_calc_; }

    virtual af::shared<FloatType> const& observables() const { return observables_; }

    virtual af::shared<FloatType> const& weights() const { return weights_; }

    virtual bool has_design_matrix() const {
      return build_design_matrix && built;
    }

  protected:
    af::versa<FloatType, af::c_grid<2> > const& design_matrix() const {
      return design_matrix_of_float_type<FloatType, StoreType>::get(design_matrix_);
    }
    /// the store itself, whatever it is held in
    af::versa<StoreType, af::c_grid<2> > const& design_matrix_store() const {
      return design_matrix_;
    }
#if defined(_OPENMP)
    #include "least_squares_omp.h"
#endif
    struct chunk {
      const int idx, size;
      chunk()
        : idx(0), size(0)
      {}
      chunk(int idx, int size)
        : idx(idx), size(size)
      {}
    };
    struct Scheduler {
      int count, current;
      boost::mutex mtx;
      Scheduler(int count)
        : count(count),
        current(0)
      {}
      chunk next() {
        boost::mutex::scoped_lock lock(mtx);
        int left = count - current;
        if (left == 0) {
          return chunk();
        }
        int sz = std::min(left, 256);
        chunk rv(current, sz);
        current += sz;
        return rv;
      }
      void reset() {
        current = 0;
      }
    };

    /// the shared one, defined next to f_calc_function_base which uses it
    typedef least_squares::flattened_jacobian_transpose<FloatType>
      flattened_jacobian_transpose;

    /** @brief The gradients of the observable with respect to the parameters,
        written into somewhere the caller already has.

    The functor gives its gradients in the basis of the structure factor's own
    parameters, which the reparametrisation's Jacobian transpose then maps onto
    the refined ones -- or gives them in that basis already, if it has applied
    the Jacobian itself.

    The destination is a bare pointer rather than a vector so that the
    accumulator's own row may be passed, which is where they would have to be
    copied to otherwise. The Jacobian is a pointer too, and null when the
    functor's gradients need none applying -- there is then nothing to flatten
    and nothing built.
    */
    static void fill_gradients(f_calc_function_base_t const &f_calc_function,
      flattened_jacobian_transpose const *jacobian_transpose,
      FloatType *destination, int n_params)
    {
      if (f_calc_function.raw_gradients()) {
        SMTBX_ASSERT(jacobian_transpose != 0);
        /* Fused when the functor can: the gradients of the observable are
           converted from the complex ones a component at a time, as the
           Jacobian reaches each column, and never assembled into a vector of
           their own. Whoever cannot do that says so and is served the long way.
         */
        if (f_calc_function.apply_jacobian_to_grad_observable(
              *jacobian_transpose, destination))
        {
          return;
        }
        jacobian_transpose->apply(f_calc_function.get_grad_observable(),
          destination);
        return;
      }
      // already in the basis of the refined parameters, so nothing to apply
      af::const_ref<FloatType> g = f_calc_function.get_grad_observable();
      SMTBX_ASSERT(g.size() == n_params)(g.size())(n_params);
      std::copy(g.begin(), g.end(), destination);
    }

    /** @brief One row of sqrt(W) J.

    Stored already multiplied by sqrt(w), so that the rows are those of
    sqrt(W) J and a matrix-vector product against them is a plain gemv. Scaling
    afterwards instead is a read, a multiply and a write over the whole matrix
    for nothing, the store being bandwidth bound either way. The accumulator
    sees the unscaled gradients, which is what it wants: it applies w itself.
    */
    static void store_design_row(
      af::versa<StoreType, af::c_grid<2> > &design_matrix,
      int i_h, FloatType const *gradients, int n_params, FloatType weight)
    {
      FloatType const sqrt_weight = std::sqrt(weight);
      StoreType *row = &design_matrix(i_h, 0);
      for (int gi = 0; gi < n_params; gi++) {
        // narrowed here and only here, once the row is final
        row[gi] = static_cast<StoreType>(sqrt_weight*gradients[gi]);
      }
    }

    /// Accumulate from reflections whose indices are
    /// returned by scheduler
    template<class NormalEquations>
    struct accumulate_reflection_chunk {
      builder_base<FloatType>& parent;
      Scheduler& scheduler;
      boost::scoped_ptr<smtbx::error> exception_;
      boost::shared_ptr<NormalEquations> normal_equations_ptr;
      NormalEquations &normal_equations;
      cctbx::xray::observations<FloatType> const &reflections;
      MaskData<FloatType> const& f_mask_data;
      twinning_processor<FloatType> const& twp;
      IWeightingScheme<FloatType> const &weighting_scheme;
      boost::optional<FloatType> scale_factor;
      boost::shared_ptr<f_calc_function_base_t> f_calc_function_ptr;
      f_calc_function_base_t &f_calc_function;
      scitbx::sparse::matrix<FloatType> const
        &jacobian_transpose_matching_grad_fc;
      boost::shared_ptr<cctbx::xray::fc_correction<FloatType> > fc_cr;
      bool objective_only, compute_grad;
      af::ref<std::complex<FloatType> > f_calc;
      af::ref<FloatType> observables;
      af::ref<FloatType> weights;
      af::versa<StoreType, af::c_grid<2> > &design_matrix;
      /* Held by reference, not copied per worker. It is read-only and indexed
         by i_h, so unlike f_calc_function and fc_cr there is nothing in it to
         fork - and a copy would refcount its af::shared members, which is not
         atomic, so workers destroying their copies would race. The builder owns
         it and outlives every worker.
       */
      ml_data<FloatType> const& ml;
      bool running;
      accumulate_reflection_chunk(
        builder_base<FloatType>& parent,
        Scheduler& scheduler,
        boost::shared_ptr<NormalEquations> const& normal_equations_ptr,
        cctbx::xray::observations<FloatType> const &reflections,
        MaskData<FloatType> const& f_mask_data,
        twinning_processor<FloatType> const& twp,
        IWeightingScheme<FloatType> const &weighting_scheme,
        boost::optional<FloatType> scale_factor,
        boost::shared_ptr<f_calc_function_base_t> const &f_calc_function_ptr,
        scitbx::sparse::matrix<FloatType> const
          &jacobian_transpose_matching_grad_fc,
        boost::shared_ptr<cctbx::xray::fc_correction<FloatType> > const &fc_cr,
        bool objective_only,
        af::ref<std::complex<FloatType> > f_calc,
        af::ref<FloatType> observables,
        af::ref<FloatType> weights,
        af::versa<StoreType, af::c_grid<2> > &design_matrix,
        ml_data<FloatType> const& ml)
      : parent(parent),
        scheduler(scheduler),
        normal_equations_ptr(normal_equations_ptr), normal_equations(*normal_equations_ptr),
        reflections(reflections), f_mask_data(f_mask_data), twp(twp),
        weighting_scheme(weighting_scheme),
        scale_factor(scale_factor),
        f_calc_function_ptr(f_calc_function_ptr), f_calc_function(*f_calc_function_ptr),
        jacobian_transpose_matching_grad_fc(jacobian_transpose_matching_grad_fc),
        fc_cr(fc_cr),
        objective_only(objective_only), compute_grad(!objective_only),
        f_calc(f_calc), observables(observables), weights(weights),
        design_matrix(design_matrix),
        ml(ml),
        running(true)
      {}

      void operator()() {
        running = true;
        try {
          af::shared<FloatType> gradients;
          int n_params = jacobian_transpose_matching_grad_fc.n_rows();
          if (compute_grad) {
            gradients.resize(n_params);
          }
          /* Flattened per worker rather than shared between them. It costs one
             pass over the non-zeros, which is nothing against the reflections a
             worker then applies it to, and it saves passing another reference
             down through every constructor. Not built at all for a pass which
             wants no derivatives, since then it is never applied to anything --
             nor for a functor which hands back gradients that are already in
             the basis of the refined parameters, the dynamical electron
             diffraction one being the case in point, which likewise never
             applies it.
           */
          boost::scoped_ptr<flattened_jacobian_transpose> jt;
          if (compute_grad && f_calc_function.raw_gradients()) {
            jt.reset(new flattened_jacobian_transpose(
              jacobian_transpose_matching_grad_fc));
          }
          /** @brief Whether nothing stands between the Jacobian product and the
              normal equations.

          The twinning processor, the Fc correction and the radial f'/f''
          correction each get a say in the gradient vector, and in an ordinary
          refinement none of the three has anything to say -- yet each is asked,
          through a virtual call, once per reflection. All three are properties
          of the build and not of the reflection, so they are settled here, once,
          and the loop below takes a path with none of them in it.

          What that path buys is not only the calls: with nothing left to modify
          the gradients after they are computed, they can be written straight
          into the row the accumulator was going to copy them into anyway.
          */
          bool const fast = compute_grad
            && twp.is_trivial()
            && fc_cr->is_trivial()
            && f_calc_function.get_dispersion_correction() == 0;
          /* The fused conversion is only reachable from the fast path, so the
             functor may stop assembling the gradients of the observable
             exactly when that path is taken -- and only when the Jacobian is
             applied at all, the other case reading them straight out.
           */
          f_calc_function.set_defer_grad_observable(
            fast && f_calc_function.raw_gradients());
          while (true) {
            chunk ch = scheduler.next();
            if (ch.size == 0) {
              break;
            }
            for (int i = 0; i < ch.size; i++) {
              int i_h = ch.idx + i;
              if (!parent.OnProgress(scheduler.count, i_h)) {
                return;
              }
              miller::index<> const& h = reflections.index(i_h);
              const twin_fraction<FloatType>* fraction = reflections.fraction(i_h);
              if (f_mask_data.size()) {
                f_calc_function.compute(h, f_mask_data.find(h), fraction, compute_grad);
              }
              else {
                f_calc_function.compute(h, boost::none, fraction, compute_grad);
              }
              f_calc[i_h] = f_calc_function.get_f_calc();
              /* The gradients live either in the accumulator's own row, on the
                 fast path, or in the local vector when something downstream is
                 going to modify them. Either way this points at them.
               */
              FloatType *grad = 0;
              if (fast) {
                grad = normal_equations.open_equation();
                fill_gradients(f_calc_function, jt.get(), grad, n_params);
              }
              else if (compute_grad) {
                grad = gradients.begin();
                fill_gradients(f_calc_function, jt.get(), grad, n_params);
                /* The radial correction of f' and f'' keeps its gradients to
                   itself -- they are not per-scatterer, so the Jacobian above
                   has nothing to say about them and leaves their slots at zero.
                   Assigning them here, rather than beside the Fc correction
                   below, puts them in before twinning_processor scales the
                   whole vector by a batch scale factor, which is what makes
                   HKLF 2 data work without any further thought.
                 */
                write_dispersion_gradients(f_calc_function, gradients);
              }
              FloatType observable;
              if (fast) {
                observable = f_calc_function.get_observable();
              }
              else {
                // sort out twinning
                observable = twp.process(i_h, f_calc_function, gradients);
                // Fc correction
                FloatType fc_k = fc_cr->compute(
                  reflections.has_wavelengths() ? reflections.wavelength(i_h) : 0,
                  h, observable, compute_grad);
                if (fc_k != 1) {
                  observable *= fc_k;
                  f_calc[i_h] *= std::sqrt(fc_k);
                }
              }
              observables[i_h] = observable;

              FloatType weight, y_obs = reflections.fo_sq(i_h);
              if (ml.active()) {
                /* Maximum likelihood: the accumulator receives an effective
                   observation and weight in place of the measured intensity
                   and the weighting scheme, so that the equations it forms are
                   those of the likelihood. See smtbx/refinement/ml_target.h.

                   A reflection the target cannot use is given zero weight
                   rather than skipped, keeping observables and weights at one
                   entry per reflection for everything downstream.
                 */
                FloatType yo_eff, w_eff;
                if (ml.is_free(i_h)) {
                  /* Held out of the target sum: alpha and beta are estimated
                     on the free set, so including it here would make that
                     estimate self-referential. Zero weight and zero residual.
                   */
                  yo_eff = observable;
                  w_eff = 0;
                }
                else if (ml.is_intensity()) {
                  // the intensity target needs the measured sigma, which is
                  // the whole of what distinguishes it from the amplitude one
                  mli_effective_observation(
                    reflections.fo_sq(i_h), reflections.sig(i_h), observable,
                    ml.alpha[i_h], ml.beta[i_h],
                    ml.scale_factor(),
                    ml.epsilon[i_h], ml.centric[i_h],
                    yo_eff, w_eff);
                }
                else {
                  /* ml.scale_factor(), not the builder's: alpha carries the
                     scale it was estimated against, and the builder's is the
                     crystallographic one R1 and the difference map read.
                   */
                  ml_effective_observation(
                    reflections.fo_sq(i_h), reflections.sig(i_h), observable,
                    ml.alpha[i_h], ml.beta[i_h],
                    ml.scale_factor(),
                    ml.epsilon[i_h], ml.centric[i_h],
                    yo_eff, w_eff);
                }
                y_obs = yo_eff;
                weight = w_eff;
              }
              else {
                weight = weighting_scheme(reflections.fo_sq(i_h),
                  reflections.sig(i_h), observable, h, scale_factor);
              }
              weights[i_h] = weight;
              if (objective_only) {
                normal_equations.add_residual(observable, y_obs, weight);
              }
              else if (fast) {
                if (build_design_matrix) {
                  store_design_row(design_matrix, i_h, grad, n_params, weight);
                }
                /* After the design matrix and not before: committing the
                   equation weights the row in place, and the row is the one
                   the gradients were written into.
                 */
                normal_equations.commit_equation(observable, grad,
                  y_obs, weight);
              }
              else {
                if (fc_cr->grad) {
                  int grad_idx = fc_cr->get_grad_index();
                  af::const_ref<FloatType> fc_cr_grads = fc_cr->get_gradients();
                  SMTBX_ASSERT(grad_idx >= 0 &&
                    grad_idx+fc_cr_grads.size() <= gradients.size());
                  FloatType grad_m = fc_cr->get_grad_Fc_multiplier();
                  if (grad_m != 1) {
                    for (int gi = 0; gi < gradients.size(); gi++) {
                     gradients[gi] *= grad_m;
                    }
                  }
                  for (int gi = 0; gi < fc_cr_grads.size(); gi++) {
                    gradients[grad_idx + gi] = fc_cr_grads[gi];
                  }
                }
                normal_equations.add_equation(observable,
                  gradients.ref(), y_obs, weight);
                if (build_design_matrix) {
                  store_design_row(design_matrix, i_h, grad, n_params, weight);
                }
              }
            }
          }
        }
        catch (smtbx::error const &e) {
          exception_.reset(new smtbx::error(e));
        }
        catch (std::exception const &e) {
          exception_.reset(new smtbx::error(e.what()));
        }
        running = false;
      }
    };


  private:
    struct null_deleter {
      void operator()(void const *) const {}
    };

  protected:
    cctbx::xray::observations<FloatType> const& reflections_;
    MaskData<FloatType> const& f_mask_data;
    boost::optional<FloatType> scale_factor;
    f_calc_function_base_t& f_calc_function;
    scitbx::sparse::matrix<FloatType> const&
      jacobian_transpose_matching_grad_fc;
    cctbx::xray::fc_correction<FloatType> const& fc_cr;
    bool objective_only,
      may_parallelise,
      use_openmp,
      built;
    int max_memory;
    /* Empty unless this is a maximum-likelihood build. Held by value, but the
       arrays inside it are references owned by the caller for the duration of
       the build - the same lifetime `reflections_` has.
     */
    ml_data<FloatType> ml;

    af::shared<std::complex<FloatType> > f_calc_;
    af::shared<FloatType> observables_;
    af::shared<FloatType> weights_;
    af::versa<StoreType, af::c_grid<2> > design_matrix_;
  };

  /** \brief Build normal equations for the given data, model, weighting
  and constraints.

  The constraints is performed with a reparametrisation whose Jacobian
  transpose is passed as an argument.
  */
  template <typename FloatType>
  struct build_normal_equations
    : public build_design_matrix_and_normal_equations<FloatType, false>
  {
    typedef build_design_matrix_and_normal_equations<FloatType, false> parent_t;
    build_normal_equations(
      cctbx::xray::observations<FloatType> const& reflections,
      MaskData<FloatType> const& f_mask_data,
      boost::optional<FloatType> scale_factor,
      f_calc_function_base<FloatType>& f_calc_function,
      scitbx::sparse::matrix<FloatType> const
      & jacobian_transpose_matching_grad_fc,
      cctbx::xray::fc_correction<FloatType> const& fc_cr,
      bool objective_only = false,
      bool may_parallelise = false,
      bool use_openmp = false)
      : parent_t(
        reflections, f_mask_data, scale_factor, f_calc_function,
        jacobian_transpose_matching_grad_fc, fc_cr,
        objective_only, may_parallelise, use_openmp)
    {}

    template<class NormalEquations>
    build_normal_equations(
       NormalEquations &normal_equations,
       cctbx::xray::observations<FloatType> const &reflections,
       MaskData<FloatType> const& f_mask_data,
       IWeightingScheme<FloatType> const &weighting_scheme,
       boost::optional<FloatType> scale_factor,
       f_calc_function_base<FloatType> &f_calc_function,
       scitbx::sparse::matrix<FloatType> const
       &jacobian_transpose_matching_grad_fc,
       cctbx::xray::fc_correction<FloatType> const &fc_cr,
       bool objective_only = false,
       bool may_parallelise = false,
       bool use_openmp = false,
       int max_memory = 300,
       ml_data<FloatType> const& ml = ml_data<FloatType>())
       : parent_t(
        normal_equations,
        reflections, f_mask_data, weighting_scheme, scale_factor, f_calc_function,
        jacobian_transpose_matching_grad_fc, fc_cr,
        objective_only, may_parallelise, use_openmp, max_memory, ml)
    {}
     virtual af::versa<FloatType, af::c_grid<2> > const& design_matrix() const {
       SMTBX_NOT_IMPLEMENTED();
       return parent_t::design_matrix_;
     }
  };

  /** \brief Build only thed esign matrix for the given data, model, weighting
  and constraints and the buld the design matrix

  The constraints is performed with a reparametrisation whose Jacobian
  transpose is passed as an argument.
  */

  template <typename FloatType, typename StoreType = FloatType>
  struct build_design_matrix
    : public build_design_matrix_and_normal_equations<FloatType, true, StoreType>
  {
    typedef build_design_matrix_and_normal_equations<FloatType, true, StoreType>
      parent_t;
    build_design_matrix(
      cctbx::xray::observations<FloatType> const& reflections,
      MaskData<FloatType> const& f_mask_data,
      boost::optional<FloatType> scale_factor,
      f_calc_function_base<FloatType>& f_calc_function,
      scitbx::sparse::matrix<FloatType> const
      & jacobian_transpose_matching_grad_fc,
      cctbx::xray::fc_correction<FloatType> const& fc_cr,
      bool objective_only = false,
      bool may_parallelise = false,
      bool use_openmp = false)
      : parent_t(
        reflections, f_mask_data, scale_factor, f_calc_function,
        jacobian_transpose_matching_grad_fc, fc_cr,
        objective_only, may_parallelise, use_openmp)
    {}

    template<class NormalEquations>
    build_design_matrix(
       NormalEquations &normal_equations,
       cctbx::xray::observations<FloatType> const &reflections,
      MaskData<FloatType> const& f_mask_data,
       IWeightingScheme<FloatType> const &weighting_scheme,
       boost::optional<FloatType> scale_factor,
       f_calc_function_base<FloatType> &f_calc_function,
       scitbx::sparse::matrix<FloatType> const
       &jacobian_transpose_matching_grad_fc,
       cctbx::xray::fc_correction<FloatType> const &fc_cr,
       bool objective_only = false,
       bool may_parallelise = false,
       bool use_openmp = false,
       int max_memory = 300,
       ml_data<FloatType> const& ml = ml_data<FloatType>())
       : parent_t(
        normal_equations,
        reflections, f_mask_data, weighting_scheme, scale_factor, f_calc_function,
        jacobian_transpose_matching_grad_fc, fc_cr,
        objective_only, may_parallelise, use_openmp, max_memory, ml)
    {}

    virtual af::versa<FloatType, af::c_grid<2> > const& design_matrix() const {
      return design_matrix_of_float_type<FloatType, StoreType>::get(
        parent_t::design_matrix_);
    }

    /// Whether the store is narrower than the arithmetic, i.e. mixed precision
    static bool is_mixed_precision() {
      return sizeof(StoreType) != sizeof(FloatType);
    }

    /// @name The products the conjugate gradients take against the design
    /// matrix, without a copy of it having to leave here
    /** The rows carry sqrt(w), so these are the products against sqrt(W) J.
        Doing them here rather than handing the matrix out avoids copying it, and
        avoids the second copy of it existing at once.

        When the store is narrower than FloatType these go through
        matrix_*_mixed, which accumulates in FloatType regardless, so halving the
        matrix costs almost nothing in accuracy. No BLAS offers that combination,
        hence the hand-written kernels; when the store *is* FloatType the BLAS
        path is used as before.
     */
    //@{
    af::shared<FloatType> times(af::const_ref<FloatType> const &x) const {
      af::versa<StoreType, af::c_grid<2> > const &a = parent_t::design_matrix_;
      int const m = static_cast<int>(a.accessor()[0]),
                n = static_cast<int>(a.accessor()[1]);
      SMTBX_ASSERT(x.size() == static_cast<std::size_t>(n))(x.size())(n);
      af::shared<FloatType> y(m, af::init_functor_null<FloatType>());
      times_(m, n, a.begin(), x.begin(), y.begin());
      return y;
    }

    af::shared<FloatType>
    transpose_times(af::const_ref<FloatType> const &x) const {
      af::versa<StoreType, af::c_grid<2> > const &a = parent_t::design_matrix_;
      int const m = static_cast<int>(a.accessor()[0]),
                n = static_cast<int>(a.accessor()[1]);
      SMTBX_ASSERT(x.size() == static_cast<std::size_t>(m))(x.size())(m);
      af::shared<FloatType> y(n, af::init_functor_null<FloatType>());
      transpose_times_(m, n, a.begin(), x.begin(), y.begin());
      return y;
    }
    //@}

  private:
    /* Overloaded rather than branched on is_mixed_precision(), because the
       BLAS calls only typecheck when the store and the vectors agree.
     */
    static void times_(int m, int n, FloatType const *a,
                       FloatType const *x, FloatType *y) {
      scitbx::matrix::matrix_vector(m, n, a, x, y);
    }
    template <typename S>
    static void times_(int m, int n, S const *a,
                       FloatType const *x, FloatType *y) {
      scitbx::matrix::matrix_vector_mixed(m, n, a, x, y);
    }
    static void transpose_times_(int m, int n, FloatType const *a,
                                 FloatType const *x, FloatType *y) {
      scitbx::matrix::matrix_transposed_vector(m, n, a, x, y);
    }
    template <typename S>
    static void transpose_times_(int m, int n, S const *a,
                                 FloatType const *x, FloatType *y) {
      scitbx::matrix::matrix_transposed_vector_mixed(m, n, a, x, y);
    }
  };

}}}

extern "C" DllExport void SetRefinementProgressListener(ProgressListener listener);

#endif // GUARD
