#include <scitbx/lstbx/normal_equations.h>
#include <scitbx/constants.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/simple_io.h>
#include <boost/lexical_cast.hpp>
#include <scitbx/random/mersenne_twister.h>
#include <boost/random/uniform_real.hpp>
#include <scitbx/random/variate_generator.h>
#include <tbxx/time_accu.hpp>

using namespace scitbx::lstbx::normal_equations;
namespace af = scitbx::af;
using scitbx::constants::pi;

template<typename T, template<typename> class SumOfRank1Updates>
void exercise(int m, int n) {
  typedef scitbx::boost_random::mt19937 engine_t;
  typedef boost::uniform_real<T> distribution_t;
  typedef scitbx::random::variate_generator<engine_t, distribution_t> variate_t;

  engine_t engine(0);
  distribution_t distribution(-10., 10.);
  variate_t variate(engine, distribution);

  af::shared<T> grad_yc(n);
  non_linear_ls_with_separable_scale_factor<T,
                                            SumOfRank1Updates> nls(n);
  tbxx::time_accu building, solving;
  for(int i=0; i<m; i++) {
    af::ref<T> g = grad_yc.ref();
    for(int j=0; j<n; j++) g[j] = variate();
    //SCITBX_EXAMINE(g);
    T yc = i, yo = 2*i;
    building.set_mark();
    nls.add_equation(yc, g, yo, 1.);
    building.accumulate();
  }
  building.set_mark();
  nls.finalise();
  building.accumulate();
  std::cout << "scale factor = " << nls.optimal_scale_factor() << "\n";
  std::cout << "chi^2 = " << nls.chi_sq() << "\n";
  linear_ls<T> step = nls.step_equations();
  solving.set_mark();
  step.solve();
  solving.accumulate();
  af::const_ref<T> s = step.solution().const_ref();
  std::cout << "step = [ " << s[0];
  for(int j=1; j<n; j++) {
    if(j % 50) continue;
    std::cout << " .. " << s[j];
  }
  std::cout << " ]\n";
  std::cout << "\n*** "
            << 8*sizeof(T) << "-bit floats, "
            << m << "x" << n << ", "
            << "building: " << building.as_double()
            << ", solving: " << solving.as_double() << " ***\n";
}

void help() {
  std::cerr << "bench [BLAS-2 | BLAS-3] [single | double] #data #parameters\n";
}

int main(int argc, char * argv[]) {
  using boost::lexical_cast;
  using boost::bad_lexical_cast;

  if(argc != 5) {
    help();
    return 1;
  }
  std::string lvl(argv[1]), precision(argv[2]);
  bool single_precision;
  if(precision == "single") single_precision = true;
  else if (precision == "double") single_precision = false;
  else {
    help();
    return 1;
  }
  try {
    int m = lexical_cast<int>(argv[3]);
    int n = lexical_cast<int>(argv[4]);
    if(lvl == "BLAS-2") {
      std::cout << "Level 2 BLAS implementation\n";
      std::cout << "===========================\n";
      if(single_precision)
        exercise<float, scitbx::matrix::sum_of_symmetric_rank_1_updates>(m, n);
      else
        exercise<double, scitbx::matrix::sum_of_symmetric_rank_1_updates>(m, n);
    }
    else if (lvl == "BLAS-3"){
      std::cout << "Level 3 BLAS implementation\n";
      std::cout << "===========================\n";
      if(single_precision)
        exercise<float, scitbx::matrix::rank_n_update>(m, n);
      else
        exercise<double, scitbx::matrix::rank_n_update>(m, n);
    }
    else {
      help();
      return 1;
    }
  }
  catch(const bad_lexical_cast &) {
    help();
    return 1;
  }
  return 0;
}
