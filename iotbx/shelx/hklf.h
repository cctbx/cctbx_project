#ifndef IOTBX_SHELX_HKLF_SPIRIT_H
#define IOTBX_SHELX_HKLF_SPIRIT_H

#include <scitbx/array_family/shared.h>
#include <cctbx/miller.h>
#ifdef IOTBX_SHELX_HKLF_SPIRIT_DEBUG
  #define BOOST_SPIRIT_DEBUG
#endif
#include <iotbx/boost_spirit/fortran_numerics.h>
#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_if.hpp>
#include <boost/spirit/include/classic_assign_actor.hpp>
#include <boost/spirit/include/classic_push_back_actor.hpp>
#include <boost/spirit/include/phoenix1_primitives.hpp>
#include <boost/spirit/include/phoenix1_binders.hpp>
#include <boost/spirit/include/phoenix1_operators.hpp>

namespace iotbx { namespace shelx {

class hklf_reader
{
  public:
    typedef char char_t;
    typedef cctbx::miller::index<> miller_t;

    template <class IteratorType>
    hklf_reader(IteratorType const &first, IteratorType const &last,
                bool strict=true) {
      init(first, last, strict);
    }

    hklf_reader(std::string const &content, bool strict) {
      init(content.begin(), content.end(), strict);
    }

    template <class IteratorType>
    void init(IteratorType const &first, IteratorType const &last, bool strict)
    {
      using namespace iotbx::boost_spirit;
      typedef scanner<IteratorType> scanner_t;
      typedef rule<scanner_t> rule_t;
      typedef typename miller_t::value_type index_t;

      fortran_int_parser<int, /*Width=*/4, /*StrictWidth=*/true> i4_p;
      fortran_real_parser<
        double, /*Width=*/8, /*FracDigits=*/2, /*StrictWidth=*/true> f8_2_p;

      miller_t h;
      subrule<0> hkl;
      subrule<1> datum;
      subrule<2> sigma;
      subrule<3> extra;
      subrule<4> perhaps_more;
      subrule<5> strict_line;
      subrule<6> line;
      subrule<7> content;
      rule_t start = (
        content = +line,
        line = strict_line >> perhaps_more,
        strict_line = hkl >> datum >> sigma >> !extra,
        hkl = (    i4_p[ assign_a(h[0]) ]
                >> i4_p[ assign_a(h[1]) ]
                >> i4_p[ assign_a(h[2]) ]
                >> eps_p(!phoenix::bind(&miller_t::is_zero)(phoenix::var(h)))
              )[ push_back_a(indices_, h) ],
        datum  = f8_2_p[ push_back_a(data_) ],
        sigma  = f8_2_p[ push_back_a(sigmas_) ],
        extra  = i4_p[ push_back_a(extra_) ],
        perhaps_more = if_p (phoenix::val(strict)) [ eol_p ]
                      .else_p [ *(anychar_p - eol_p) >> eol_p ]
      );
      #ifdef IOTBX_SHELX_HKLF_SPIRIT_DEBUG
        BOOST_SPIRIT_DEBUG_RULE(start);
      #endif

      parse_info<IteratorType> info = parse(first, last, start);
      std::runtime_error not_hklf("Not a SHELX hklf file.");
      std::runtime_error empty_hklf("No data in SHELX hklf file.");
      if (!info.full && !h.is_zero()) throw not_hklf;
      if (indices_.size() == 0) throw empty_hklf;
    }

    scitbx::af::shared<miller_t> indices() { return indices_; };

    scitbx::af::shared<double> data() { return data_; }

    scitbx::af::shared<double> sigmas() { return sigmas_; }

    scitbx::af::shared<int> alphas() { return extra_; }

    scitbx::af::shared<int> batch_numbers() { return extra_; }

  private:
    scitbx::af::shared<miller_t> indices_;
    scitbx::af::shared<double> data_, sigmas_;
    scitbx::af::shared<int> extra_; // batch numbers or phases
};

}} // iotbx::shelx

#endif // GUARD
