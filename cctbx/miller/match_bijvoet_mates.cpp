#include <cctbx/miller/match_bijvoet_mates.h>
#include <cctbx/error.h>
#include <algorithm>
#include <vector>

namespace cctbx { namespace miller {

  namespace {

    struct keyed_index
    {
      index<> h;
      std::size_t i;
    };

    struct keyed_index_less
    {
      bool operator()(keyed_index const& a, keyed_index const& b) const
      {
        return fast_less_than<>()(a.h, b.h);
      }
    };

    bool
    same_index(index<> const& a, index<> const& b)
    {
      return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
    }

    /* For every position i the position of the last occurrence of
       -miller_indices[i], or npos if there is none.

       The indices are sorted once in the lexicographic order of
       fast_less_than. That order reverses under negation, so the negated
       indices in ascending order are the sorted array read backwards, and
       one linear merge of the array with itself read backwards finds every
       h, -h match. The last occurrence wins for duplicate indices, as it
       did with the std::map this replaces (its insertions overwrote).
     */
    std::vector<std::size_t>
    mates_of_negated_indices(
      af::const_ref<index<> > const& miller_indices,
      bool assert_is_unique_set_under_symmetry)
    {
      const std::size_t npos = static_cast<std::size_t>(-1);
      std::size_t n = miller_indices.size();
      std::vector<keyed_index> sorted(n);
      for(std::size_t i=0;i<n;i++) {
        sorted[i].h = miller_indices[i];
        sorted[i].i = i;
      }
      std::sort(sorted.begin(), sorted.end(), keyed_index_less());
      if (assert_is_unique_set_under_symmetry) {
        for(std::size_t p=1;p<n;p++) {
          if (same_index(sorted[p].h, sorted[p-1].h)) {
            throw CCTBX_ERROR("miller array is not a unique set under symmetry");
          }
        }
      }
      std::vector<std::size_t> mate(n, npos);
      fast_less_than<> less;
      std::size_t p = 0;
      std::size_t q = n;
      while (p < n && q > 0) {
        index<> hp = sorted[p].h;
        index<> hq = -sorted[q-1].h;
        if (less(hp, hq)) {
          p++;
        }
        else if (less(hq, hp)) {
          q--;
        }
        else {
          // the run of hp at [p, p_end) mates with the run of -hp at
          // [q_begin, q)
          std::size_t p_end = p + 1;
          while (p_end < n && same_index(sorted[p_end].h, hp)) p_end++;
          std::size_t q_begin = q - 1;
          while (q_begin > 0
                 && same_index(sorted[q_begin-1].h, sorted[q-1].h)) q_begin--;
          std::size_t last = sorted[q_begin].i;
          for(std::size_t t=q_begin+1;t<q;t++) {
            if (sorted[t].i > last) last = sorted[t].i;
          }
          for(std::size_t t=p;t<p_end;t++) mate[sorted[t].i] = last;
          p = p_end;
          q = q_begin;
        }
      }
      return mate;
    }

  } // namespace anonymous

  void
  match_bijvoet_mates::match_(sgtbx::reciprocal_space::asu const& asu,
                              bool assert_is_unique_set_under_symmetry)
  {
    const std::size_t npos = static_cast<std::size_t>(-1);
    std::vector<std::size_t> mate = mates_of_negated_indices(
      miller_indices_.const_ref(), assert_is_unique_set_under_symmetry);
    std::vector<bool> paired_already(miller_indices_.size(), false);
    for(std::size_t i=0;i<miller_indices_.size();i++) {
      if (paired_already[i]) continue;
      if (miller_indices_[i].is_zero()) {
        singles_[0].push_back(i);
      }
      else {
        int asu_which = asu.which(miller_indices_[i]);
        CCTBX_ASSERT(asu_which != 0);
        std::size_t j = mate[i];
        if (j == npos) {
          if (asu_which > 0) {
            singles_[0].push_back(i);
          }
          else {
            singles_[1].push_back(i);
          }
        }
        else {
          if (asu_which > 0) {
            pairs_.push_back(af::tiny<std::size_t, 2>(i, j));
          }
          else {
            pairs_.push_back(af::tiny<std::size_t, 2>(j, i));
          }
          paired_already[j] = true;
        }
      }
    }
  }

  af::shared<std::size_t>
  match_bijvoet_mates::singles(char plus_or_minus) const
  {
    std::size_t j = plus_or_minus_index_(plus_or_minus);
    return singles_[j];
  }

  void
  match_bijvoet_mates::size_assert_intrinsic() const
  {
    CCTBX_ASSERT(miller_indices_.size() == size_processed());
  }

  void
  match_bijvoet_mates::size_assert(std::size_t sz) const
  {
    size_assert_intrinsic();
    CCTBX_ASSERT(sz == size_processed());
  }

  af::shared<std::size_t>
  match_bijvoet_mates::pairs_hemisphere_selection(char plus_or_minus) const
  {
    std::size_t j = plus_or_minus_index_(plus_or_minus);
    af::const_ref<pair_type> pairs_ref = pairs_.const_ref();
    af::shared<std::size_t> result((af::reserve(pairs_ref.size())));
    for(std::size_t i=0;i<pairs_ref.size();i++) {
      result.push_back(pairs_[i][j]);
    }
    return result;
  }

  af::shared<index<> >
  match_bijvoet_mates::miller_indices_in_hemisphere(char plus_or_minus) const
  {
    std::size_t j = plus_or_minus_index_(plus_or_minus);
    af::shared<index<> > result((af::reserve(pairs_.size())));
    for(std::size_t i=0;i<pairs_.size();i++) {
      result.push_back(miller_indices_[pairs_[i][j]]);
    }
    return result;
  }

  std::size_t
  match_bijvoet_mates::plus_or_minus_index_(char plus_or_minus) const
  {
    CCTBX_ASSERT(plus_or_minus == '+' || plus_or_minus == '-');
    size_assert_intrinsic();
    if (plus_or_minus == '-') return 1;
    return 0;
  }

}} // namespace cctbx::miller
