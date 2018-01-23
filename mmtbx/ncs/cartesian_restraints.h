#ifndef MMTBX_NCS_CARTESIAN_RESTRAINTS_H
#define MMTBX_NCS_CARTESIAN_RESTRAINTS_H

#include <mmtbx/error.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/auto_array.h>
#include <map>
#include <vector>

namespace mmtbx { namespace ncs { namespace cartesian_restraints {

  namespace af = scitbx::af;

  struct pair_registry : boost::noncopyable
  {
    private:
      // Terminology:
      // i_seq - i_seq of atoms belong to master
      // j_seq - i_seq of atoms belong to copies
      // j_ncs - index number of NCS copy

      // 2D array for matching i_seq and j_ncs with j_seq:
      // tab_i_seqs_[i_seq][j_ncs] contains j_seq
      std::vector<scitbx::auto_array<unsigned> > tab_i_seqs_;
      // Contains  for each j_seq
      std::vector<unsigned> tab_ref_i_seqs_;
      // Array with j_ncs for each j_seq
      std::vector<unsigned> tab_i_ncs_;
      // Array for counting number of atoms in each ncs copy
      std::vector<unsigned> counts_;
    public:

    unsigned number_of_additional_isolated_sites;

    pair_registry(unsigned n_seq, unsigned n_ncs)
    :
      tab_i_seqs_(n_seq),
      tab_ref_i_seqs_(n_seq, n_seq),
      tab_i_ncs_(n_seq, n_ncs),
      counts_(n_ncs, 0),
      number_of_additional_isolated_sites(0),
      init_n_ncs(n_ncs)
    {}

    std::size_t
    n_seq() const { return tab_i_seqs_.size(); }
    unsigned init_n_ncs;

    void
    register_additional_isolated_sites(unsigned number)
    {
      number_of_additional_isolated_sites += number;
    }

    af::tiny<unsigned, 2>
    enter(unsigned i_seq, unsigned j_seq, unsigned j_ncs)
    {
      unsigned n_seq = static_cast<unsigned>(tab_i_seqs_.size());
      MMTBX_ASSERT(i_seq < n_seq);
      MMTBX_ASSERT(j_seq < n_seq);
      unsigned n_ncs = static_cast<unsigned>(counts_.size());
      MMTBX_ASSERT(j_ncs > 0);
      MMTBX_ASSERT(j_ncs < n_ncs);
      if (   tab_ref_i_seqs_[j_seq] == i_seq
          && tab_i_ncs_[j_seq] == j_ncs) {
        return af::tiny<unsigned, 2>(0u, 1u);
      }
      if (tab_i_ncs_[j_seq] != n_ncs) {
        return af::tiny<unsigned, 2>(1u, tab_i_ncs_[j_seq]);
      }
      if (tab_ref_i_seqs_[i_seq] != n_seq) {
        return af::tiny<unsigned, 2>(2u, tab_i_ncs_[i_seq]);
      }
      if (   tab_i_seqs_[i_seq].get() != 0
          && tab_i_seqs_[i_seq][j_ncs] != n_seq
          && tab_i_seqs_[i_seq][j_ncs] != j_seq) {
        return af::tiny<unsigned, 2>(3u, tab_i_seqs_[i_seq][j_ncs]);
      }
      MMTBX_ASSERT(tab_ref_i_seqs_[j_seq] == n_seq);
      tab_ref_i_seqs_[j_seq] = i_seq;
      tab_i_ncs_[j_seq] = j_ncs;
      if (!tab_i_seqs_[i_seq]) {
        tab_i_seqs_[i_seq].reset(new unsigned[n_ncs]);
        std::fill_n(tab_i_seqs_[i_seq].get(), n_ncs, n_seq);
        tab_i_seqs_[i_seq][0] = i_seq;
      }
      tab_i_seqs_[i_seq][j_ncs] = j_seq;
      counts_[j_ncs]++;
      return af::tiny<unsigned, 2>(0u, 0u);
    }

    std::auto_ptr<pair_registry>
    proxy_select(af::const_ref<std::size_t> const& iselection) const
    {
      unsigned n_seq = static_cast<unsigned>(tab_i_seqs_.size());
      unsigned n_ncs = static_cast<unsigned>(counts_.size());
      unsigned n_addl = number_of_additional_isolated_sites;
      unsigned result_n_seq = 0;
      unsigned result_n_addl = 0;
      for(std::size_t i=0;i<iselection.size();i++) {
        if (iselection[i] < n_seq) result_n_seq++;
        else                       result_n_addl++;
      }
      std::auto_ptr<pair_registry> result(
        new pair_registry(result_n_seq, counts_.size()));
      unsigned* result_counts = &*result->counts_.begin();
      result->number_of_additional_isolated_sites = result_n_addl;
      af::shared<std::size_t> reindexing_array(n_seq+n_addl, n_seq+n_addl);
      std::size_t* ra = reindexing_array.begin();
      for(std::size_t result_i_seq=0;
                      result_i_seq<iselection.size();
                      result_i_seq++) {
        MMTBX_ASSERT(iselection[result_i_seq] < n_seq+n_addl);
        ra[iselection[result_i_seq]] = result_i_seq;
        MMTBX_ASSERT(   result_i_seq < result_n_seq
                     || iselection[result_i_seq] >= n_seq);
      }
      for(unsigned i_seq=0;i_seq<n_seq;i_seq++) {
        unsigned result_i_seq = ra[i_seq];
        if (result_i_seq == n_seq+n_addl) continue;
        unsigned* i_seqs = tab_i_seqs_[i_seq].get();
        if (i_seqs != 0) {
          MMTBX_ASSERT(i_seqs[0] == i_seq);
          result->tab_i_seqs_[result_i_seq].reset(new unsigned[n_ncs]);
          unsigned* result_i_seqs = result->tab_i_seqs_[result_i_seq].get();
          std::fill_n(result_i_seqs, n_ncs, result_n_seq);
          result_i_seqs[0] = result_i_seq;
          for(unsigned j_ncs=1;j_ncs<n_ncs;j_ncs++) {
            unsigned j_seq = i_seqs[j_ncs];
            if (j_seq == n_seq) continue;
            unsigned result_j_seq = ra[j_seq];
            if (result_j_seq == n_seq+n_addl) continue;
            result_i_seqs[j_ncs] = result_j_seq;
          }
        }
        unsigned j_ncs = tab_i_ncs_[i_seq];
        if (tab_ref_i_seqs_[i_seq] == n_seq) {
          MMTBX_ASSERT(j_ncs == n_ncs);
          result->tab_ref_i_seqs_[result_i_seq] = result_n_seq;
          result->tab_i_ncs_[result_i_seq] = n_ncs;
        }
        else {
          MMTBX_ASSERT(j_ncs != n_ncs);
          unsigned result_j_seq = ra[tab_ref_i_seqs_[i_seq]];
          if (result_j_seq == n_seq+n_addl) continue;
          result->tab_ref_i_seqs_[result_i_seq] = result_j_seq;
          result->tab_i_ncs_[result_i_seq] = j_ncs;
          result_counts[j_ncs]++;
        }
      }
      return result;
    }

    std::vector<af::tiny<af::shared<std::size_t>, 2> >
    selection_pairs() const
    {
      unsigned n_ncs = static_cast<unsigned>(counts_.size());
      std::vector<af::tiny<af::shared<std::size_t>, 2> > result;
      result.reserve(n_ncs-1);
      for(unsigned j_ncs=1;j_ncs<n_ncs;j_ncs++) {
        result.resize(j_ncs);
        af::tiny<af::shared<std::size_t>, 2>& result_j = result[j_ncs-1];
        unsigned n_pairs = counts_[j_ncs];
        result_j[0].reserve(n_pairs);
        result_j[1].reserve(n_pairs);
      }
      unsigned n_seq = static_cast<unsigned>(tab_i_seqs_.size());
      for(unsigned i_seq=0;i_seq<n_seq;i_seq++) {
        unsigned* i_seqs = tab_i_seqs_[i_seq].get();
        if (i_seqs == 0) continue;
        MMTBX_ASSERT(i_seqs[0] == i_seq);
        for(unsigned j_ncs=1;j_ncs<n_ncs;j_ncs++) {
          unsigned j_seq = i_seqs[j_ncs];
          if (j_seq == n_seq) continue;
          af::tiny<af::shared<std::size_t>, 2>& result_j = result[j_ncs-1];
          result_j[0].push_back(i_seq);
          result_j[1].push_back(j_seq);
        }
      }
      for(unsigned j_ncs=1;j_ncs<n_ncs;j_ncs++) {
        af::tiny<af::shared<std::size_t>, 2>& result_j = result[j_ncs-1];
        unsigned n_pairs = counts_[j_ncs];
        MMTBX_ASSERT(result_j[0].size() == n_pairs);
        MMTBX_ASSERT(result_j[1].size() == n_pairs);
      }
      return result;
    }

    double
    adp_iso_residual_sum(
      double weight,
      double average_power,
      af::const_ref<double> const& u_isos,
      double u_average_min,
      af::ref<double> const& gradients)
    {
      MMTBX_ASSERT(   u_isos.size()
                   ==   tab_i_seqs_.size()
                      + number_of_additional_isolated_sites);
      MMTBX_ASSERT(   gradients.size() == 0
                   ||    gradients.size()
                      ==   tab_i_seqs_.size()
                         + number_of_additional_isolated_sites);
      unsigned n_seq = static_cast<unsigned>(tab_i_seqs_.size());
      unsigned n_ncs = static_cast<unsigned>(counts_.size());
      double unweighted_residual_sum = 0;
      for(unsigned i_seq=0;i_seq<n_seq;i_seq++) {
        if (!tab_i_seqs_[i_seq]) continue;
        const unsigned* i_seqs = tab_i_seqs_[i_seq].get();
        unsigned n = 0;
        double u_sum = 0;
        for(unsigned i=0;i<n_ncs;i++) {
          unsigned i_seq = i_seqs[i];
          if (i_seq == n_seq) continue;
          n++;
          u_sum += u_isos[i_seq];
        }
        double u_ave = u_sum / n;
        if (u_ave >= u_average_min) {
          double u_ave_pow = std::pow(u_ave, average_power);
          double u_diff_sum = 0;
          double u_diff_sum_sq = 0;
          for(unsigned i=0;i<n_ncs;i++) {
            unsigned i_seq = i_seqs[i];
            if (i_seq == n_seq) continue;
            double u_diff = u_ave - u_isos[i_seq];
            u_diff_sum_sq += u_diff * u_diff;
            if (gradients.size() != 0) {
              u_diff_sum += u_diff;
            }
          }
          unweighted_residual_sum += u_diff_sum_sq / u_ave_pow;
          if (gradients.size() != 0) {
            double n_ave_pow = n * u_ave_pow;
            double f1 = weight * 2 / n_ave_pow;
            double t2 = weight * average_power / (n_ave_pow * u_ave)
                      * u_diff_sum_sq;
            for(unsigned i=0;i<n_ncs;i++) {
              unsigned i_seq = i_seqs[i];
              if (i_seq == n_seq) continue;
              double u_diff = u_isos[i_seq] - u_ave;
              gradients[i_seq] += (f1 * (u_diff * n - u_diff_sum)) - t2;
            }
          }
        }
      }
      return weight * unweighted_residual_sum;
    }
  };

}}} // namespace mmtbx::ncs::cartesian_restraints

#endif // MMTBX_NCS_CARTESIAN_RESTRAINTS_H
