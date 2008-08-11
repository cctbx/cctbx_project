#include <iotbx/pdb/input.h>
#include <scitbx/array_family/sort.h>
#include <boost/scoped_array.hpp>

namespace iotbx { namespace pdb {

  namespace {

    void
    append_residue_group(
      const detail::input_atom_labels* iall,
      hierarchy::atom* atoms,
      hierarchy::chain& chain,
      bool link_to_previous,
      std::map<str4, std::vector<unsigned> >& altloc_resname_indices,
      bool residue_group_post_processing)
    {
      hierarchy::residue_group rg(
        iall->resseq_small().elems,
        iall->icode_small().elems,
        link_to_previous);
      chain.append_residue_group(rg);
      unsigned n_ag = static_cast<unsigned>(altloc_resname_indices.size());
      rg.pre_allocate_atom_groups(n_ag);
      typedef std::map<str4, std::vector<unsigned> >::const_iterator ari_it;
      boost::scoped_array<ari_it> ari_iters(new ari_it[n_ag]);
      boost::scoped_array<unsigned> first_indices(new unsigned[n_ag]);
      ari_it ari_end = altloc_resname_indices.end();
      unsigned i = 0;
      for(ari_it ari=altloc_resname_indices.begin(); ari!=ari_end; ari++,i++) {
        ari_iters[i] = ari;
        first_indices[i] = (ari->second.size() ? ari->second[0] : 0);
      }
      af::shared<std::size_t> permutation = af::sort_permutation(
        af::const_ref<unsigned>(first_indices.get(), n_ag));
      const std::size_t* perm = permutation.begin();
      char altloc[2];
      altloc[1] = '\0';
      for(i=0;i<n_ag;i++) {
        ari_it ari = ari_iters[perm[i]];
        altloc[0] = ari->first.elems[0];
        hierarchy::atom_group ag(altloc, ari->first.elems+1);
        rg.append_atom_group(ag);
        ag.pre_allocate_atoms(ari->second.size());
        typedef std::vector<unsigned>::const_iterator i_it;
        i_it i_end = ari->second.end();
        for(i_it i=ari->second.begin();i!=i_end;i++) {
          ag.append_atom(atoms[*i]);
        }
      }
      altloc_resname_indices.clear();
      if (residue_group_post_processing) {
        rg.edit_blank_altloc();
      }
    }

  } // namespace <anonymous>

  hierarchy::root
  input::construct_hierarchy(
    bool residue_group_post_processing)
  {
    af::const_ref<str8>
      model_ids = model_ids_.const_ref();
    af::const_ref<std::vector<unsigned> >
      chain_indices = chain_indices_.const_ref();
    const std::size_t* break_index = break_indices_.begin();
    const std::size_t* break_indices_end = break_indices_.end();
    unsigned next_break_index = static_cast<unsigned>(
      break_index == break_indices_end ?
        atoms_.size() : *break_index++);
    SCITBX_ASSERT(chain_indices.size() == model_ids.size());
    hierarchy::root result;
    result.pre_allocate_models(model_ids.size());
    const detail::input_atom_labels* iall = input_atom_labels_list_.begin();
    hierarchy::atom* atoms = atoms_.begin();
    unsigned next_chain_range_begin = 0;
    for(unsigned i_model=0;i_model<model_ids.size();i_model++) {
      hierarchy::model model(model_ids[i_model].elems);
      result.append_model(model);
      model.pre_allocate_chains(chain_indices[i_model].size());
      range_loop<unsigned> ch_r(
        chain_indices[i_model], next_chain_range_begin);
      for(unsigned i_chain=0;ch_r.next();i_chain++) {
        hierarchy::chain chain(iall[ch_r.begin].chain_small().elems);
        model.append_chain(chain);
        std::map<str4, std::vector<unsigned> > altloc_resname_indices;
        unsigned rg_start = ch_r.begin;
        bool link_to_previous = false;
        const char* prev_resid = 0;
        const char* prev_resname = 0;
        bool open_resname_run_has_blank_altloc = false;
        for (unsigned i_atom=ch_r.begin; i_atom!=ch_r.end; i_atom++) {
          bool is_first_after_break = (i_atom == next_break_index);
          if (is_first_after_break) {
            next_break_index = static_cast<unsigned>(
              break_index == break_indices_end ?
                atoms_.size() : *break_index++);
          }
          detail::input_atom_labels const& ial = iall[i_atom];
          const char* resid = ial.resid_begin();
          const char* resname = ial.resname_begin();
          bool curr_blank_altloc = (ial.altloc_begin()[0]==blank_altloc_char);
          if (prev_resid != 0) {
            bool is_boundary = (std::memcmp(prev_resid, resid, 5U) != 0);
            if (!is_boundary && std::memcmp(prev_resname, resname, 3U) != 0) {
              if (open_resname_run_has_blank_altloc || curr_blank_altloc) {
                is_boundary = true;
              }
              else {
                for (unsigned j_atom=i_atom+1; j_atom!=ch_r.end; j_atom++) {
                  detail::input_atom_labels const& fwd_ial = iall[j_atom];
                  const char* fwd_resid = fwd_ial.resid_begin();
                  const char* fwd_resname = fwd_ial.resname_begin();
                  if (std::memcmp(resname, fwd_resname, 3U) != 0) break;
                  if (std::memcmp(resid, fwd_resid, 5U) != 0) break;
                  if (fwd_ial.altloc_begin()[0] == blank_altloc_char) {
                    is_boundary = true;
                    break;
                  }
                }
              }
            }
            if (is_boundary) {
              append_residue_group(
                iall+rg_start,
                atoms+rg_start,
                chain,
                link_to_previous,
                altloc_resname_indices,
                residue_group_post_processing);
              rg_start = i_atom;
              link_to_previous = !is_first_after_break;
              open_resname_run_has_blank_altloc = false;
            }
            else if (is_first_after_break) {
              char buf[64];
              std::sprintf(buf,
                "Misplaced BREAK record (%s line %u).",
                source_info_.size()
                  ? (source_info_ + ",").c_str()
                  : "input",
                break_record_line_numbers[
                  break_index - break_indices_.begin() - 1]);
              throw std::runtime_error(buf);
            }
          }
          prev_resid = resid;
          prev_resname = resname;
          if (curr_blank_altloc) open_resname_run_has_blank_altloc = true;
          altloc_resname_indices[ial.confid_small()].push_back(
            i_atom-rg_start);
        }
        if (prev_resid != 0) {
          append_residue_group(
            iall+rg_start,
            atoms+rg_start,
            chain,
            link_to_previous,
            altloc_resname_indices,
            residue_group_post_processing);
        }
        if (residue_group_post_processing) {
          chain
            .merge_disconnected_residue_groups_with_pure_altloc();
        }
      }
      next_chain_range_begin = ch_r.end;
    }
    SCITBX_ASSERT(break_index == break_indices_end);
    return result;
  }

  void
  input_atoms_with_labels_generator::run(input const& inp)
  {
    af::const_ref<str8>
      model_ids = inp.model_ids_small().const_ref();
    af::const_ref<std::vector<unsigned> >
      chain_indices = inp.chain_indices().const_ref();
    SCITBX_ASSERT(chain_indices.size() == model_ids.size());
    const std::size_t* break_index = inp.break_indices().begin();
    const std::size_t* break_indices_end = inp.break_indices().end();
    unsigned next_break_index = static_cast<unsigned>(
      break_index == break_indices_end ?
        inp.atoms().size() : *break_index++);
    const detail::input_atom_labels*
      iall = inp.input_atom_labels_list().begin();
    const hierarchy::atom* atoms = inp.atoms().begin();
    unsigned next_chain_range_begin = 0;
    for(unsigned i_model=0;i_model<model_ids.size();i_model++) {
      str8 const& model_id = model_ids[i_model];
      if (!process_model(model_id)) return;
      range_loop<unsigned> ch_r(
        chain_indices[i_model], next_chain_range_begin);
      for(unsigned i_chain=0;ch_r.next();i_chain++) {
        bool is_first_in_chain = true;
        for (unsigned i_atom=ch_r.begin; i_atom!=ch_r.end; i_atom++) {
          detail::input_atom_labels const& ial = iall[i_atom];
          bool is_first_after_break = (i_atom == next_break_index);
          if (is_first_after_break) {
            next_break_index = static_cast<unsigned>(
              break_index == break_indices_end ?
                inp.atoms().size() : *break_index++);
            if (!process_break()) return;
          }
          if (!process_atom(hierarchy::atom_with_labels(
                 atoms[i_atom],
                 model_id.elems,
                 ial.chain_small().elems,
                 ial.resseq_small().elems,
                 ial.icode_small().elems,
                 ial.altloc_small().elems,
                 ial.resname_small().elems,
                 is_first_in_chain,
                 is_first_after_break))) return;
          is_first_in_chain = false;
        }
        if (!process_ter()) return;
      }
      if (!process_endmdl(model_id)) return;
      next_chain_range_begin = ch_r.end;
    }
    SCITBX_ASSERT(break_index == break_indices_end);
    process_end();
  }

}} // namespace iotbx::pdb
