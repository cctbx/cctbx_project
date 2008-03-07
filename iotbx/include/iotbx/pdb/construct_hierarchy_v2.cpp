#include <iotbx/pdb/input.h>
#include <scitbx/misc/fill_ranges.h>
#include <scitbx/array_family/sort.h>
#include <boost/format.hpp>
#include <boost/scoped_array.hpp>

namespace iotbx { namespace pdb {

  af::shared<hierarchy_v2::atom>
  input::atoms_v2()
  {
    unsigned n_atoms = static_cast<unsigned>(atoms_.size());
    if (atoms_v2_.size() == 0 && atoms_.size() != 0) {
      SCITBX_ASSERT(atom_serial_number_strings_.size() == atoms_.size());
      atoms_v2_.reserve(n_atoms);
      hierarchy_v1::atom* a = atoms_.begin();
      std::string* s = atom_serial_number_strings_.begin();
      for(unsigned i=0;i<n_atoms;i++,a++,s++) {
        const hierarchy_v1::atom_data* d = a->data.get();
        atoms_v2_.push_back(hierarchy_v2::atom(
          d->name.elems, d->segid.elems,
          d->element.elems, d->charge.elems, s->c_str(),
          d->xyz, d->sigxyz,
          d->occ, d->sigocc,
          d->b, d->sigb,
          d->uij, d->siguij,
          d->hetero));
      }
    }
    return atoms_v2_;
  }

  namespace {

    void
    append_residue_group(
      const input_atom_labels* iall,
      hierarchy_v2::atom* atoms,
      hierarchy_v2::chain& chain,
      std::map<str4, std::vector<unsigned> >& altloc_resname_indices)
    {
      hierarchy_v2::residue_group rg = chain.new_residue_group(
        iall->resseq_small().elems,
        iall->icode_small().elems);
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
        hierarchy_v2::atom_group ag = rg.new_atom_group(
          altloc, ari->first.elems+1);
        ag.pre_allocate_atoms(ari->second.size());
        typedef std::vector<unsigned>::const_iterator i_it;
        i_it i_end = ari->second.end();
        for(i_it i=ari->second.begin();i!=i_end;i++) {
          ag.append_atom(atoms[*i]);
        }
      }
      altloc_resname_indices.clear();
    }

  } // namespace <anonymous>

  hierarchy_v2::root
  input::construct_hierarchy_v2(
    bool residue_group_post_processing)
  {
    af::const_ref<int>
      model_numbers = model_numbers_.const_ref();
    af::const_ref<std::vector<unsigned> >
      chain_indices = chain_indices_.const_ref();
    SCITBX_ASSERT(chain_indices.size() == model_numbers.size());
    hierarchy_v2::root result;
    result.new_models(model_numbers.size());
    const input_atom_labels* iall = input_atom_labels_list_.begin();
    atoms_v2(); // to fill array
    hierarchy_v2::atom* atoms = atoms_v2_.begin();
    unsigned next_chain_range_begin = 0;
    for(unsigned i_model=0;i_model<model_numbers.size();i_model++) {
      hierarchy_v2::model model = result.models()[i_model];
      model.data->id = (boost::format("%4d") % model_numbers[i_model]).str();
      model.new_chains(chain_indices[i_model].size());
      range_loop<unsigned> ch_r(
        chain_indices[i_model], next_chain_range_begin);
      for(unsigned i_chain=0;ch_r.next();i_chain++) {
        hierarchy_v2::chain chain = model.chains_begin()[i_chain];
        chain.data->id = iall[ch_r.begin].chain();
        // convert break_indices to break_range_ids
        boost::scoped_array<unsigned>
          break_range_ids(new unsigned[ch_r.size]);
        scitbx::misc::fill_ranges(
          ch_r.begin, ch_r.end,
          break_indices_.begin(), break_indices_.end(),
          break_range_ids.get());
        unsigned rg_start = ch_r.begin;
        std::map<str4, std::vector<unsigned> > altloc_resname_indices;
        const char* prev_resid = 0;
        for (unsigned i_atom=ch_r.begin; i_atom!=ch_r.end; i_atom++) {
          const char* resid = iall[i_atom].resid_begin();
          if (   prev_resid != 0
              && std::memcmp(prev_resid, resid, 5U) != 0) {
            append_residue_group(
              iall+rg_start, atoms+rg_start, chain, altloc_resname_indices);
            rg_start = i_atom;
          }
          prev_resid = resid;
          altloc_resname_indices[iall[i_atom].altloc_resname_small()]
            .push_back(i_atom-rg_start);
        }
        if (prev_resid != 0) {
          append_residue_group(
            iall+rg_start, atoms+rg_start, chain, altloc_resname_indices);
        }
        if (residue_group_post_processing) {
          chain \
            .merge_disconnected_residue_groups_with_pure_altloc();
          chain \
            .split_residue_groups_with_mixed_resnames_but_only_blank_altloc();
        }
      }
      next_chain_range_begin = ch_r.end;
    }
    return result;
  }

}} // namespace iotbx::pdb
