from __future__ import division, print_function, absolute_import
from scitbx.matrix import sqr, row
from rstbx.sublattice_support.change_basis import sublattice_change_of_basis
# for debug: from cctbx.crystal_orientation import crystal_orientation
from cctbx.miller import set as miller_set
from cctbx.crystal import symmetry
import copy, logging
from dials.array_family import flex
from dials.algorithms.indexing.stills_indexer import calc_2D_rmsd_and_displacements

SL = sublattice_change_of_basis(max_modulus=2)
SLT = list(SL.yield_transformations_ascending_modulus())
# reorder ersatz yield from the toolbox to conform with organized published list
# Set up the phil file so common sense lookup gives the desired transformation
SLT = [SLT[3], SLT[1], SLT[0], SLT[5], SLT[4], SLT[2], SLT[6]]
assert SLT[0].matS().elems == (2,0,0,0,1,0,0,0,1) # Double the a-axis
assert SLT[1].matS().elems == (1,0,0,0,2,0,0,0,1) # Double the b-axis
assert SLT[2].matS().elems == (1,0,0,0,1,0,0,0,2) # Double the c-axis
assert SLT[3].matS().elems == (2,1,0,0,1,0,0,0,1) # C-face centering
assert SLT[4].matS().elems == (2,0,1,0,1,0,0,0,1) # B-face centering
assert SLT[5].matS().elems == (1,0,0,0,2,1,0,0,1) # A-face centering
assert SLT[6].matS().elems == (2,1,1,0,1,0,0,0,1) # Body centering

logger = logging.getLogger(__name__)

def integrate_coset(self, experiments, indexed):
        TRANS = self.params.integration.coset.transformation

        # here get a deepcopy that we are not afraid to modify:
        experiments_local = copy.deepcopy(experiments)

        logger.info("*" * 80)
        logger.info("Coset Reflections for modeling or validating the background")
        logger.info("*" * 80)
        from dials.algorithms.profile_model.factory import ProfileModelFactory
        from dials.algorithms.integration.integrator import create_integrator

        # XXX Fixme later implement support for non-primitive lattices NKS
        base_set = miller_set( crystal_symmetry = symmetry(
            unit_cell = experiments_local[0].crystal.get_unit_cell(),
            space_group = experiments_local[0].crystal.get_space_group()),
            indices = indexed["miller_index"]
          )
        triclinic = base_set.customized_copy(
          crystal_symmetry=symmetry(unit_cell = experiments_local[0].crystal.get_unit_cell(),space_group="P1"))

        # ================
        # Compute the profile model
        # Predict the reflections
        # Create the integrator
        # This creates a reference to the experiment, not a copy:
        experiments_local = ProfileModelFactory.create(self.params, experiments_local, indexed)
        # for debug SLT[TRANS].show_summary()

        for e in experiments_local:
          e.crystal.set_space_group(triclinic.space_group())
          Astar = e.crystal.get_A()
          # debug OriAstar = crystal_orientation(Astar,True)
          # debug OriAstar.show(legend="old ")
          Astarprime = sqr(Astar)* ( sqr(SLT[TRANS]._reindex_N).transpose().inverse() )
          e.crystal.set_A(Astarprime)
          # debug OriAstarprime = crystal_orientation(Astarprime,True)
          # debug OriAstarprime.show(legend="new ")

        logger.info("Predicting coset reflections")
        logger.info("")
        predicted = flex.reflection_table.from_predictions_multi(
            experiments_local,
            dmin=self.params.prediction.d_min,
            dmax=self.params.prediction.d_max,
            margin=self.params.prediction.margin,
            force_static=self.params.prediction.force_static,
        )
        logger.info("sublattice total predictions %d"%len(predicted))

        # filter the sublattice, keep only the coset indices
        miller = predicted["miller_index"]
        # coset of modulus 2, wherein there is a binary choice
        # see Sauter & Zwart, Acta D (2009) 65:553, Table 1; select the valid coset using eqn(5).
        coset_select_algorithm_2 = flex.bool()
        M_mat = SLT[TRANS].matS() # the transformation
        M_p = M_mat.inverse()
        for idx in miller:
          H_row = row(idx)
          h_orig_setting = H_row * M_p
          on_coset=False
          for icom in h_orig_setting.elems:
            if icom.denominator() > 1: on_coset=True; break
          coset_select_algorithm_2.append(on_coset)
        predicted = predicted.select(coset_select_algorithm_2)
        logger.info("of which %d are in coset %d"%(len(predicted), TRANS))

        logger.info("")
        integrator = create_integrator(self.params, experiments_local, predicted)

        # Integrate the reflections
        integrated = integrator.integrate()

        # Delete the shoeboxes used for intermediate calculations, if requested
        if self.params.integration.debug.delete_shoeboxes and "shoebox" in integrated:
            del integrated["shoebox"]

        if self.params.output.composite_output:
            if (
                self.params.output.coset_experiments_filename
                or self.params.output.coset_filename
            ):
                assert (
                    self.params.output.coset_experiments_filename is not None
                    and self.params.output.coset_filename is not None
                )
                n = len(self.all_coset_experiments)
                self.all_coset_experiments.extend(experiments_local)
                for i, experiment in enumerate(experiments_local):
                    refls = integrated.select(integrated["id"] == i)
                    refls["id"] = flex.int(len(refls), n)
                    del refls.experiment_identifiers()[i]
                    refls.experiment_identifiers()[n] = experiment.identifier
                    self.all_coset_reflections.extend(refls)
                    n += 1
        else:
            # Dump experiments to disk
            if self.params.output.coset_experiments_filename:

                experiments_local.as_json(self.params.output.coset_experiments_filename)

            if self.params.output.coset_filename:
                # Save the reflections
                self.save_reflections(
                    integrated, self.params.output.coset_filename
                )

        rmsd_indexed, _ = calc_2D_rmsd_and_displacements(indexed)
        log_str = "coset RMSD indexed (px): %f\n" % (rmsd_indexed)
        log_str += "integrated %d\n"%len(integrated)
        for i in range(6):
            bright_integrated = integrated.select(
                (
                    integrated["intensity.sum.value"]
                    / flex.sqrt(integrated["intensity.sum.variance"])
                )
                >= i
            )
            if len(bright_integrated) > 0:
                rmsd_integrated, _ = calc_2D_rmsd_and_displacements(bright_integrated)
            else:
                rmsd_integrated = 0
            log_str += (
                "N reflections integrated at I/sigI >= %d: % 4d, RMSD (px): %f\n"
                % (i, len(bright_integrated), rmsd_integrated)
            )

        for crystal_model in experiments_local.crystals():
            if hasattr(crystal_model, "get_domain_size_ang"):
                log_str += (
                    ". Final ML model: domain size angstroms: %f, half mosaicity degrees: %f"
                    % (
                        crystal_model.get_domain_size_ang(),
                        crystal_model.get_half_mosaicity_deg(),
                    )
                )

        logger.info(log_str)
        logger.info("")
