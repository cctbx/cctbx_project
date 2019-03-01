from __future__ import absolute_import, division, print_function
from scitbx import matrix
from cctbx.uctbx import unit_cell
from cctbx.sgtbx import space_group as SG


class CrystalFactory(object):
    @staticmethod
    def from_dict(d, t=None):
        """Convert the dictionary to a crystal model

        Params:
            d The dictionary of parameters
            t The template dictionary to use

        Returns:
            The crystal model

        """
        # If None, return None
        if d == None:
            if t == None:
                return None
            else:
                return from_dict(t, None)
        elif t != None:
            d = dict(t.items() + d.items())

        # Create the model from the dictionary
        if "ML_half_mosaicity_deg" in d:
            assert "ML_domain_size_ang" in d
            if d["ML_half_mosaicity_deg"] is None or d["ML_domain_size_ang"] is None:
                assert (
                    d["ML_half_mosaicity_deg"] is None
                    and d["ML_domain_size_ang"] is None
                )
            else:
                if "mosaicity" in d and d["mosaicity"] > 0:
                    print(
                        "Warning, two kinds of mosaicity found. Using Sauter2014 model"
                    )
                from dxtbx.model import MosaicCrystalSauter2014

                return MosaicCrystalSauter2014.from_dict(d)
        if "mosaicity" in d:
            from dxtbx.model import MosaicCrystalKabsch2010

            return MosaicCrystalKabsch2010.from_dict(d)
        else:
            from dxtbx.model import Crystal

            return Crystal.from_dict(d)

    @staticmethod
    def from_mosflm_matrix(
        mosflm_A_matrix, unit_cell=None, wavelength=None, space_group=None
    ):
        """
        Create a crystal_model from a Mosflm A matrix (a*, b*, c*); N.B. assumes
        the mosflm coordinate frame::

                                                         /!
                            Y-axis                      / !
                              ^                        /  !
                              !                       /   !
                              !                      /    !
                              !   /                 /  Xd !
                              !  /                 / * ^  !
                              ! /                  ! 3 !  !
                              !/      X-ray beam   !   !  !
                              /------------------------/--!---->X-axis
                             /                     !  / *1!
                          <-/-                     ! /    !
                           /  \+ve phi             ! Yd  /
                          /   /                    ! 2  /
                         /                         ! * /
                        Z-axis                  Ys ^ _/
                      Rotation                     ! /| Xs
                       axis                        !/
                                                   O

        Also assume that the mosaic spread is 0.

        :param mosflm_A_matrix: The A matrix in Mosflm convention.
        :type mosflm_A_matrix: tuple of floats
        :param unit_cell: The unit cell parameters which are used to determine the
                          wavelength from the Mosflm A matrix.
        :type unit_cell: cctbx.uctbx.unit_cell
        :param wavelength: The wavelength to scale the A matrix
        :type wavelength: float
        :param space_group: If the space group is None then the space_group will
                            be assigned as P1
        :type space_group: cctbx.sgtbx.space_group
        :returns: A crystal model derived from the given Mosflm A matrix
        :rtype: :py:class:`crystal_model`
        """
        from dxtbx.model import Crystal

        if not space_group:
            space_group = SG("P1")

        A_star = matrix.sqr(mosflm_A_matrix)
        A = A_star.inverse()

        if unit_cell:
            unit_cell_constants = unit_cell.parameters()
            a = matrix.col(A.elems[0:3])
            wavelength = unit_cell_constants[0] / a.length()
            A *= wavelength
        elif wavelength:
            A *= wavelength
        else:
            # assume user has pre-scaled
            pass

        a = A.elems[0:3]
        b = A.elems[3:6]
        c = A.elems[6:9]
        rotate_mosflm_to_imgCIF = matrix.sqr((0, 0, 1, 0, 1, 0, -1, 0, 0))
        _a = rotate_mosflm_to_imgCIF * a
        _b = rotate_mosflm_to_imgCIF * b
        _c = rotate_mosflm_to_imgCIF * c

        return Crystal(_a, _b, _c, space_group=space_group)
