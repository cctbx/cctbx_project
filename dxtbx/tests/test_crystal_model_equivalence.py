from __future__ import absolute_import, division, print_function
from dxtbx.model.experiment_list import ExperimentListFactory
from scitbx.matrix import sqr
from dxtbx.model import Crystal
from cctbx import crystal_orientation
from libtbx.test_utils import approx_equal

"""
This test exercises and demos several pieces of cctbx/dxtbx/nanoBragg
functionality:

1) Convert the known orientation matrix used to generate a nanoBragg image to
a dxtbx crystal object
2) Convert a dxtbx crystal object to a cctbx crystal_orientation object
3) Verify that these conversions are equivalent
4) Rotate an indexing solution onto the known orientation, accounting for
the indexing happening using an alternate symmetry operation
5) Verify the mis-orientation extent between the known and calculated
orientations

"""

from cctbx.uctbx import unit_cell
from cctbx.crystal import symmetry

GFC = symmetry(
    unit_cell=unit_cell((5.876, 7.299, 29.124, 90.00, 95.79, 90.00)),
    space_group_symbol="C 1 2/c 1",
)
CB_OP_C_P = GFC.change_of_basis_op_to_primitive_setting()  # from C to P

permute = sqr((0, 0, 1, 0, 1, 0, -1, 0, 0))

json_is = """{
  "__id__": "ExperimentList",
  "experiment": [
    {
      "__id__": "Experiment",
      "identifier": "",
      "beam": 0,
      "detector": 0,
      "crystal": 0,
      "imageset": 0
    }
  ],
  "imageset": [
    {
      "__id__": "ImageSet",
      "images": [
        "step5_000009.img.gz"
      ],
      "mask": "",
      "gain": "",
      "pedestal": "",
      "dx": "",
      "dy": "",
      "params": {}
    }
  ],
  "beam": [
    {
      "direction": [
        -0.0,
        -0.0,
        1.0
      ],
      "transmission": 1.0,
      "polarization_normal": [
        0.0,
        1.0,
        0.0
      ],
      "divergence": 0.0,
      "polarization_fraction": 0.999,
      "flux": 0.0,
      "sigma_divergence": 0.0,
      "wavelength": 1.30432
    }
  ],
  "detector": [
    {
      "hierarchy": {
        "origin": [
          0.0,
          0.0,
          0.0
        ],
        "fast_axis": [
          1.0,
          0.0,
          0.0
        ],
        "name": "",
        "raw_image_offset": [
          0,
          0
        ],
        "slow_axis": [
          0.0,
          1.0,
          0.0
        ],
        "material": "",
        "mask": [],
        "thickness": 0.0,
        "mu": 0.0,
        "gain": 1.0,
        "trusted_range": [
          0.0,
          0.0
        ],
        "image_size": [
          0,
          0
        ],
        "px_mm_strategy": {
          "type": "SimplePxMmStrategy"
        },
        "identifier": "",
        "type": "",
        "children": [
          {
            "panel": 0
          }
        ],
        "pixel_size": [
          0.0,
          0.0
        ]
      },
      "panels": [
        {
          "origin": [
            -97.24,
            97.13000000000001,
            -50.0
          ],
          "fast_axis": [
            1.0,
            0.0,
            0.0
          ],
          "name": "Panel",
          "raw_image_offset": [
            0,
            0
          ],
          "slow_axis": [
            0.0,
            -1.0,
            0.0
          ],
          "material": "",
          "mask": [],
          "thickness": 0.0,
          "mu": 0.0,
          "gain": 1.0,
          "trusted_range": [
            -9.0,
            65525.0
          ],
          "image_size": [
            1765,
            1765
          ],
          "px_mm_strategy": {
            "type": "SimplePxMmStrategy"
          },
          "identifier": "",
          "type": "SENSOR_CCD",
          "pixel_size": [
            0.11,
            0.11
          ]
        }
      ]
    }
  ],
  "goniometer": [],
  "scan": [],
  "crystal": [
    {
      "__id__": "crystal",
      "real_space_a": [
        -1.6148276313609897,
        -5.437716408522997,
        -1.4401497447540712
      ],
      "real_space_b": [
        -0.3944198624524956,
        1.9800095805737359,
        -7.033859173305951
      ],
      "real_space_c": [
        28.45936960940352,
        -4.754762088402217,
        -2.934294620730966
      ],
      "space_group_hall_symbol": "-C 2yc",
      "ML_half_mosaicity_deg": 0.3692704083777256,
      "ML_domain_size_ang": 416.04559892901716
    }
  ]
}"""


def test_compare_example():
    experiments = ExperimentListFactory.from_json(json_is, check_format=False)

    for experiment in experiments:
        header = {
            "DIM": "2",
            "DENZO_X_BEAM": "97.185",
            "RANK": "0",
            "PREFIX": "step5_000009",
            "BEAM_CENTER_Y": "97.13",
            "BEAM_CENTER_X": "97.13",
            "WAVELENGTH": "1.30432",
            "OSC_START": "0",
            "ADC_OFFSET": "10",
            "BYTE_ORDER": "little_endian",
            "DIRECT_SPACE_ABC": "2.7790989649304656,3.721227037283121,0.616256870237976,-4.231398741367156,1.730877297917864,1.0247061633019547,2.9848502901631253,-4.645083818143041,28.595825588147285",
            "OSC_RANGE": "0",
            "DIALS_ORIGIN": "-97.185,97.185,-50",
            "MOSFLM_CENTER_X": "97.13",
            "MOSFLM_CENTER_Y": "97.13",
            "CLOSE_DISTANCE": "50",
            "BEAMLINE": "fake",
            "TWOTHETA": "0",
            "ADXV_CENTER_Y": "96.965",
            "ADXV_CENTER_X": "97.185",
            "HEADER_BYTES": "1024",
            "DETECTOR_SN": "000",
            "DISTANCE": "50",
            "PHI": "0",
            "SIZE1": "1765",
            "SIZE2": "1765",
            "XDS_ORGX": "884",
            "XDS_ORGY": "884",
            "DENZO_Y_BEAM": "97.185",
            "TIME": "1",
            "TYPE": "unsigned_short",
            "PIXEL_SIZE": "0.11",
        }
        # the header of the simulated image, containing ground truth orientation

        rsabc = (
            sqr([float(v) for v in header["DIRECT_SPACE_ABC"].split(",")])
            * permute.inverse()
        )
        rsa = rsabc[0:3]
        rsb = rsabc[3:6]
        rsc = rsabc[6:9]
        header_cryst = Crystal(rsa, rsb, rsc, "P1").change_basis(CB_OP_C_P.inverse())
        header_cryst.set_space_group(experiment.crystal.get_space_group())
        print("Header crystal")
        print(header_cryst)

        expt_crystal = experiment.crystal
        print("Integrated crystal")
        print(expt_crystal)

        header_ori = crystal_orientation.crystal_orientation(
            header_cryst.get_A(), crystal_orientation.basis_type.reciprocal
        )
        expt_ori = crystal_orientation.crystal_orientation(
            expt_crystal.get_A(), crystal_orientation.basis_type.reciprocal
        )

        print("Converted to cctbx")
        header_ori.show()
        expt_ori.show()

        # assert the equivalence between dxtbx crystal object and the cctbx crystal_orientation object
        assert approx_equal(header_cryst.get_U(), header_ori.get_U_as_sqr())
        assert approx_equal(header_cryst.get_A(), header_ori.reciprocal_matrix())
        assert approx_equal(
            header_cryst.get_B(),
            sqr(header_ori.unit_cell().fractionalization_matrix()).transpose(),
        )

        assert approx_equal(expt_crystal.get_U(), expt_ori.get_U_as_sqr())
        assert approx_equal(expt_crystal.get_A(), expt_ori.reciprocal_matrix())
        assert approx_equal(
            expt_crystal.get_B(),
            sqr(expt_ori.unit_cell().fractionalization_matrix()).transpose(),
        )

        cb_op_align = sqr(expt_ori.best_similarity_transformation(header_ori, 50, 1))

        print("XYZ angles", cb_op_align.r3_rotation_matrix_as_x_y_z_angles(deg=True))
        aligned_ori = expt_ori.change_basis(cb_op_align)

        U_integrated = aligned_ori.get_U_as_sqr()
        U_ground_truth = header_ori.get_U_as_sqr()

        missetting_rot = U_integrated * U_ground_truth.inverse()
        print("determinant", missetting_rot.determinant())
        assert approx_equal(missetting_rot.determinant(), 1.0)
        assert missetting_rot.is_r3_rotation_matrix()

        angle, axis = missetting_rot.r3_rotation_matrix_as_unit_quaternion().unit_quaternion_as_axis_and_angle(
            deg=True
        )
        print("Angular offset is %13.10f deg." % (angle))

        assert approx_equal(angle, 0.2609472065)
