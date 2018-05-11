from __future__ import absolute_import, division, print_function

def id_conversion(path):
  return path.replace('/', '-')

smv_images = [
    "image_examples/ALS_501/als501_q4_1_001.img",
    "image_examples/SPring8_BL26B1_SaturnA200/A200_000001.img",
    "image_examples/SPring8_BL26B1_SaturnA200/A200_000002.img",
    "image_examples/APS_19ID/q315_unbinned_a.0001.img",
    "image_examples/ALS_821/q210_lyso_1_101.img",
    "image_examples/MLFSOM_simulation/fake_00001.img",
    "image_examples/ALS_831/q315r_lyso_001.img",
    "image_examples/DESY_ID141/q210_2_001.img",
    "image_examples/SSRL_bl91/q315_1_001.img",
    "image_examples/APS_24IDC/q315_1_001.img",
    "image_examples/ALS_1231/q315r_lyso_1_001.img",
    "image_examples/SRS_142/q4_1_001.img",
    "image_examples/ALS_422/lyso_041013a_1_001.img",
    "image_examples/APS_17ID/q210_1_001.img",
    "image_examples/saturn/lyso_00001.img",
]
smv_image_ids = map(id_conversion, smv_images)

tiff_images = [
    "image_examples/SPring8_BL32XU_MX225HS/ds_000045.img",
    "image_examples/SPring8_BL32XU_MX225HS/ds_000001.img",
    "image_examples/SPring8_BL44XU_MX300HE/bl44xu_lys_000002.img",
    "image_examples/SPring8_BL44XU_MX300HE/bl44xu_lys_000001.img",
    "image_examples/SPring8_BL32XU/rayonix225hs_0001.img",
    "image_examples/SPring8_BL32XU/rayonix225_0001.img",
    "image_examples/SLS_X06SA/mar225_2_001.img",
    "image_examples/CLS1_08ID1/mar225_2_E0_0001.img",
    "image_examples/SPring8_BL38B1_MX225HE/bl38b1_001.img",
    "image_examples/SPring8_BL38B1_MX225HE/bl38b1_090.img",
    "image_examples/SRS_101/mar225_001.img",
    "image_examples/SPring8_BL12B2_MX225HE/lys001_000091.img",
    "image_examples/SPring8_BL12B2_MX225HE/lys001_000001.img",
    "image_examples/SPring8_BL26B2_MX225/2sec_Al200um_000001.img",
    "image_examples/SPring8_BL26B2_MX225/2sec_Al200um_000090.img",
]
tiff_image_ids = map(id_conversion, tiff_images)

cbf_images = [
    "image_examples/ESRF_ID29/trypsin_1_0001.cbf",
    "image_examples/xia2/merge2cbf_averaged_0001.cbf",
    "image_examples/dials-190/whatev1_01_00002.cbf",
    "image_examples/dials-190/whatev1_02_00001.cbf",
    "image_examples/dials-190/whatev1_02_00002.cbf",
    "image_examples/dials-190/whatev1_03_00001.cbf",
    "image_examples/dials-190/whatev1_01_00001.cbf",
    "image_examples/dials-190/whatev1_03_00002.cbf",
    "image_examples/APS_24IDE_test/thaum-12_1_0005.cbf",
    "image_examples/APS_24IDE_test/thaum-12_1_0004.cbf",
    "image_examples/APS_24IDE_test/thaum-12_1_0006.cbf",
    "image_examples/APS_24IDE_test/thaum-12_1_0008.cbf",
    "image_examples/APS_24IDE_test/thaum-12_1_0009.cbf",
    "image_examples/APS_24IDE_test/thaum-12_1_0010.cbf",
    "image_examples/APS_24IDE_test/thaum-12_1_0007.cbf",
    "image_examples/APS_24IDE_test/thaum-12_1_0002.cbf",
    "image_examples/APS_24IDE_test/thaum-12_1_0001.cbf",
    "image_examples/APS_24IDE_test/thaum-12_1_0003.cbf",
    "image_examples/APS_24IDC/pilatus_1_0001.cbf",
    "image_examples/SLS_Eiger_16M_as_CBF/insu_with_bs_labelit_0901.cbf",
    "image_examples/SLS_Eiger_16M_as_CBF/insu_with_bs_labelit_0001.cbf",
    "image_examples/SPring8_BL41XU_PILATUS3_6M/data1_000001.cbf",
    "image_examples/SPring8_BL41XU_PILATUS3_6M/data1_000901.cbf",
    "image_examples/DLS_I02/X4_wide_M1S4_1_0001.cbf",
    "image_examples/DLS_I23/I23_P12M_alpha_0001.cbf",
    "image_examples/DLS_I23/germ_13KeV_0001.cbf",
    "image_examples/SLS_X06SA/pilatus6m_1_00001.cbf",
    "image_examples/SPring8_ADSC_SN916/Xtal17-2phi_3_015.cbf",
    "image_examples/DLS_I19/I19_P300k_00001.cbf",
    "image_examples/ED_From_TIFF/170112330001.cbf",
]
cbf_image_ids = map(id_conversion, cbf_images)

cbf_multitile_images = [
    "stills_test_data/hit-20111202210224984.cbf",
    "stills_test_data/hit-s00-20140306002935980.cbf",
    "stills_test_data/hit-s00-20140306002857363.cbf",
]
cbf_multitile_image_ids = map(id_conversion, cbf_multitile_images)

hdf5_images = [
    "image_examples/putative_imgCIF_HDF5_mapping/minicbf.h5",
]
hdf5_image_ids = map(id_conversion, hdf5_images)
