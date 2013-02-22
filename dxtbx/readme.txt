Diffraction Experiment Toolbox
------------------------------

A cctbx-style toolbox to describe single-crystal diffraction experiments, where
a monochromatic beam is used to illuminate a sample which is rotated during
the exposure and diffraction recorded on a flat area detector.

This toolbox will include code for:

 - reading image headers
 - transforming contents of image header to standard (i.e. imgCIF) frame
 - python models of experiment
 - reading a sweep into memory using existing cctbx image reading tools in iotbx

Initially implementing to support xia2 development and make more general for
supporting other detectors.
