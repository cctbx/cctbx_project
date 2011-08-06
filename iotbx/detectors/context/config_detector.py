import re,time,copy

#------------------ALS-----------------------
'''
Prior to March 2004:
Beamline 5.0.1 ADSC Q4U  s/n 401
Beamline 5.0.2 ADSC Q210 s/n 447
Beamline 5.0.3 ADSC Q4R  s/n 423
Beamline 8.2.1 ADSC Q210 s/n 445
Beamline 8.2.2 ADSC Q315 s/n 905
After March 2004:
Beamline 5.0.1 ADSC Q210 s/n 447
Beamline 5.0.2 ADSC Q315 s/n 913
Beamline 5.0.3 ADSC Q4R  s/n 423
Beamline 8.2.1 ADSC Q210 s/n 445
Beamline 8.2.1 ADSC Q315 s/n 925
Beamline 8.2.2 ADSC Q315 s/n 905
Beamline 8.3.1 ADSC Q210 s/n 442
Beamline 12.3.2 ADSC Q315 s/n 907
After Aug 2006:
Beamline 5.0.3 ADSC Q315r  s/n 923
Installed Nov 2009:
Beamline 5.0.1 ADSC Q315r s/n 931'''

known_als_detectors = [401,423,445,447,905,913,923,925,913]
known_als831_detectors = [907,442]

def als_beam_rules(iobj): #take an ADSC image object

  beam_center_convention = 1
  #The only use for this default of 1 is for the swap_beam simulation
  # in the pipeline server.  All other known ALS images will adhere to
  # the convention of 0, if the following DENZO overrides are used:

  for tag,search,datatype in [
          ('BEAM_CENTER_X','DENZO_BEAM_CENTER_X',float),
          ('BEAM_CENTER_Y','DENZO_BEAM_CENTER_Y',float),
          ('BEAM_CENTER_X','DENZO_XBEAM',float),
          ('BEAM_CENTER_Y','DENZO_YBEAM',float),
          ]:
          pattern = re.compile(search+'='+r'(.*);')
          matches = pattern.findall(iobj.header)
          if len(matches)>0:
            beam_center_convention = 0
            iobj.parameters[tag] = datatype(matches[-1])
  return beam_center_convention

def als_beamline831_rules(iobj):
  beam_center_convention = 5
  return beam_center_convention

#-------------------Other synchrotrons-------------------------
'''
Denzo = header has special records for Denzo beam
Beam = 0 is LABELIT default
Regression = there is a corresponding entry in LABELIT regression suite

Successful Indexing:
ADSC s/n 402  Web Denzo Beam5 Regression
ADSC s/n 403  SSRL BL11-3 Beam5
ADSC s/n 406  CHESS F1, Marian Szebenyi, reverse phi, Beam5
ADSC s/n 409  CHESS F1 upper upside-down detector, NOT reverse phi, Beam5
ADSC s/n 410  NSLS X9B
ADSC s/n 411  SSRL BL1-5 Beam5
ADSC s/n 413  Web [ESRF ID14-2 Q4] Beam0 Beam_very_close_to_center Regression
ADSC Q4  414  CHESS F3, Marian Szebenyi, reverse phi or Pringle-Shen, Beam 5
ADSC s/n 415  APS: BIOCARS 14-BM-C Beam0 (Web case 83995)
ADSC s/n 418  NSLS X4A
ADSC s/n 420  Web [ESRF ID14-3 Q4R] Beam0 Regression
ADSC Q4R 428  Web Beam0 Regression ESRF
ADSC s/n 429  Web Beam0 (no strong evidence for 0) Regression
ADSC 210 441  CHESS A1, Marian Szebenyi, reverse phi after 10/2004, Beam5
ADSC 210 443  APS IMCA-CAT 17ID, Xaiochun Yang, Beam0
ADSC 210 446  BNL X6A, Vivian Stojanoff, Beam5
ADSC 210 448  CHESS F2, Marian Szebenyi, reverse phi, Beam5
ADSC 210r 457 Australian Synchrotron Tom Caradoc-Davies, reverse phi, beam on center
ADSC 270 471  CHESS F1, Marian Szebenyi, reverse phi, Beam5
ADSC 210 901  SSRL BL9-2 Beam5
ADSC s/n 902  SSRL BL11-1 Beam5
ADSC s/n 903  Web unknown location Beam5
ADSC s/n 904  APS NE-CAT [24-ID-C or 24-BM-B] Beam0 (JCSG 2003Nov beam not really tested; close to center)
ADSC s/n 908  SSRL BL9-1 Beam5
ADSC s/n 910  APS BioCARS 14-BM-C installed before August 2007.
ADSC s/n 911  APS [24-ID-C or 24-BM-B] (NE-CAT), Beam0
ADSC s/n 914  APS ID19
ADSC s/n 916  APS 24-ID-E (NE-CAT), Beam0
ADSC s/n 917  [ESRF ID23-1 Q315]
ADSC s/n 918  [ESRF ID29   Q315]
ADSC s/n 919  [ESRF ID14-4 Q315]
ADSC s/n 924  ESRF BM30:A French Beamline. Reverse Phi
ADSC 315r 928 Australian Synchrotron 3ID microfocus, Tom Caradoc-Davis, reverse phi, beam on center
              Looking from the source towards the detector our goniometer is horizontal
              and on the left-hand side of the beam. Facing the goniometer, a positive
              rotation turns the air-bearing clockwise.
Unsuccessful:
ADSC s/n 415  Web Regression Submitted virus images have wrong center - off by mm
                  Impossible to determine convention.  Subsequent images indexed
                  correctly; see above.
ADSC s/n 416  Web Beam center off by many cm; no indexing possible
ADSC s/n 444  Web Error upon image read
ADSC s/n 444  [ESRF ID14-1 Q210]
ADSC s/n 910  APS ID19 Beam0 beam perfectly centered, but it looks like
                  the phi axis is vertical. Not tested. Beginning in 2005,
                  ID19 had detector s/n 914.

Olof Svenson: all ESRF beamlines have the same coordinate system.
detector serial numbes listed at
http://www.esrf.fr/UsersAndScience/Experiments/MX/Software/PXSOFT/Denzo/
ID23-2 Mar 225, serial #5
ID13   Mar 133, serial #8
'''

def ADSC910_at_BioCARS(iobj):
  if iobj.serial_number != 910: return False
  record_date = iobj.parameters["DATE"]
  record_tse = time.mktime(time.strptime(record_date))
  #pending further information, assume 910 at BioCARS beginning in 2007
  cutoff_this = time.mktime(time.strptime("Mon Jan 01 00:00:00 2007"))
  return record_tse > cutoff_this

def other_beamlines(iobj,passthru_convention):
  beam5 = [402,403,406,409,410,411,414,418,441,446,448,471,901,902,903,908]
  beam0 = [413,415,420,428,429,443,444,457,904,914,916,917,918,919,924,928]
  alld = beam5+beam0+known_als_detectors+known_als831_detectors
  if iobj.serial_number in beam5:
    beam_center_convention = 5
  elif ADSC910_at_BioCARS(iobj):
    beam_center_convention = 0
  elif iobj.serial_number not in alld:
    print "WARNING (possibly fatal): new beamline; coordinate system unknown. Please contact the authors"
    beam_center_convention = passthru_convention
  else:
    beam_center_convention = passthru_convention
  return beam_center_convention

reference_information_ADSC_detectors = [
{'type':'Q4U', 'pixels_unbinned':10616832, 'pixels_binned':2654208},
{'type':'Q4R', 'pixels_unbinned':10616832, 'pixels_binned':2654208},
{'type':'Q210', 'pixels_unbinned':33554432, 'pixels_binned':8388608},
{'type':'Q315', 'pixels_unbinned':75497472, 'pixels_binned':18874368},
]

def set_convention(value,phil_params):
  if phil_params.convention_override != None:
     phil_params.spot_convention = copy.copy(
      phil_params.convention_override)
  else:
    phil_params.spot_convention = value

def beam_center_convention_from_image_object(imageobject,phil_params):

    if imageobject.vendortype == "ADSC":
      if imageobject.serial_number in known_als_detectors:
        beam_center_convention = als_beam_rules(imageobject)
      else:
        beam_center_convention = 0
      set_convention(0,phil_params)
      if imageobject.serial_number in known_als831_detectors:
        beam_center_convention = als_beamline831_rules(imageobject)

      beam_center_convention = other_beamlines(imageobject,
        passthru_convention = beam_center_convention)

    elif imageobject.vendortype == "CBF":
      beam_center_convention = 0
      set_convention(0,phil_params)

    elif imageobject.vendortype == "MacScience":
      beam_center_convention = 0
      set_convention(0,phil_params)

    elif imageobject.vendortype == "Bruker Proteus CCD":
      beam_center_convention = 0
      set_convention(0,phil_params)

    elif imageobject.vendortype == "RigakuSaturn":
      beam_center_convention = 5
      set_convention(0,phil_params)

    elif imageobject.vendortype=="MARCCD":
      '''Explanation: there are only two test datasets:
         /net/adder/raid1/sauter/marccd/brunzelle  4 5 6 (serial#=='1')
         /net/adder/raid1/sauter/marccd/flav  1 3 (no serial #)
           For orientation, When looking at the detector from the source,
           the beamstop shadow comes in from the right of the image
           and the beam is about 15 pixels below the center.
           This is a mosaic image from our 225cm 3x3 detector - It is 3072x3072
           pixels.The format is exactly the same as the 2k marCCD images -
           just more pixels. The crystals is flavodoxin and the images were
           taken at the ESRF BM14 by Martin Walsh
         ...listed with the beam_center conventions that support these data.
         Have to choose one, so choose 1.  Later figure out if these
         were collected at different beamlines'''
      beam_center_convention = 1
      set_convention(0,phil_params)

      '''For SSRL's MarCCD, serial number 11, there is a different beam
      convention.  We assume here SSRL has the only MarCCD with this
      serial number.  There isn't any information from other beamlines
      to tell if the SSRL convention is general.'''
      try:
        from labelit.detectors.mar import CompleteMarHeader
        C = CompleteMarHeader(imageobject)
        #C.dumpHeader()
        if C.get_serial_number() in ['11']:
          beam_center_convention = 5
          #print "Using SSRL Mar CCD beam_center system"
      except Exception:pass

    elif imageobject.vendortype=="MARIP":
      beam_center_convention = 0
      set_convention(0,phil_params)
      #but note coordinate transformation going to mosflm
    elif imageobject.vendortype=="RAXIS":
      beam_center_convention = 2
      set_convention(2,phil_params)

    elif imageobject.vendortype=="npy_raw":
      beam_center_convention = 2
      set_convention(2,phil_params)

    else: beam_center_convention = None

    if imageobject.vendortype == "CBF" and \
       imageobject.size1==2527 and imageobject.size2==2463:
       imageobject.vendortype = "Pilatus-6M"
    if imageobject.vendortype == "Pilatus" and \
       imageobject.size1==1679 and imageobject.size2==1475:
       imageobject.vendortype = "Pilatus-2M"
    if imageobject.vendortype in ["Pilatus-6M","Pilatus-2M"]:
       beam_center_convention = 0
       set_convention(0,phil_params)
       if phil_params.distl.minimum_signal_height==None:
          phil_params.distl.minimum_signal_height=2.5
       if phil_params.distl.minimum_spot_area==None or \
          phil_params.distl.minimum_spot_area > 5:
          phil_params.distl.minimum_spot_area=5

    return beam_center_convention
