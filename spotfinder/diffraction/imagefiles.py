import os, re
from iotbx.detectors import ImageFactory, url_support
from iotbx.detectors.beam_center_convention import convert_beam_instrument_to_imageblock

# Contain information about a single file name
  # Pattern1:  valid image file names conform to the regular expression
  # fileroot_[0-9][0-9][0-9].ext
  # valid for adsc images *.img and renamed MarCCD images *.tif
pattern1a = re.compile(r'\A(?P<fileroot>.*)_(?P<number>[0-9]{3,6})\.(?P<ext>.*)\Z')
pattern1b = re.compile(r'\A(?P<fileroot>.*)\.(?P<number>[0-9]{3,6})\.(?P<ext>.*)\Z')
pattern1c = re.compile(r'\A(?P<fileroot>.*)(?P<number>[0-9]{3,6})\.(?P<ext>.*)\Z')
  # Pattern2:  valid image file names conform to the regular expression
  # fileroot.[0-9][0-9][0-9][0-9]  (between 3 & 5 numerals in the extension)
  # valid for MarCCD images from Blum's ftp site
pattern2 = re.compile(r'\A(?P<fileroot>.*)\.(?P<number>[0-9]{3,5})\Z')
pattern_small2grid = re.compile(r'\A(?P<fileroot>.*)_(?P<number1>[0-9]{1,3})_(?P<number2>[0-9]{1,3})\.(?P<ext>.*)\Z')

class FileName:
  exts = ["img","tif","tiff","image","mccd",
          "mar1200","mar1800","mar1600","mar2400","mar2000","mar3000","mar2300","mar3450",
          "cbf","osc","ipf","sfrm","edf"
         ] #Permissible filename extensions for pattern 1
  """attributes of this class are:
     base = the file name without directory path
     cwd = just the directory path
     fileroot = the dataset root, like lyso_1 in lyso_1_002.img
     number = the image number, like 2 in lyso_1_002.img
     ext = the image extension, like img in lyso_1_002.img
     template = the dataset naming pattern, like lyso_1_###.img in lyso_1_002.img
     pattern = 1: name has an alphabetic file extension.  2: the number is in the file extension
     full_url = URL, if input as such
  """
  def __init__(self,dirname,fn):
    self.base = fn
    self.cwd = dirname
  def fullpath(self):
    if self.__dict__.has_key("full_url"):
      return self.full_url
    return os.path.join(os.path.abspath(self.cwd),self.base)
  def isImageFileName(self):
    match1a = pattern1a.match(self.base)
    match1b = pattern1b.match(self.base)
    match1c = pattern1c.match(self.base)
    match_small2grid = pattern_small2grid.match(self.base)
    match2 = pattern2.match(self.base)
    if match1a!=None:
      d = match1a.groupdict()
      self.fileroot = d['fileroot']
      self.number = int(d['number'])
      self.ext = d['ext']
      self.template = "%s_%s.%s"%(self.fileroot,'#'*len(d['number']),self.ext)
      self.pattern = 1
      if self.ext.lower() in FileName.exts: return 1
    elif match1b!=None:
      d = match1b.groupdict()
      self.fileroot = d['fileroot']
      self.number = int(d['number'])
      self.ext = d['ext']
      self.template = "%s.%s.%s"%(self.fileroot,'#'*len(d['number']),self.ext)
      self.pattern = 1
      if self.ext in FileName.exts: return 1
    elif match1c!=None:
      d = match1c.groupdict()
      self.fileroot = d['fileroot']
      self.number = int(d['number'])
      self.ext = d['ext']
      self.template = "%s%s.%s"%(self.fileroot,'#'*len(d['number']),self.ext)
      self.pattern = 1
      if self.ext in FileName.exts: return 1
    elif match_small2grid!=None:
      d = match_small2grid.groupdict()
      self.fileroot = d['fileroot']
      self.number = 1000*int(d['number1'])+int(d['number2'])
      self.ext = d['ext']
      self.template = "%s_%s_%s.%s"%(self.fileroot,'#'*len(d['number1']),'#'*len(d['number2']),self.ext)
      self.pattern = 1
      if self.ext.lower() in FileName.exts: return 1
    if match2!=None:
      d = match2.groupdict()
      self.fileroot = d['fileroot']
      self.number = int(d['number'])
      self.numberlength = len(d['number'])
      self.ext = None
      self.template = "%s."%self.fileroot + "#"*self.numberlength
      self.pattern = 2
      return 1
    return 0
  def __repr__(self):
    return "FileName object(%s)"%self.fullpath()

class file_names:
  #Note: the arg_module used to be simply "sys"; this use was Deprecated so
  #      that the indexing functionality would be identical when used from
  #      the api.
  def __init__(self,arg_module):
    self.arg_module = arg_module
    self.FN = []
    if arg_module==None: return #added for directory analysis
    #file names come from either command line or current directory
    if len(self.arg_module.argv)==1:
      # Interface 1. Current directory
      self.interface1_directory(os.getcwd())
    elif os.path.isdir(self.arg_module.argv[1]):
      if len(self.arg_module.argv)==2:
      # Interface 1. Look for all files in the given directory
        self.interface1_directory(self.arg_module.argv[1])
      else:
      # Interface 2. argv gives directory plus image numbers
        self.interface2_directory_and_frames()
    else:
     if '#' in self.arg_module.argv[1]:
       self.interface4_template()
     else:
      # Interface 3. File pathnames given on command line
      # if images are taken from command line, must recalculate
      #  DISTL_pickle because images might be different each time
      self.interface3_parse_command()

  def interface3_FN_factory(self,absfile,error_message):
    cwd = os.path.dirname(absfile)
    item = os.path.basename(absfile)
    VF = FileName(cwd,item)
    if VF.isImageFileName():
      self.FN.append(VF)
    else:
      raise Exception("Input error: "+error_message)
    return VF

  def interface3_parse_command(self):
    #The assumption is that there will be one or two regular files specified
    # on command line that are valid image file names
    # If there are two, root and ext must be the same in each case.
    # In Unix, wildcards are permitted because they are expanded by the shell.

    for file in self.arg_module.argv[1:]:
      if os.path.isfile(file):
        self.interface3_FN_factory(os.path.abspath(file),error_message="File name not accepted")
      else:
        A = url_support.potential_url_request(file)
        if A.is_url_request():
          VF = self.interface3_FN_factory(A.file,error_message="URL %s not accepted"%file)
          VF.full_url = A.text
        elif file=="data_in_object":
          VF = FileName("cxi_data","present_image_0000001")
          VF.number = 1
          VF.template = "present_image_#######"
          VF.fileroot = "present_image"
          self.FN.append(VF)
        else:
          raise Exception("File not found: %s"%file)

  def __call__(self):
    return [item.fullpath() for item in self.FN]

  def frames(self):
    return [item.number for item in self.FN]

class image_files:
  def __init__(self,arg_module,verbose=True):
    self.verbose = verbose
    self.filenames = file_names(arg_module)
    self.images = []
    for indx,name in enumerate(self.filenames()):
        A = ImageFactory(name)
        self.images.append(A)

  def frames(self,wedgelimit=None): # gives the frame numbers
    if wedgelimit == None:return [item.number for item in self.filenames.FN]
    import inspect
    print "image_files.frames deprecated usage called by %s line %d; contact nksauter@lbl.gov"%(
      inspect.currentframe().f_back.f_code.co_name,
      inspect.currentframe().f_back.f_lineno)
    return [item.number for item in self.filenames.FN[0:wedgelimit]]

  def imageindex(self,indexnumber): # gives the actual image
    for s in xrange(len(self.filenames.frames())):
      if self.filenames.frames()[s]==indexnumber:
        return self.images[s]

  def imagepath(self,indexnumber): #convenience function for finding filename
    for s in xrange(len(self.filenames.frames())):
      if self.filenames.frames()[s]==indexnumber:
        return self.filenames()[s]

class spotfinder_image_files(image_files):
  def __init__(self,arg_module,phil_params,verbose=True):
    self.verbose = verbose
    self.filenames = file_names(arg_module)
    self.phil_params = phil_params
    self.images = []
    for indx,name in enumerate(self.filenames()):
        A = ImageFactory(name)
        self.site_modifications(A,self.filenames.FN[indx])
        self.images.append(A)
    self.acceptable_use_tests_basic()

  def acceptable_use_tests_basic(self):
    if self.images[0].parameters.has_key('TWOTHETA'):
      if abs(self.images[0].twotheta) < 0.02:  #round off to zero and
                                               #retain legacy behavior
        for ik in xrange(len(self.images)):
          self.images[ik].parameters['TWOTHETA']=0.0

  def site_modifications(self,imageobject,filenameobject):

    from iotbx.detectors.context.config_detector\
      import beam_center_convention_from_image_object

    beam_center_convention = beam_center_convention_from_image_object(imageobject,self.phil_params)

    #we may elect to override the beam position globally for LABELIT.
    #Case I.  The user has provided a tuple of floats, superceding all else
    if self.phil_params.autoindex_override_beam != None:
      imageobject.parameters['BEAM_CENTER_X'],\
      imageobject.parameters['BEAM_CENTER_Y']=\
      self.phil_params.autoindex_override_beam
      imageobject.beam_center_reference_frame = "imageblock"

    #Case II.  An XY convention has been defined.
    elif beam_center_convention != 0:
      convert_beam_instrument_to_imageblock(imageobject,beam_center_convention)

    if self.phil_params.autoindex_override_distance != None:
      imageobject.parameters['DISTANCE']=self.phil_params.autoindex_override_distance

    if self.phil_params.autoindex_override_wavelength != None:
      imageobject.parameters['WAVELENGTH']=self.phil_params.autoindex_override_wavelength

    if self.phil_params.autoindex_override_deltaphi != None:
        if self.verbose:
          print "Overriding deltaphi not fully supported: contact authors"
        imageobject.parameters['OSC_RANGE']=self.phil_params.autoindex_override_deltaphi

    # override twotheta angle
    if self.phil_params.autoindex_override_twotheta != None:
      imageobject.parameters['TWOTHETA']=\
        self.phil_params.autoindex_override_twotheta

    if self.phil_params.image_specific_osc_start != None:
        imageobject.parameters['OSC_START']= \
          eval("(%s)(%d)"%(
          self.phil_params.image_specific_osc_start,filenameobject.number))

    #take care of unbinned Quantum 315
    if (self.phil_params.distl_permit_binning and \
      imageobject.size1 > 4000) or \
      self.phil_params.distl_force_binning:
      imageobject.setBin(2)
      self.phil_params.distl.minimum_spot_area = min(
        self.phil_params.distl.minimum_spot_area,
        self.phil_params.distl_binned_image_spot_size)

    if imageobject.vendortype=="MARCCD":
      #This section corrects for the fact that ESRF writes the mar ccd header
      #  with beam center in mm instead of pixels.
      detector_center_in_mm = 0.5*imageobject.size1*imageobject.pixel_size
      one_tenth_error = 0.1*detector_center_in_mm

      #offset between given beam and detector center
      import math
      def distance(a,b):
        return math.sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1]))
      offset1=distance( (detector_center_in_mm,detector_center_in_mm),
                        (imageobject.beamx,imageobject.beamy) )
      if offset1>one_tenth_error:
        newx = imageobject.beamx/imageobject.pixel_size
        newy = imageobject.beamy/imageobject.pixel_size
        #offset between corrected beam and detector center
        offset2=distance( (detector_center_in_mm,detector_center_in_mm),
                        (newx,newy) )
        if offset2<one_tenth_error:
          imageobject.parameters['BEAM_CENTER_X'] = newx
          imageobject.parameters['BEAM_CENTER_Y'] = newy
          #Furthermore the x and y are transposed in the one example we've been given
          convert_beam_instrument_to_imageblock(imageobject,
            beam_center_convention,force=True)
          if self.verbose:
            print "Mar CCD image appears to have beam center %.2f %.2f in mm instead of pixels"%(
            imageobject.beamx,imageobject.beamy)

class Spotspickle_argument_module:  #almost verbatim copy from procedure.py
  def __init__(self,directory,framelist=[]):
    self.argv = ['SP_argument']
    self.argv.append(directory)
    for item in framelist:
      self.argv.append('%d'%item)

def quick_image(filepath):
  '''A convenience factory function to return a single image object
     given a single file path, with implicit execution of all the file
     checking machinery of the ImageFiles class.  As used here it is
     an almost-thin wrapper around iotbx ImageFactory.'''
  argument_module = Spotspickle_argument_module(filepath)
  frames = image_files(argument_module)
  return frames.images[0]
