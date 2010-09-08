import os, re
from iotbx.detectors import ImageFactory, url_support

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
    return os.path.join(self.cwd,self.base)
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
    if '--index_only' in arg_module.argv:
      procedure_preferences.index_only = 1
      arg_module.argv.remove('--index_only')
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
