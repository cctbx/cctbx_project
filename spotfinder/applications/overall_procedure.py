from spotfinder.applications.practical_heuristics import sf2
from spotfinder.diffraction.imagefiles import Spotspickle_argument_module

SP_argument_module = Spotspickle_argument_module

def spotfinder_no_pickle(frames,s3_passthru,spot_convention):

    local_frames=frames.frames(1)

    A = frames.images[0]
    pd = {'directory':frames.filenames.FN[0].cwd,
          'template': frames.filenames.FN[0].template,
          'identifier':frames.filenames.FN[0].fileroot,
          'vendortype':A.vendortype,
          'binning':'%d'%A.bin,
          'distance':'%f'%A.distance,
          'wavelength':'%f'%A.wavelength,
          'deltaphi':'%f'%A.deltaphi,
          }

    #temp values for getting coordinate convention
    pd['pixel_size']='%f'%A.pixel_size
    pd['size1']='%f'%A.size1
    pd['ybeam'] = '%f'%A.beamy
    pd['xbeam'] = '%f'%A.beamx
    try:
      pd['twotheta'] = '%f'%A.twotheta
    except Exception:
      pd['twotheta'] = '0.0'
    pd['s3_passthru']=s3_passthru
    pd['spot_convention']=spot_convention

    Spotfinder = sf2(pd)

    for framenumber in local_frames:
      try:
        assert Spotfinder.images.has_key(framenumber)
      except Exception:
        Spotfinder.register_frames(framenumber,frames)

    return Spotfinder
