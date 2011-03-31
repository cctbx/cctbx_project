
def pretty_filename(spotfinder,key):
  nwildcard = spotfinder.pd['template'].count('#')
  template_f='#'*nwildcard
  extensn_f='%%0%dd'%nwildcard
  path = spotfinder.pd['template'].replace(template_f,extensn_f%key)
  return path

def optionally_add_saturation(canonical_info,image):
  if image.has_key('saturation'):
    sat_info = ("%6s",image['saturation'].message(),
                      image['saturation'].format())
    canonical_info.append(sat_info)

    sat_info = ("%6s",image['saturation'].message2(),
                      image['saturation'].format2())
    canonical_info.append(sat_info)

def optionally_add_saturation_webice(canonical_info,image):
  if image.has_key('saturation'):
    sat_info = ("%5s %%",image['saturation'].message().split("%")[1],
                      "%.1f"%(100*image['saturation'].p_saturation))
    canonical_info.append(sat_info)

    sat_info = ("%7s","In-resolution overloaded spots",
                      image['saturation'].format2())
    canonical_info.append(sat_info)

def key_adaptor(mapping,key,idx=None):
  if mapping.has_key(key):
    if idx == None:
      return mapping[key]
    else: return mapping[key][idx]
  else:
    return None

def key_safe_items(image):
  return [
      ("%6d","Spot Total",key_adaptor(image,'N_spots_total')),
      ("%6d","Remove Ice",key_adaptor(image,'N_ice_free_resolution_spots')),
      ("%6d","In-Resolution Total",key_adaptor(image,'N_spots_resolution')),
      ("%6d","Good Bragg Candidates",key_adaptor(image,'N_spots_inlier')),
      ("%6d","Ice Rings",key_adaptor(image,'ice-ring_impact')),
      ("%6.2f","Method 1 Resolution",key_adaptor(image,'distl_resolution')),
      ("%6.2f","Method 2 Resolution",key_adaptor(image,'resolution')),
      ("%6.1f","Maximum unit cell",key_adaptor(image,'maxcel')),
      ("%7.3f ","<Spot model eccentricity>",key_adaptor(image,'eccen',0)),
  ]

def key_safe_items_webice(image):
  return [
      ("%7d","Initial spot picks (yellow/green)",key_adaptor(image,'N_spots_total')),
      ("%7d","Good Bragg candidates (green)",key_adaptor(image,'N_spots_inlier')),
      ("%7.3f ","Average spot model eccentricity",image['eccen'][0]),
      ("%7d","Ice rings (orange)",key_adaptor(image,'ice-ring_impact')),
      ("%5.1f &#197","Resolution estimate before indexing",key_adaptor(image,'resolution')),
      ("%5.0f &#197","Maximum unit cell edge",key_adaptor(image,'maxcel')),
  ]

def pretty_image_stats(Spotfinder,key):
    print
    image = Spotfinder.images[key]

    canonical_info = [
      ("%s","File",pretty_filename(Spotfinder,key)),
    ]

    canonical_info.extend(key_safe_items(image))

    optionally_add_saturation(canonical_info,image)

    for item in canonical_info:
      if item[2]==None:
        print "%25s : None"%item[1]
      else:
        print "%25s : %s"%(item[1],item[0]%item[2])

def webice_image_stats(Spotfinder,key):
    import StringIO
    image = Spotfinder.images[key]
    canonical_info = []
    canonical_info.extend(key_safe_items_webice(image))
    optionally_add_saturation_webice(canonical_info,image)
    g = StringIO.StringIO()
    for item in canonical_info:
      if item[2]==None:
        print >>g,"%35s: None"%item[1]
      else:
        print >>g,"%35s: %s"%(item[1],item[0]%item[2])
    ibinfac = int(Spotfinder.pd['binning'])
    if ibinfac > 1:
      print >>g,"\nImage was processed with %1dx%1d pixel binning\n to increase Viewer speed."%(ibinfac,ibinfac)
    return g.getvalue()

def notes(Spotfinder,key):
    print
    image = Spotfinder.images[key]
    if image.has_key('resolution_divisor'):
      print "Bin population cutoff for method 2 resolution: %.0f%%"%(
        100./image['resolution_divisor'])

if __name__=='__main__':
  Spotfinder = unpickle_spotfinder()

  for key in Spotfinder.pd['osc_start'].keys():
    pretty_image_stats(Spotfinder,key)
  notes(Spotfinder,Spotfinder.pd['osc_start'].keys()[0])
