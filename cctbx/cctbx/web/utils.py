import exceptions
class FormatError(exceptions.Exception): pass

def show_input_symbol(sgsymbol, convention, label="Input"):
  if (sgsymbol != ""):
    print label, "space group symbol:", sgsymbol
    print "Convention:",
    if   (convention == "A1983"):
      print "International Tables for Crystallography, Volume A 1983"
    elif (convention == "I1952"):
      print "International Tables for Crystallography, Volume I 1952"
    elif (convention == "Hall"):
      print "Hall symbol"
    else:
      print "Default"
  print

def interpret_skip_columns(skip_columns):
  result = int(skip_columns)
  if (result < 0):
    raise ValueError, "Negative number for columns to skip."
  return result

def interpret_coordinate_line(line, skip_columns):
  flds = line.split()
  if (len(flds) < skip_columns + 3): raise FormatError, line
  coordinates = [0,0,0]
  for i in xrange(3):
    try: coordinates[i] = float(flds[skip_columns + i])
    except: raise FormatError, line
  return " ".join(flds[:skip_columns]), coordinates
