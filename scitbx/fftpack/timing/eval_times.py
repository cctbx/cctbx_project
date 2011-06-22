import fileinput, re

class empty: pass

def run():
  idCodes = {}
  for line in fileinput.input():
    if (not line.startswith("timing_") and not line.startswith("Time")):
      flds = re.sub(r"[\[\],]", "", line).split()
      idCode = flds[0]
      data = empty()
      data.unit_cell = [float(u) for u in flds[1:7]]
      data.resolution = float(flds[7])
      data.grid = [int(i) for i in flds[8:11]]
      data.memory = int(flds[11])
      data.iter = int(flds[12])
      data.time = float(flds[13])
      try: idCodes[idCode].append(data)
      except Exception: idCodes[idCode] = [data]
  for idCode, data_list in idCodes.items():
    assert len(data_list) == 4
    assert data_list[0].iter == 0
    assert data_list[1].iter == 1
    assert data_list[2].iter == 0
    assert data_list[3].iter == 1
    t_fftw = data_list[1].time - data_list[0].time
    t_fftpack = data_list[3].time - data_list[2].time
    mega_byte = data_list[0].memory / (1024. * 1024)
    if (t_fftw < 0 or t_fftpack < 0):
      print "#", idCode, mega_byte, t_fftw, t_fftpack
    else:
      print idCode, mega_byte, t_fftw, t_fftpack, t_fftpack/t_fftw

if (__name__ == "__main__"):
  run()
