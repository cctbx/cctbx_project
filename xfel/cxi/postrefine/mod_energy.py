import sys
import numpy 

from cctbx.array_family import flex

class energy_handler(object):
  '''
  classdocs
  '''

  def __init__(self):
    '''
    Constructor
    '''
    self.img_filename = ''
    self.pulses = 0
    self.mean_energy = 0
    self.std_energy = 0
    self.overall_intensity = 0
    self.max_intensity = 0
    self.energy_at_max_intensity = 0
    self.min_intensity = 0
    self.energy_at_min_intensity = 0
    self.mean_wavelength = 0
    self.std_wavelength = 0
    self.wavelength_at_max_intensity = 0
    self.wavelength_at_min_intensity = 0
    
    self.hc_const = 1.23984193*(10**-6) #eVm
    self.m2ang = 1E10
        
  def get_energy_info(self, file_name_in_img, file_name_in_energy, pickle_filename):
    
    #find the image name for this pickle file
    
    file_img = open(file_name_in_img, 'r')
    data_img = file_img.read().split('\n')
    for line_img in data_img:
      data_img = line_img.split(' ')
      if pickle_filename.find(data_img[0]) > 0:
        self.img_filename = data_img[1]
        break
    
    if self.img_filename != '':
      file_energy = open(file_name_in_energy,'r')
      data_energy=file_energy.read().split('\n')
      for line_energy in data_energy:
        data_energy = line_energy.split(',')
        if len(data_energy) >= 12:
          if self.img_filename.find(data_energy[3]) > 0:
            self.pulses = float(data_energy[4])
            self.mean_energy = float(data_energy[5])
            self.std_energy = float(data_energy[6])
            self.overall_intensity = float(data_energy[7])
            self.max_intensity = float(data_energy[8])
            self.energy_at_max_intensity = float(data_energy[9])
            self.min_intensity = float(data_energy[10])
            self.energy_at_min_intensity = float(data_energy[11])
            
            
            #caculate wavelength based on different energy level
            self.mean_wavelength = (self.hc_const/self.mean_energy)*self.m2ang
            self.std_wavelength = self.mean_wavelength - ((self.hc_const/(self.mean_energy+self.std_energy))*self.m2ang)
            self.wavelength_at_max_intensity = (self.hc_const/self.energy_at_max_intensity)*self.m2ang
            self.wavelength_at_min_intensity = (self.hc_const/self.energy_at_min_intensity)*self.m2ang
            break
      
      #get 10 peaks for energy spectrum
      _img_filename = self.img_filename.split('/')
      _img_filename_only = _img_filename[len(_img_filename)-1]
      img_filename_only = _img_filename_only.split('.')
      
      file_ev_full = open('energy_spectrum/'+img_filename_only[0]+'.txt', 'r')
      data_ev_full = file_ev_full.read().split('\n')
      photon_counts = flex.double()
      ev_at_counts = flex.double()
      wavelength_at_counts = flex.double()
      for data_row in data_ev_full:
        data_col = data_row.split(',')
        if len(data_col) == 2:
          photon_counts.append(float(data_col[0]))
          ev_at_counts.append(float(data_col[1]))
          wavelength_at_counts.append((self.hc_const/float(data_col[1]))*self.m2ang)
      
      self.photon_counts = photon_counts
      self.photon_counts_normalized = photon_counts/ max(photon_counts)
      self.ev_at_counts = ev_at_counts
      self.wavelength_at_counts = wavelength_at_counts
      
      perm = flex.sort_permutation(photon_counts, reverse=True)
      #start adding peak (has to be far apart enough)
      ev_add_thres = 6
      photon_counts_sort = photon_counts.select(perm)
      ev_at_counts_sort = ev_at_counts.select(perm)
      
      photon_counts_sel = flex.double()
      ev_at_counts_sel = flex.double()
      n_max_peak = 5
      cn_peak = 0
      for i in range(len(photon_counts_sort)):
        if len(photon_counts_sel) == 0:
          photon_counts_sel.append(photon_counts_sort[i])
          ev_at_counts_sel.append(ev_at_counts_sort[i])
          cn_peak += 1
        else:
          flag_add = True
          for ev in ev_at_counts_sel:
            if abs(ev_at_counts_sort[i] - ev) < ev_add_thres:
              flag_add = False
              break
          
          if flag_add:
            photon_counts_sel.append(photon_counts_sort[i])
            ev_at_counts_sel.append(ev_at_counts_sort[i])
            cn_peak += 1
               
        if cn_peak == n_max_peak:
          break
      
      #calculate weight based on photon count
      max_photon = photon_counts_sel[0]
      peak_weight_sel = photon_counts_sel / max_photon
      
      self.photon_counts_list = photon_counts_sel
      self.ev_at_counts_list = ev_at_counts_sel
      self.wavelength_at_counts_list = (self.hc_const/ev_at_counts_sel)*self.m2ang
      self.ev_weight_list = peak_weight_sel
      
      
        
        
        
        
        
        
        
        
        
        
        
        
