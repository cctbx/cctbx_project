from __future__ import absolute_import, division, print_function
# PDB manipulation tools for ensemble models
# Tom Burnley

import sys, re, pickle
from six.moves import input

def remove_HOH(start_pdb, noHOH_pdb):

  print("Removing HOH...")

  openpdb = open(start_pdb,"r")
  noHOHpdb = open(noHOH_pdb, 'w')
  find_HOH = re.compile('HOH')

  for line in openpdb:
    if not find_HOH.search(line):
      noHOHpdb.write(line)

def find_number_model(open_pdb):
  find_MODEL = re.compile('MODEL')
  max_model_number = 0
  for line in open_pdb:
    if find_MODEL.match(line):
      model_number = int(line.split()[-1])
      if model_number > max_model_number:
        max_model_number = model_number

  print("Number of models in PDB file : ", max_model_number)
  return max_model_number

def find_header(open_pdb, num_pdb):
  print("Getting header info...................")
  find_CRYST1 = re.compile('CRYST')
  find_SCALE = re.compile('SCALE')
  for line in open_pdb:
    if find_CRYST1.search(line):
      num_pdb.write(line)
    if find_SCALE.search(line):
      num_pdb.write(line)

def model_parse(open_pdb, num_pdb, model_num, new_model_num, start_pdb):
  model_num = str(model_num)
  find_MODEL = re.compile('MODEL')
  find_ENDMDL = re.compile('ENDMDL')
  open_pdb = open(start_pdb,"r")
  for line in open_pdb:
    if find_MODEL.search(line):
      items = line.split()
      if items[1] == model_num:
        newN = str(new_model_num)
        if new_model_num < 10:
          new_line = str('MODEL        '+ newN+'\n')
        elif new_model_num < 100:
          new_line = str('MODEL       '+ newN+'\n')
        elif new_model_num < 1000:
          new_line = str('MODEL      '+ newN+'\n')
        num_pdb.write(new_line)
        for line in open_pdb:
          if find_ENDMDL.search(line):
            num_pdb.write(line)
            return
          else:
            num_pdb.write(line)

def remove_specific_HOH(open_pdb,open_water_list,wat_pdb):
  water_list = []
  for wat_num in open_water_list:
    x = wat_num.split()
    water_list.append(x[0])

  find_HOH = re.compile('HOH')

  for line in open_pdb:
    if not find_HOH.search(line):
      wat_pdb.write(line)
    else:
      linesplit = line.split()
      if (len(linesplit) > 5):
        if (linesplit[4] in water_list) or (linesplit[5] in water_list):
          print(linesplit)
          wat_pdb.write(line)

def parse_bfactor_occ_infomation(open_pdb, start_pdb):
  print("Parsing B and Occ information")
  find_MODEL = re.compile('MODEL')
  find_ENDMDL = re.compile('ENDMDL')
  #
  b_total_array = []
  b_model_array = []
  q_total_array = []
  q_model_array = []
  for line in open_pdb:
    if find_MODEL.search(line):
      if len(b_model_array) > 0:
        b_total_array.append(b_model_array)
        b_model_array=[]
        q_total_array.append(q_model_array)
        q_model_array=[]
    items = line.split()
    if items[0] == 'ATOM' or items[0] == 'HETATM':
      b_model_array.append(float(items[-3]))
      q_model_array.append(float(items[-4]))
  #Save last model b's
  b_total_array.append(b_model_array)
  q_total_array.append(q_model_array)
  b_and_q={'b_array' : b_total_array, 'q_array' : q_total_array}
  #Pickle b_and_q information
#  pickle_name = start_pdb+"_b_q_pickle.pkl"
#  pickle_jar = open(pickle_name, 'wb')
#  pickle.dump(b_and_q, pickle_jar)
#  pickle_jar.close()
  print(b_and_q['b_array'])

def parse_specific_pdb(open_pdb, start_pdb, num_pdb, last_model, new_model_num=1):
  model_num = str(last_model)
  find_MODEL = re.compile('MODEL')
  find_ENDMDL = re.compile('ENDMDL')
  open_pdb = open(start_pdb,"r")
  for line in open_pdb:
    if find_MODEL.search(line):
      items = line.split()
      if items[1] == model_num:
        newN = str(new_model_num)
        if new_model_num < 10:
          new_line = str('MODEL        '+ newN+'\n')
        elif new_model_num < 100:
          new_line = str('MODEL       '+ newN+'\n')
        elif new_model_num < 1000:
          new_line = str('MODEL      '+ newN+'\n')

        num_pdb.write(new_line)
        for line in open_pdb:
          if find_ENDMDL.search(line):
            num_pdb.write(line)
            return
          else:
            num_pdb.write(line)

def col_swap(open_pdb, open_second_pdb, get_col, replace_col, new_pdb):
  if get_col == 'q':
    get_range = [54,60]
  else:
    get_range = [60,66]

  if replace_col == 'q':
    replace_range = [54,60]
  else:
    replace_range = [60,66]

  get_col_info=[]
  for line in open_second_pdb:
    if len(line) > 66:
      get_col_info.append(line[get_range[0]:get_range[1]])

  n = 0
  for line in open_pdb:
    if len(line) > 66:
      new_pdb.write((line[0:(replace_range[0])] + get_col_info[n] + line[(replace_range[1]):]))
      n+=1
    else:
      new_pdb.write(line)

  print("Replaced : ",  replace_col, " with : ", get_col)

def remove_anisou(start_pdb, fin_pdb):
  print("Removing ANISOU...")
  openpdb = open(start_pdb,"r")
  finpdb = open(fin_pdb, 'w')
  find_ANISOU = re.compile('ANISOU')

  for line in openpdb:
    if not find_ANISOU.search(line):
      finpdb.write(line)

def phxgro_to_phx(start_pdb, fin_pdb):
  openpdb = open(start_pdb,"r")
  finpdb = open(fin_pdb, 'w')

  #Change res names
  residue_to_change = {'CY2':'CYS', 'HIB':'HIS', 'LYH':'LYS', 'TRY':'TRP'}

  for line in openpdb:
    if len(line) > 66:
      residue_name =  line[17:20]
      if residue_name in residue_to_change:
        line = line[0:17] + residue_to_change[residue_name] + line[20:]
      if residue_name == 'HOH':
        water_atom_name = line[13:16]
        new_w_a_n = water_atom_name.replace('W', "")
        new_w_a_n = new_w_a_n + " "
        line = line[0:13] + new_w_a_n + line[16:]

    finpdb.write(line)

def phx_to_phxgro(start_pdb, fin_pdb):
  openpdb = open(start_pdb,"r")
  finpdb = open(fin_pdb, 'w')

  #Change res names
  residue_to_change = {'CYS':'CY2', 'HIS':'HIB', 'LYS':'LYH', 'TRP':'TRY'}

  for line in openpdb:
    if len(line) > 66:
      residue_name =  line[17:20]
      if residue_name in residue_to_change:
        line = line[0:17] + residue_to_change[residue_name] + line[20:]
      if residue_name == 'HOH':
        water_atom_name = line[13:16]
        if water_atom_name == "O  ":
          new_w_a_n = "OW "
        if water_atom_name == "H1 ":
          new_w_a_n = "HW1"
        if water_atom_name == "H2 ":
          new_w_a_n = "HW2"
        line = line[0:13] + new_w_a_n + line[16:]

    finpdb.write(line)

def extra_conformation(start_pdb, fin_pdb):
  open_pdb = open(start_pdb,"r")
  fin_pdb = open(fin_pdb, 'w')

  find_ATOM = re.compile('ATOM')

  for line in open_pdb:
    if len(line) > 66 and find_ATOM.search(line):
      print(line)
      new_lineA = line[0:16] + 'A' + line[17:56] + '0.50' + line[60:]
      new_lineB = line[0:16] + 'B' + line[17:56] + '0.50' + line[60:]
      print(new_lineA)
      print(new_lineB)
      fin_pdb.write(new_lineA)
      fin_pdb.write(new_lineB)
    else:
      fin_pdb.write(line)
def main():

  print("\n\n\
=========================== Ensemble PDB Display Tools =========================")
  if len(sys.argv) < 2:
    print("\n\nSpecify pdb file at command line:\n")
    print("    python ens_pdb.py my.pdb\n")
    return
  else:
    print("1. Downsample ensemble")
    print("2. Remove explicit solvent atoms")
    print("3. Downsample and remove solvent")
#    print "4. Remove specific water from list"
#    print "5. Parse B-factor and Occupancy info"
#    print "6. Last ensemble model to single pdb"
#    print "7. Swap B and Q cols"
#    print "8. Remove ANISOU lines"
#    print "9. Convert phx_gro to phx"
#    print "10. Convert phx to phx_gro"
#    print "11. Add extra conformation"
    option = str(input("Input option number : "))

    if option == str("1") or option == str("3"):
      start_pdb = sys.argv[1]
      open_pdb = open(start_pdb,"r")

      number_models = find_number_model(open_pdb)

      target_number = input("Number models for output (default 25) : ")
      try:
        target_number = int(target_number)
        div = int(number_models / target_number)
      except ValueError:
        target_number = 25
        div = int(number_models / 25)

      fin_pdb = str(str(target_number) + "_dwnsmp_" + start_pdb)
      num_pdb = open(fin_pdb, 'w')

      open_pdb = open(start_pdb,"r")
      find_header(open_pdb = open_pdb, num_pdb = num_pdb)

      model_num = int(1)
      new_model_num = int(1)

      while new_model_num <= target_number:
        print("Parsing model  :", model_num)
        model_parse(open_pdb, num_pdb, model_num, new_model_num, start_pdb)
        model_num = model_num + div
        new_model_num += 1

      num_pdb.close()
      if option == str("1"):
        return

    if option == str("2"):
      start_pdb = sys.argv[1]
      fin_pdb = str("noHOH_" + start_pdb)
      remove_HOH(start_pdb, fin_pdb)
      return

    if option == str("3"):
      noHOH_pdb = str("noHOH_" + fin_pdb)
      remove_HOH(fin_pdb, noHOH_pdb)
      return

    if option == str("4"):
      print("Remove specific HOH")
      if len(sys.argv) < 3:
        print("\n\nERROR - specify water list at command line after PDB file!\n\n")
        main()
      start_pdb = sys.argv[1]
      open_pdb = open(start_pdb,"r")
      water_list = sys.argv[2]
      open_water_list = open(water_list,"r")
      wat_pdb = open('goodHOH.pdb', 'w')
      remove_specific_HOH(open_pdb,open_water_list,wat_pdb)

    if option == str("5"):
      start_pdb = sys.argv[1]
      open_pdb = open(start_pdb,"r")
      parse_bfactor_occ_infomation(open_pdb, start_pdb)

    if option == str("6"):
      start_pdb = sys.argv[1]
      open_pdb = open(start_pdb,"r")
      last_model = input("Parse x model)...     ")
      fin_pdb = str(last_model + "_" + start_pdb)
      num_pdb = open(fin_pdb, 'w')
      find_header(open_pdb, num_pdb)
      parse_specific_pdb(open_pdb, start_pdb, num_pdb, last_model, new_model_num=1)

    if option == str("7"):
      start_pdb = sys.argv[1]
      open_pdb = open(start_pdb,"r")
      if len(sys.argv) > 2:
        second_pdb = sys.argv[2]
      else:
        second_pdb = input("Second pdb import infomation from..... : ")
      open_second_pdb = open(second_pdb,"r")
      get_col = input("Import 'b' or 'q' column from second pdb  : ")
      if get_col != 'q':
        get_col = 'b'
      print("Importing col : ", get_col)
      replace_col = input("Replace 'b' or 'q' column from first pdb  : ")
      if replace_col != 'q':
        replace_col = 'b'
      new_pdb_name = start_pdb + "_swap_" + replace_col + "_with_" + get_col + ".pdb"
      new_pdb = open(new_pdb_name, 'w')
      col_swap(open_pdb, open_second_pdb, get_col, replace_col, new_pdb)

    if option == str("8"):
      start_pdb = sys.argv[1]
      fin_pdb = str("noANISOU_" + start_pdb)
      remove_anisou(start_pdb, fin_pdb)
      return

    if option == str("9"):
      start_pdb = sys.argv[1]
      fin_pdb = str("phxgro_to_phx_" + start_pdb)
      phxgro_to_phx(start_pdb, fin_pdb)
      return

    if option == str("10"):
      start_pdb = sys.argv[1]
      fin_pdb = str("phx_to_phxgro" + start_pdb)
      phx_to_phxgro(start_pdb, fin_pdb)
      return

    if option == str("11"):
      start_pdb = sys.argv[1]
      fin_pdb = str("extra_con" + start_pdb)
      extra_conformation(start_pdb, fin_pdb)
      return

    else:
      print("\n\nPlease choose from list : ")
      main()

if __name__ == '__main__':
  main()
