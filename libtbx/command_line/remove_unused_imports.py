from __future__ import division
from __future__ import print_function

from libtbx import group_args

def run(args):
  if len(args) != 1:
    print("phenix.python remove_unused_args <file_with_instructions>")
    return
  file_name = args[0]

  instruction_list = get_instruction_list(file_name)
  for file_instructions in instruction_list:
    edit_one_file(file_instructions)

def edit_one_file(file_instructions):
  file_name = file_instructions.file_name
  lines = []
  for l in open(file_name).readlines():
    lines.append(l.rstrip())
  ok = True
  did_something = False
  for instruction in file_instructions.instructions:
    word = instruction.key
    i = instruction.line_number - 1
    spl = get_split_text(lines[i])
    if not spl:
      ok = False
      continue
    if spl[-1] == "\\":
      ok = False
      continue
    if not word in spl:
      ok = False
      continue
    info = remove_word(lines[i], word)
    if info.removed:
      lines[i] = info.new_line
      did_something = True
  new_text = "\n".join(lines)
  if did_something:
    f = open(file_name,'w')
    f.write(new_text + "\n")
    f.close()
    if ok:
      print("Edited %s" %(file_name))
    else:
      print("Edited %s with errors" %(file_name))
  elif ok:
    print("Nothing done with %s" %(file_name))
  else:
    print("Errors in working with %s" %(file_name))

def remove_word(line, word):
  ''' from x import word, other_word
      import word, other_word
      from word import word
  '''
  info = group_args(group_args_type = 'remove_word info',
    word = word,
    line = line,
    new_line = line,
    removed = None,)
  spl = get_split_text(line)
  if not word in spl:
    return info
  spl = remove_word_after_key(spl, word, key = 'as')
  spl = remove_word_after_key(spl, word, key = 'import')
  starting = get_blanks_at_start(line)
  new_line = starting + " ".join(spl)
  new_line = new_line.replace(", ,",",").replace("import ,", "import").replace("as,", "as")
  if new_line.endswith(","):
    new_line = new_line[:-1].rstrip()
  if new_line.endswith("as") or new_line.endswith("import"): # nothing left
    info.new_line = ""
    info.removed = True
    return info
  info.new_line = new_line
  info.removed = True
  return info

def remove_word_after_key(spl, word, key = None):
  new_spl = []
  found_key = (key is None)
  for w in spl:
    if w == key:
      found_key = True
    if (not found_key) or (w != word):
      new_spl.append(w)
    else:
      pass # skip word after key is found
  return new_spl

def get_blanks_at_start(line):
  starting = ""
  for a in line:
    if a == " ":
      starting += a
    else:
      return starting
  return starting
def get_split_text(line):
  return line.strip().replace(","," , ").split()

def get_instruction_list(file_name):
  all_instructions = []
  for line in open(file_name).readlines():
    if line.startswith("In file"):
      line = line.replace(":","")
      fn = line.split()[-1]
      file_instructions=group_args(group_args_type = 'instructions',
        file_name = fn,
        instructions= [],
       )
      all_instructions.append(file_instructions)
    elif line.find("imported at line") > -1:
       spl = line.split()
       key = spl[0]
       line_number = int(spl[-1])
       single_import = group_args(group_args_type = 'single import',
         key = key,
         line_number = line_number,
         )
       file_instructions.instructions.append(single_import)
  return all_instructions


if __name__ == "__main__":
  import sys
  run(sys.argv[1:])
