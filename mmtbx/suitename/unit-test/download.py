import urllib.request
import os, sys

workDir = os.path.dirname(sys.argv[0]).strip()
if workDir:
    os.chdir(workDir)
    # else we must be there already
print("working directory=", os.getcwd())


def download(molecule):
  url = f"https://files.rcsb.org/download/{molecule}.pdb"
  try:
      print('downloading' , molecule)
      name, result = urllib.request.urlretrieve(url, f"./{molecule}.pdb")
      print(name)
      return True
  except urllib.error.HTTPError as err:
      if err.code == 404:
          print(f"{molecule} does not exist")
      else:
          print(f"{molecule} failed with error {err.code}: {err.reason}")
      return False
        
if __name__ == '__main__':
    molecule = sys.argv[1]
    download(molecule)
    
    
# example how it works:
# url = 'https://files.rcsb.org/download/1Q9A.pdb'
# urllib.request.urlretrieve(url, './1q9a.pdb')
