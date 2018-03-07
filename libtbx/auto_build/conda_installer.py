import urllib
import subprocess
import os
import platform
import math

class conda_installer():
  def __init__(self, miniconda=None):
    if miniconda == None:
      self.miniconda='newconda/'
    self.conda_bin=self.miniconda+'bin/conda'

  def miniconda_installer(self, dry_run=True):
    if dry_run:
      print 'This is a dry run. Miniconda installer script will be downloaded and program will terminate'
    urllib.urlretrieve("https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh", "Miniconda2-latest-Linux-x86_64.sh")
    subprocess.call(['chmod', '+x', 'Miniconda2-latest-Linux-x86_64.sh'])
    if not dry_run:
      subprocess.call(['sh','./Miniconda2-latest-Linux-x86_64.sh', '-b', '-p', self.miniconda])
      if os.path.isdir(os.path.join(os.getcwd(), 'base')):
        os.unlink(os.path.join(os.getcwd(), 'base'))

      subprocess.call([self.conda_bin, 'env','remove','-y', '-n','myEnv'])
      subprocess.call([self.conda_bin, 'create','-y', '-n','myEnv'])
      subprocess.call(['ln', '-s', '%senvs/myEnv'%self.miniconda,'base'])

  def psanaconda_lcls_installer(self):

    # Not sure if platform.linux_distribution is supported beyond python 3.6
    rhel_version = platform.linux_distribution()
    rhel_version = int(float(rhel_version[1]))
    print 'INSTALLING PSANA-CONDA on LINUX VERSION = ',rhel_version
    subprocess.call([self.conda_bin, 'clean','--index-cache'])
    subprocess.call([self.conda_bin, 'install','-y', '--verbose', '--dry-run','--channel', 'lcls-rhel%s'%rhel_version, 'psana-conda'])

  def packages_installer(self, pkgs):
    subprocess.call([self.conda_bin, 'clean','--index-cache'])
    for pkg in self.pkgs:
      print 'PACKAGE BEING INSTALLED = ',pkg
      subprocess.call([self.conda_bin, 'install', '-y', '--verbose', '--dry-run', pkg])   

if __name__ == '__main__':
  conda_helper = conda_installer()
  conda_helper.miniconda_installer(dry_run=True)

