import urllib
import subprocess
import os
import platform
import math

class conda_installer():
  def __init__(self, miniconda=None, myEnv=None):
    if miniconda == None:
      self.miniconda='newconda/'
    self.conda_bin=self.miniconda+'bin/conda'
    if myEnv == None:
      self.myEnv = 'myEnv'

  def miniconda_installer(self, dry_run=True):
    if dry_run:
      print 'This is a dry run. Miniconda installer script will be downloaded and program will terminate'
    urllib.urlretrieve("https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh", "Miniconda2-latest-Linux-x86_64.sh")
    subprocess.call(['chmod', '+x', 'Miniconda2-latest-Linux-x86_64.sh'])
    if not dry_run:
      subprocess.call(['sh','./Miniconda2-latest-Linux-x86_64.sh', '-b', '-p', self.miniconda])
      if os.path.isdir(os.path.join(os.getcwd(), 'base')):
        os.unlink(os.path.join(os.getcwd(), 'base'))

      subprocess.call([self.conda_bin, 'env','remove','-y', '-n',self.myEnv])
      subprocess.call([self.conda_bin, 'create','-y', '-n',self.myEnv])
      subprocess.call(['ln', '-s', '%senvs/%s'%(self.miniconda, self.myEnv),'base'])

  def psanaconda_lcls_installer(self, dry_run=True, verbose=True):

    # Not sure if platform.linux_distribution is supported beyond python 3.6
    rhel_version = platform.linux_distribution()
    rhel_version = int(float(rhel_version[1]))
    print 'INSTALLING PSANA-CONDA on LINUX VERSION = ',rhel_version
    opts = [self.conda_bin, 'install','-y']
    if dry_run:
      opts +=['--dry-run']
    if verbose:
      opts +=['--verbose']
    opts+=['--name', self.myEnv, '--channel', 'lcls-rhel%s'%rhel_version, 'psana-conda']
    subprocess.call([self.conda_bin, 'clean','--index-cache'])
    subprocess.call(opts)

  def packages_installer(self, pkgs, dry_run=True,verbose=True):
    subprocess.call([self.conda_bin, 'clean','--index-cache'])
    opts = [self.conda_bin, 'install','-y']
    if dry_run:
      opts +=['--dry-run']
    if verbose:
      opts +=['--verbose']
    opts+=['--name', self.myEnv]
    for pkg in pkgs:
      print 'PACKAGE BEING INSTALLED = ',pkg
      subprocess.call(opts+[pkg])   

if __name__ == '__main__':
  conda_helper = conda_installer()
  conda_helper.miniconda_installer(dry_run=False)

