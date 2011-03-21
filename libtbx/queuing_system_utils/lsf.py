import pickle
from subprocess import *
#import sys
import os

class Job:
   def __init__(self, name, execObj, modules=[], pythonExec='python'):
        '''name must be unique, execObj must be pickle_able and have a run method
        pythonExec is the python command to use e.g. phenix.python'''
        self.name=name
        self.execObj=execObj
        self.submitted=False
        self.pythonExec=pythonExec
        self.modules=modules

   def start(self):
        self.pickleInputFileName=self.name+'_input.pkl'
        self.pickleOutputFileName=self.name+'_output.pkl'
        pickleFile=open(self.pickleInputFileName,'wb')
        pickle.dump(self.execObj,pickleFile)
        pickleFile.close()
        self.scriptFileName=self.name+'.py'
        scriptFile=open(self.scriptFileName,'w')
        print >> scriptFile, 'import pickle'
        for m in self.modules:
           print >> scriptFile, 'from %s import *' % m
        print >> scriptFile, 'f=open("%s","rb")' % self.pickleInputFileName
        print >> scriptFile, 'execObj=pickle.load(f)'
        print >> scriptFile, 'f.close()'
        print >> scriptFile, 'result=execObj.run()'
        print >> scriptFile, 'f=open("%s","wb")' % self.pickleOutputFileName
        print >> scriptFile, 'pickle.dump(result,f)'
        print >> scriptFile, 'f.close()'
        scriptFile.close()
        cmd='bsub -K %s %s' % (
                self.pythonExec,
                self.scriptFileName)

        #print cmd,'libtbx'
        self.process = Popen(cmd,stdout=PIPE,stderr=PIPE,shell=True)
        self.submitted=True

   def isAlive(self):
        return self.process.poll() is None

   def get(self):
        if self.isAlive():
            raise Exception('Job not finished')
        f=open(self.pickleOutputFileName,'rb')
        result=pickle.load(f)
        f.close()
        if self.process.poll() != 0:
            raise Exception("process %s failed" % self.name)
        os.remove(self.pickleInputFileName)
        os.remove(self.pickleOutputFileName)
        os.remove(self.scriptFileName)
        return result

class testObj:
   def __init__(self):
        self.data=1
        self.result=[]

   def run(self):
        import math
        for i in range(10000):
                self.result.append(math.sin(i+self.data))
        return self

if __name__=='__main__':
   import time
   o=testObj()
   j=Job('j1',o,modules=['lsf'],pythonExec='phenix.python')
   j.start()
   print j.isAlive()
   while j.isAlive():
      time.sleep(10)
   r=j.get()
   print r.result[0]
   print len(r.result)
