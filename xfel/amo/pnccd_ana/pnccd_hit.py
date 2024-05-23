from __future__ import absolute_import, division, print_function
import numpy as np
from libtbx.mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

class hit(object):


    def __init__(self):
        pass

    def endrun(self):
        obj={'endrun':True}
        comm.send(obj,dest=0,tag=rank)


    def send(self, nr, image = None, saxs = None, c2 = None, ind = None):


        if image is not None :
           obj_image={'nr':nr,'shape':image.shape,'endrun':False}
           comm.send(obj_image,dest=0,tag=rank)
           comm.Send(image,dest=0,tag=rank)


        if saxs is not None :
           obj_saxs={'nr':nr,'shape':saxs.shape,'endrun':False}
           comm.send(obj_saxs,dest=0,tag=rank)
           comm.Send(saxs,dest=0,tag=rank)


        if c2 is not None :
           obj_c2={'nr':nr,'shape':c2.shape,'endrun':False}
           comm.send(obj_c2,dest=0,tag=rank)
           comm.Send(c2,dest=0,tag=rank)


        if ind is not None :
           obj_ind={'nr':nr,'shape':ind.shape,'endrun':False}
           comm.send(obj_ind,dest=0,tag=rank)
           comm.Send(ind,dest=0,tag=rank)



    def recv(self):


        status     = MPI.Status()
        self.myobj = comm.recv(source=MPI.ANY_SOURCE,tag=MPI.ANY_TAG,status=status)
        recvRank   = status.Get_source()


        if not hasattr(self, 'total_ave'):
           adim  = None
        else:
           adim  = self.total_ave[0].shape

        if not hasattr(self, 'total_c2'):
           cdim  = None
        else:
           cdim  = self.total_c2[0].shape

        if not hasattr(self, 'total_saxs'):
           sdim  = None
        else:
           sdim  = self.total_saxs[0].shape

        if not hasattr(self, 'total_ind'):
           idim  = None
        else:
           idim  = self.total_ind[0].shape


        if self.myobj['endrun'] == False :
           # Average image
           if self.myobj['shape'] == adim :
               self.total_ev_a[recvRank-1]     = int(self.myobj['nr'])
               self.total_ave[recvRank-1]      = np.empty(self.myobj['shape'])
               comm.Recv(self.total_ave[recvRank-1],source=recvRank,tag=MPI.ANY_TAG)
           # C2 data
           elif self.myobj['shape'] == cdim :
               self.total_ev_c[recvRank-1]     = int(self.myobj['nr'])
               self.total_c2[recvRank-1]       = np.empty(self.myobj['shape'])
               comm.Recv(self.total_c2[recvRank-1],source=recvRank,tag=MPI.ANY_TAG)
           # SAXS data
           elif self.myobj['shape'] == sdim :
               self.total_ev_s[recvRank-1]     = int(self.myobj['nr'])
               self.total_saxs[recvRank-1]     = np.empty(self.myobj['shape'])
               comm.Recv(self.total_saxs[recvRank-1],source=recvRank,tag=MPI.ANY_TAG)
           # Ind
           elif self.myobj['shape'] == idim :
               self.total_ev_i[recvRank-1]     = int(self.myobj['nr'])
               self.total_ind[recvRank-1]      = np.empty(self.myobj['shape'])
               comm.Recv(self.total_ind[recvRank-1],source=recvRank,tag=MPI.ANY_TAG)



        return (self.myobj['endrun'])
