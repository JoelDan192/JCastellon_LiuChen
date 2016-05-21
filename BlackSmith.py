import numpy as np
from Hamiltonian import HamiltonianWorker


##A small utility
import time
import sys

def timeit(method):

    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        print '%r (%r, %r) %2.2f sec' % \
              (method.__name__, args, kw, te-ts)
        return result

    return timed

class BlackSmith(object):

    '''The One who does annealing'''
    def __init__(self,schedule=None,q=2,n=None,desc_thresh=0.98,T_Period=10,workline=None,AdjM=np.array([])):
        self.q=q
        self.n=n
        self.desc_thresh = desc_thresh
        self.T_Period = T_Period
        self.workline = workline
        self.eps = np.finfo(float).eps
        self.AdjMatrix = AdjM
        if schedule not in ['Cauchy','Baseline','Exp']:
            sys.exit('Unknown temperature schedule')
        self.schedule = schedule
        self.X=None
        self.E=None

    def work(self):
        'Print annealing schedule started...'

        for mean_deg,temperature in self.workline:
            E,X,trace = self.__anneal(c=mean_deg,temp=temperature)
            yield (mean_deg,temperature,E,X,trace)

    @timeit
    def __anneal(self,c=5,temp=None):
        '''One annealing run corresponds to one Graph G(n,c/n)'''
        steps = 0
        T = temp
        coloring_trace = []
        worker = HamiltonianWorker(q=self.q,n=self.n,c=c,AdjM=self.AdjMatrix)

        if not temp:
            T,dsteps = worker.bootstrap()
            steps +=dsteps

        if self.schedule == 'Exp':
            try:
                E,desc_ratio,stct_dsc_ratio = worker.walk(T,self.T_Period)
                steps += self.T_Period
                coloring_trace.append(E)
            except:
                print 'Could not finish annealing run. T hit 0.0'
                return coloring_trace

            print 'Desc ratio %.2f'%desc_ratio

            while  self.eps < (self.desc_thresh - desc_ratio):
                T = T * 1.5

                try:
                    E,desc_ratio,stct_dsc_ratio = worker.walk(T,self.T_Period)
                    steps += self.T_Period
                    coloring_trace.append(E)
                except:
                    print 'Could not finish annealing run. T hit 0.0'
                    return coloring_trace

            print 'Desc ratio %.2f'%desc_ratio

            while T > self.eps :
                T = T*0.95
                try:
                    E,desc_ratio,stct_dsc_ratio = worker.walk(T,self.T_Period)
                    steps += self.T_Period
                    coloring_trace.append(E)
                except:
                    print 'Could not finish annealing run. T hit 0.0'
                    return coloring_trace
            print 'Temperature after cooling %.2f'%T

            #Check if we get a valid coloring
            assert(int(worker.get_E_from_scratch())==int(E)), \
                "Final energy must be equal to energy of returned coloring E=%.2f, Ecol=%.2f"%(int(E),int(worker.get_E_from_scratch()))
            assert (worker.X.values<=self.q).all(), \
                "Max color value must be %d"%(self.q)
            assert (worker.X.values>=1).all(), \
                "Min color value must be %d"%(1)

            print 'E=%d'%E
            print 'Ec=%d'%int(worker.get_E_from_scratch())
            self.E=E
            self.X=worker.X.values
            self.X.shape = (self.X.size,1)
            return (self.E,self.X,coloring_trace)

        elif self.schedule == 'Baseline':

            try:
                E,desc_ratio,stct_dsc_ratio = worker.walk(T,self.T_Period)
                steps += self.T_Period
                coloring_trace.append(E)
            except:
                print 'Could not finish annealing run. T hit 0.0'
                return coloring_trace

            print 'Desc ratio %.2f'%desc_ratio

            while T > self.eps :
                T = T*0.95
                try:
                    E,desc_ratio,stct_dsc_ratio = worker.walk(T,self.T_Period)
                    steps += self.T_Period
                    coloring_trace.append(E)
                except:
                    print 'Could not finish annealing run. T hit 0.0'
                    return coloring_trace
            print 'Temperature after cooling %.2f'%T

            #Check if we get a valid coloring
            assert(int(worker.get_E_from_scratch())==int(E)), \
                "Final energy must be equal to energy of returned coloring E=%.2f, Ecol=%.2f"%(int(E),int(worker.get_E_from_scratch()))
            assert (worker.X.values<=self.q).all(), \
                "Max color value must be %d"%(self.q)
            assert (worker.X.values>=1).all(), \
                "Min color value must be %d"%(1)

            print 'E=%d'%E
            print 'Ec=%d'%int(worker.get_E_from_scratch())
            self.E=E
            self.X=worker.X.values
            self.X.shape = (self.X.size,1)
            return (self.E,self.X,coloring_trace)

        elif self.schedule == 'Cauchy' :
            try:
                E,desc_ratio,stct_dsc_ratio = worker.walk(T,self.T_Period)
                steps += self.T_Period
                coloring_trace.append(E)
            except:
                print 'Could not finish annealing run. T hit 0.0'
                return coloring_trace

            print 'Desc ratio %.2f'%desc_ratio

            while  self.eps < (self.desc_thresh - desc_ratio):
                T = T * 1.5

                try:
                    E,desc_ratio,stct_dsc_ratio = worker.walk(T,self.T_Period)
                    steps += self.T_Period
                    coloring_trace.append(E)
                except:
                    print 'Could not finish annealing run. T hit 0.0'
                    return coloring_trace

            print 'Desc ratio %.2f'%desc_ratio

            while (desc_ratio - self.desc_thresh) > self.eps :
                T = T/1.5
                try:
                    E,desc_ratio,stct_dsc_ratio = worker.walk(T,self.T_Period)
                    steps += self.T_Period
                    coloring_trace.append(E)
                except:
                    print 'Could not finish annealing run. T hit 0.0'
                    return coloring_trace

            print 'Desc ratio %.2f'%desc_ratio


            while  self.eps < (self.desc_thresh - desc_ratio):
                T = T * 1.5

                try:
                    E,desc_ratio,stct_dsc_ratio = worker.walk(T,self.T_Period)
                    steps += self.T_Period
                    coloring_trace.append(E)
                except:
                    print 'Could not finish annealing run. T hit 0.0'
                    return coloring_trace
            print 'Desc ratio %.2f'%desc_ratio


            print 'Active temperature: %.2f'%T

            while T > self.eps :
                T = T/1.5
                try:
                    E,desc_ratio,stct_dsc_ratio = worker.walk(T,self.T_Period)
                    steps += self.T_Period
                    coloring_trace.append(E)
                except:
                    print 'Could not finish annealing run. T hit 0.0'
                    return coloring_trace
            print 'Temperature after cooling %.2f'%T

            #Check if we get a valid coloring
            assert(int(worker.get_E_from_scratch())==int(E)), \
                "Final energy must be equal to energy of returned coloring E=%.2f, Ecol=%.2f"%(int(E),int(worker.get_E_from_scratch()))
            assert (worker.X.values<=self.q).all(), \
                "Max color value must be %d"%(self.q)
            assert (worker.X.values>=1).all(), \
                "Min color value must be %d"%(1)

            print 'E=%d'%E
            print 'Ec=%d'%int(worker.get_E_from_scratch())
            self.E=E
            self.X=worker.X.values
            self.X.shape = (self.X.size,1)
            return (self.E,self.X,coloring_trace)
