from numpy import array, arange
import numdifftools as nd
import scipy.integrate as integrate
import re

class SimulationPotentialConverter:
    #placeholder for now
    def __init__(self):
        self.r_cut = None
        return None
    
    #loads in the potential
    def LoadPotential(self, potential, params_val):
        self.ur_ = lambda x: potential.Potential(x, params_val)
        self.dur_ = nd.Derivative(self.ur_)
        return None
    
    #creates a tabulated version of the potential (and forces)
    def TabulatePotential(self, r_table_max=20.0, dr=0.005):
        self.r = arange(0.0, r_table_max, dr)
        self.ur = self.ur_(self.r) #array([self.ur_(x) for x in self.r])
        self.dur = self.dur_(self.r) #array([self.dur_(x) for x in self.r])
        self.num_pts = len(self.r)
        return None
    
    #optional cut and shift of potential according to the forces
    def CutShiftTabulated(self, e_max=0.1, f_max=0.01, r_max=8.0, shift=False):
        #identify an acceptable cut point
        for i in range(self.num_pts-1, -1, -1):
            if abs(self.dur[i]) > f_max:
                i_cut = i + 1
                break
        self.r_cut = self.r[i_cut]
        #error check
        if self.ur[i_cut] > e_max:
            raise Exception('The potential appears to be non-trivial at the force cutoff. This could mean the potential is being truncated improperly. On likely possibility is the r_max for tabulating the potential is to short and a finite energy low force region in the potential (e.g. a plateau) has been incorrectly identified as the cutoff.')
        if self.r_cut > r_max:    
            raise Exception('The potential cutoff has exceeded the user specified maximum.')
        #cut and shift    
        if shift:
            self.ur[0:i_cut+1] = self.ur[0:i_cut+1] - self.ur[i_cut]
            self.ur[i_cut+1:] = 0.0
            self.dur[i_cut+1:] = 0.0   
        return None
    
    #make the table file with proper format
    def MakeTable(self, filename='./table.xvg'):
        #error check
        if self.r_cut is None:
            raise AttributeError('Potential cutoff was not properly set before writing table file.')
        f = open(filename, 'w')
        for i in range(self.num_pts):
            r_, ur_, dur_ = (self.r[i], self.ur[i], self.dur[i])
            f.write('{:.10e}   {:.10e} {:.10e}   {:.10e} {:.10e}   {:.10e} {:.10e}\n'.format(r_, 0.0, 0.0, 0.0, 0.0, ur_, -1.0*dur_))
        f.close()
        return None
    
    #replace the r_cut value
    def InsertGromppCutoff(self, r_buffer=0.5, filename='./grompp.mdp'):
        #error check
        if self.r_cut is None:
            raise AttributeError('Potential cutoff was not properly set before writing grompp file.')
        #read in the data first
        f = open(filename, 'r')
        gromacs_data = f.read()
        f.close()
        #substitute the cutoff
        r_cut_buff = self.r_cut + r_buffer
        gromacs_data = re.sub(r'(rlist\s*=\s*).*', (r'\g<1>%f' % r_cut_buff), gromacs_data)
        gromacs_data = re.sub(r'(rcoulomb\s*=\s*).*', (r'\g<1>%f' % r_cut_buff), gromacs_data)
        gromacs_data = re.sub(r'(rvdw\s*=\s*).*', (r'\g<1>%f' % r_cut_buff), gromacs_data)
        #overwite the old data
        f = open(filename, 'w')
        f.write(gromacs_data)
        f.close()
        return None