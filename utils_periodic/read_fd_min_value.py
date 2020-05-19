import numpy as np
import scipy as sc
from read_tip_force_utils import read_tip_force 

if __name__=='__main__':

	fields = ['frame','c11','c12','c13','c_sum','sum','fx','fy','fz','zmin','numAtoms']
	
	rfnames = 