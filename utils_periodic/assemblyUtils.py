import numpy as np
import scipy as sc
import sys
sys.path.append('./utils/')
from atomUtils import Atom 
from elementUtils import Element
from read_data_utils import read_data, read_data_charge, read_data_hybrid, read_data_dipole, initialize_atoms
from write_data_utils import write_data, write_data_charge, write_data_dipole, write_data_hybrid


	

class Assembly(object):

	def __init__(self):

		self._atoms = []
		self._ntypes = None
		self._N = None


	@property
	def atoms(self):
		return self._atoms
	@property
	def ntypes(self):
		return self._ntypes
	
	@property
	def N(self):
		return len(self._atoms)
		
	
	
	@atoms.setter
	def atoms(self,data):
		self._atoms = data
		return
	@ntypes.setter
	def ntypes(self,data):
		self._ntypes = data
		return 


	def append(self,atoms,shift):


		for a in atoms:
			a.x += shift[0]
			a.y += shift[1]
			a.z += shift[2]

			self._atoms.append(a)

		return None


	def initialize(self,elementMap):

		initialize_atoms(self.atoms,elementMap)

		return None

	def write_data(self,wfname,flag=''):


		print('writing data')
		if flag == 'charge':
			write_data_charge(wfname,self.atoms,self.N,self.ntypes)

		elif flag == 'dipole':
			write_data_dipole(wfname,self.atoms,self.N,self.ntypes)

		elif flag == 'hybrid':
			write_data_hybrid(wfname,self.atoms,self.N,self.ntypes)

		else:

			write_data(wfname,self.atoms,self.N,self.ntypes)

		return None