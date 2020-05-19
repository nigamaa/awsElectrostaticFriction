import numpy as np
import scipy as sc
from scipy.sparse.linalg import spsolve
from copy import deepcopy
import math
import pickle
import copy
from limitUtils import Limits, compute_limits
from scipy.special import erf
from atomUtils import Atom
from read_data_utils import read_data, read_data_dipole, read_data_charge, read_data_hybrid
from read_data_utils import initialize_atoms
from write_data_utils import write_data, write_data_charge, write_data_dipole, write_data_hybrid

def empty_copy(obj):
	class Empty(obj._class__):
		def __init__(self): pass
	newCopy = Empty()
	newCopy.__class__ = obj.__class__
	return newCopy

class SourceConfig():
	def __init__(self):

		#initializing attributes
		self._eps = 8.85
		self._atoms, self._N = None, None
		self._gateZ = None			# gate in Z for mirror atoms
		self._v = None				# bias of each of the atom
		self._coeffMat = None		# coeffcient matrix of linear system of equation
		self._rhsMat = None			# rhs matrix of the linear system of equations
		self._R = None				# radius of the gaussian cloud
		self._ntypes = None
		self._qtypes = None 		# charge types
		self._dptypes = None 		# polarization types
		self._elementMap = None 	# elementMap
		self._lmbdax = None 		# periodicity in x direction
		self._lmbday = None 		# periodicity in y direction
		self._cntDia = None 		# diameter of the cnt in the source file
		self._cr = None 		# cut off radius factor



	@property
	def eps(self):
		return self._eps
	@property
	def N(self):
		return len(self._atoms)
	@property
	def atoms(self):
		return self._atoms
	@property
	def gateZ(self):
		return self._gateZ
	@property
	def v(self):
		return self._v
	@property
	def coeffMat(self):
		return self._coeffMat
	@property
	def rhsMat(self):
		return self._rhsMat
	@property
	def R(self):
		return self._R
	@property
	def ntypes(self):
		return self._ntypes
	@property
	def qtypes(self):
		return self._qtypes
	@property
	def dptypes(self):
		return self._dptypes
	@property
	def elementMap(self):
		return self._elementMap
	@property
	def lmbdax(self):
		return self._lmbdax
	@property
	def lmbday(self):
		return self._lmbday
	@property
	def cntDia(self):
		return self._cntDia
	@property
	def cr(self):
		return self._cr




	@R.setter
	def R(self,data):
		self._R = data
		return

	@eps.setter
	def eps(self,data):
		self._eps = data
		return
	@N.setter
	def N(self,data):
		self._N = data
		return
	@atoms.setter
	def atoms(self,data):
		self._atoms = data
		return
	@gateZ.setter
	def gateZ(self,data):
		self._gateZ = data
		return
	@v.setter
	def v(self,data):
		self._v = data
		return
	@coeffMat.setter
	def coeffMat(self,data):
		self._coeffMat = data
		return
	@rhsMat.setter
	def rhsMat(self,data):
		self._rhsMat = data
		return
	@ntypes.setter
	def ntypes(self,data):
		self._ntypes = data
		return
	@qtypes.setter
	def qtypes(self,data):
		self._qtypes = data
		return
	@dptypes.setter
	def dptypes(self,data):
		self._dptypes = data
		return
	@elementMap.setter
	def elementMap(self,data):
		self._elementMap = data
		return
	@lmbdax.setter
	def lmbdax(self,data):
		self._lmbdax = data
		return
	@lmbday.setter
	def lmbday(self,data):
		self._lmbday = data
		return
	@cntDia.setter
	def cntDia(self,data):
		self._cntDia = data
		return
	@cr.setter
	def cr(self,data):
		self._cr = data
		return


	def __copy__(self):
		newCopy = empty_copy(self)
		newCopy.x = self.x
		newCopy.y = self.y
		newCopy.z = self.z
		newCopy.q = self.q
		return newCopy

	def initialize_bias(self):

		for a in self.atoms:
			a.v = self.v
		return None

	def initialize_gateZ(self):

		for a in self.atoms:
			a.z += self.gateZ
		return

	def deinitialize_gateZ(self):

		for a in self.atoms:
			a.z -= self.gateZ
		return None

	def initialize(self):

		initialize_atoms(self.atoms,self.elementMap)

		self.initialize_bias()
		self.initialize_gateZ()

		return None

	def read_data(self,fname,sl=6,flag=''):

		if flag == 'charge':
			atoms = read_data_charge(fname,sl,self.qtypes)
			self.atoms = atoms

		elif flag == 'dipole':
			atoms = read_data_dipole(fname,sl,self.qtypes)
			self.atoms = atoms

		elif flag == 'hybrid':
			atoms = read_data_hybrid(fname,sl,self.qtypes)
			self.atoms = atoms

		else:
			atoms = read_data(fname,sl,self.qtypes)
			self.atoms = atoms



		return None


	def write_data(self,wfname,flag=''):

		if flag == 'charge':
			write_data_charge(wfname,self.atoms,self.N,self.ntypes)

		elif flag == 'dipole':
			write_data_dipole(wfname,self.atoms,self.N,self.ntypes)

		elif flag == 'hybrid':
			write_data_hybrid(wfname,self.atoms,self.N,self.ntypes)

		else:
			write_data(wfname,self.atoms,self.N,self.ntypes)


		return



	def computeAB(self,atomi,atomj):

		const  = 1/(4*np.pi*self.eps)
		a = 1/(2*self.R*self.R)

		A,B = 0,0

		if atomi.x == atomj.x and atomi.y == atomj.y:
			Ai = ((np.sqrt(2/np.pi)) * const)/self.R
			riIi = 2*atomi.z
			riIi = riIi/erf(riIi*np.sqrt(a))
			BiIi = const/riIi
			A = Ai
			B = BiIi

		else:

			rij = np.sqrt((atomi.x-atomj.x)**2 + (atomi.y-atomj.y)**2 + (atomi.z-atomj.z)**2)
			rij = rij/erf(rij*np.sqrt(a))
			Bij = const/rij

			riIj = np.sqrt((atomi.x-atomj.x)**2 + (atomi.y-atomj.y)**2 + (atomi.z+atomj.z)**2)
			riIj = riIj/erf(riIj*np.sqrt(a))
			BiIj = const/riIj
			A = Bij
			B = BiIj

		return A,B


	def computeCoeffMats(self):

		mat = np.zeros((self.N,self.N))
		rhs = np.zeros(self.N)

		for i in range(self.N):

			atomi = self.atoms[i]
			rhs[i] = atomi.v

			A,B = self.computeAB(atomi,atomj)
			mat[i][i] = A - B

			for j in atomi.neighList:
				atomj = self.atoms[j]
				A,B = self.computeAB(atomi,atomj)
				mat[i][j] = A - B

		self.coeffMat = mat
		self.rhsMat = rhs
		return


	def computePeriodicAtoms(self,i,j):
		atomi = copy.copy(self.atoms[i])
		atomj = copy.copy(self.atoms[j])
		neighAtoms = []

		nx = int(self.cr/self.lmbdax)
		ny = int(self.cr/self.lmbday)

		for r in range(-nx,nx+1):
			for c in range(-ny,ny+1):
				xshift = (r)*self.lmbdax
				yshift = (c)*self.lmbday

				nAtom = copy.copy(atomj)
				nAtom.x = atomj.x+xshift
				nAtom.y = atomj.y+yshift
				neighAtoms.append(nAtom)

		return neighAtoms


	def computeCoeffMatsPeriodic(self):

		# print("computeCoeffMatsPeriodic Function")

		mat = np.zeros((self.N,self.N))
		rhs = np.zeros(self.N)

		for i in range(self.N):

			atomi = self.atoms[i]
			rhs[i] = (atomi.v)

			# for each atom in the neighbor list, get its periodic image
			# compute A,B for periodic images
			# add its effect on to the matrix
			for j in atomi.neighList:
				atomj = self.atoms[j]
				periodicAtoms = self.computePeriodicAtoms(i,j)

				for atomk in periodicAtoms:
					A,B = self.computeAB(atomi,atomk)
					mat[i][j] += A - B

		# updating the coeffMat and rhsMat
		self.coeffMat = mat
		self.rhsMat = rhs

		return





	def computeQ(self,all=False):
		# print("entered computeQ function")
		self.initialize()

		if all:
			self.computeCoeffMats()
		else:

			self.computeCoeffMatsPeriodic()

		Q, flag = sc.sparse.linalg.minres(self.coeffMat,self.rhsMat)


		for i in range(self.N):
			atomi = self.atoms[i]
			atomi.q = (0.001/1.602)*Q[i]

		# print('computed charges')

		self.deinitialize_gateZ()

		return None
