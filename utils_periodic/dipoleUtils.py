import numpy as np
import scipy as sc
import math
import EUtils
from atomUtils import Atom
from read_data_utils import read_data_charge,read_data_hybrid, read_data_dipole, initialize_atoms
from write_data_utils import write_data_hybrid, write_data_dipole
class DipoleConfig():

	def __init__(self):

		#initializing attributes
		self._eps = 8.85
		self._atoms = None
		self._N = None
		self._coeffMat = None
		self._rhsMat = None
		self._source = None
		self._dptypes = None
		self._ntypes = None
		self._qtypes = None
		self._elementMap = None
		self._lmbdax = None
		self._lmbday = None
		self._cr = None


	@property
	def eps(self):
		return self._eps
	@property
	def atoms(self):
		return self._atoms
	@property
	def N(self):
		return len(self.atoms)
	@property
	def coeffMat(self):
		return self._coeffMat
	@property
	def rhsMat(self):
		return self._rhsMat
	@property
	def source(self):
		return self._source
	@property
	def dptypes(self):
		return self._dptypes
	@property
	def ntypes(self):
		return self._ntypes
	@property
	def qtypes(self):
		return self._qtypes
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
	def cr(self):
		return self._cr







	@eps.setter
	def eps(self,data):
		self._eps = data
		return
	@atoms.setter
	def atoms(self,data):
		self._atoms = data
		return
	@N.setter
	def N(self,data):
		self._N = data
		return
	@coeffMat.setter
	def coeffMat(self,data):
		self._coeffMat = data
		return
	@rhsMat.setter
	def rhsMat(self,data):
		self._rhsMat = data
		return
	@source.setter
	def source(self,data):
		self._source = data
		return
	@source.setter
	def ntypes(self,data):
		self._ntypes = data
		return
	@dptypes.setter
	def dptypes(self,data):
		self._dptypes = data
		return
	@qtypes.setter
	def qtypes(self,data):
		self._qtypes = data
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
	@cr.setter
	def cr(self,data):
		self._cr = data
		return



	def read_source_data(self,fname,sl=6,flag=''):

		if flag == 'charge':
			atoms = read_data_charge(fname,sl,self.qtypes)
			self.source = atoms
		if flag == 'dipole':
			atoms = read_data_dipole(fname,sl,self.qtypes)
			self.source = atoms
		if flag == 'hybrid':
			atoms = read_Data_hybrid(fname,sl,self.qtypes)
			self.source = atoms

		return



	def initialize(self):

		for a in self.atoms:
			a.R = self.elementMap[a.type].R
			a.alpha = self.elementMap[a.type].alpha
			a.density = self.elementMap[a.type].den

		# initialize_atoms(self.atoms,self.elementMap)


		return None

	def read_data(self,fname,sl=6,flag=''):

		if flag == 'dipole':
			atoms = read_data_dipole(fname,sl,self.dptypes)
			self.atoms = atoms
		elif flag == 'hybrid':
			atoms = read_data_hybrid(fname,sl,self.dptypes)
			self.atoms = atoms


		return None


	def write_data(self,wfname,flag=''):

		if flag == 'dipole':
			write_data_dipole(wfname,self.atoms,self.N,self.ntypes)

		elif flag == 'hybrid':
			write_data_hybrid(wfname,self.atoms,self.N,self.ntypes)

		return None


	def compute_rij(self,atomi,atomj):

		rijx = atomj.x - atomi.x
		rijy = atomj.y - atomi.y
		rijz = atomj.z - atomi.z

		rij = np.sqrt(rijx**2 + rijy**2 + rijz**2)

		return rij, rijx, rijy, rijz

	def computeE(self,periodic=True):

		for i in range(self.N):

			atomi = self.atoms[i]
			if periodic:
				E = EUtils.computeEPeriodic(self.source,atomi,self.cr,self.lmbdax,self.lmbday)
			else:

				E = EUtils.computeE(self.source,atomi)

			atomi.ex = E[0]
			atomi.ey = E[1]
			atomi.ez = E[2]

		return



	def computeCoeff_px(self,i,j,const):
		a,b,c = 0,0,0
		atomi = self.atoms[i]
		atomj = self.atoms[j]

		if i == j:
			a = 0.5*(atomi.alpha[0][0] + atomi.alpha[0][0])
			b = 0.5*(atomi.alpha[0][1] + atomi.alpha[1][0])
			c = 0.5*(atomi.alpha[0][2] + atomi.alpha[2][0])
		else:
			rij, rijx, rijy, rijz = self.compute_rij(atomi, atomj)
			a = -(3*const*rijx*rijx/(2*rij**5)) + const/(2*rij**3)
			b = -(3*const*rijx*rijy/(2*rij**5))
			c = -(3*const*rijx*rijz/(2*rij**5))

		return a,b,c

	def computeCoeff_py(self,i,j,const):
		a,b,c = 0,0,0
		atomi = self.atoms[i]
		atomj = self.atoms[j]

		if i==j:
			a = 0.5*(atomi.alpha[0][1] + atomi.alpha[1][0])
			b = 0.5*(atomi.alpha[1][1] + atomi.alpha[1][1])
			c = 0.5*(atomi.alpha[1][2] + atomi.alpha[2][2])
		else:
			rij, rijx, rijy, rijz = self.compute_rij(atomi, atomj)
			a = -(3*const*rijx*rijy/(2*rij**5))
			b = -(3*const*rijy*rijy/(2*rij**5)) + const/(2*rij**3)
			c = -(3*const*rijy*rijz/(2*rij**5))

		return a,b,c

	def computeCoeff_pz(self,i,j,const):
		a,b,c = 0,0,0
		atomi = self.atoms[i]
		atomj = self.atoms[j]

		if i==j:
			a = 0.5*(atomi.alpha[0][2] + atomi.alpha[2][0])
			b = 0.5*(atomi.alpha[1][2] + atomi.alpha[2][1])
			c = 0.5*(atomi.alpha[2][2] + atomi.alpha[2][2])
		else:
			rij, rijx, rijy, rijz = self.compute_rij(atomi, atomj)
			a = -(3*const*rijx*rijz/(2*rij**5))
			b = -(3*const*rijy*rijz/(2*rij**5))
			c = -(3*const*rijz*rijz/(2*rij**5)) + const/(2*rij**3)

		return a,b,c


	def computeCoeffMats(self):

		mat = np.zeros((3*self.N,3*self.N))
		rhs = np.zeros(3*self.N)
		const = 1/(4*np.pi*self.eps)

		for i in range(self.N):

			atomi = self.atoms[i]

			rhs[3*i+0] = atomi.ex
			rhs[3*i+1] = atomi.ey
			rhs[3*i+2] = atomi.ez

			mat[3*i+0][3*i+0], mat[3*i+0][3*i+1], mat[3*i+0][3*i+2] = self.computeCoeff_px(i,i,const)
			mat[3*i+1][3*i+0], mat[3*i+1][3*i+1], mat[3*i+1][3*i+2]  = self.computeCoeff_py(i,i,const)
			mat[3*i+2][3*i+0], mat[3*i+2][3*i+1], mat[3*i+2][3*i+2] = self.computeCoeff_pz(i,i,const)

			for j in atomi.neighList:

				mat[3*i+0][3*j+0], mat[3*i+0][3*j+1], mat[3*i+0][3*j+2]  = self.computeCoeff_px(i,j,const)
				mat[3*i+1][3*j+0], mat[3*i+1][3*j+1], mat[3*i+1][3*j+2]  = self.computeCoeff_py(i,j,const)
				mat[3*i+2][3*j+0], mat[3*i+2][3*j+1], mat[3*i+2][3*j+2] = self.computeCoeff_pz(i,j,const)


		self.coeffMat = mat
		self.rhsMat = rhs
		return None


	def computeCoeffMatsAll(self):

		# mat, rhs = [],[]
		mat = np.zeros((3*self.N,3*self.N))
		rhs = np.zeros(3*self.N)
		const = 1/(4*np.pi*self.eps)

		for i in range(self.N):

			atomi = self.atoms[i]

			rhs[3*i+0] = atomi.ex
			rhs[3*i+1] = atomi.ey
			rhs[3*i+2] = atomi.ez


			for j in range(self.N):
				atomj = self.atoms[j]
				mat[3*i+0][3*j+0], mat[3*i+0][3*j+1], mat[3*i+0][3*j+2] = self.computeCoeff_px(i,j,const)
				mat[3*i+1][3*j+0], mat[3*i+1][3*j+1], mat[3*i+1][3*j+2] = self.computeCoeff_py(i,j,const)
				mat[3*i+2][3*j+0], mat[3*i+2][3*j+1], mat[3*i+2][3*j+2] = self.computeCoeff_pz(i,j,const)

				# print("atomi,atomj: ",i,j,mat[3*i+2][3*j+0], mat[3*i+2][3*j+1], mat[3*i+2][3*j+2])

		self.coeffMat = mat
		self.rhsMat = rhs
		return None



	def computeP(self):


		self.initialize()
		self.computeE(periodic=True)
		self.computeCoeffMatsAll()

		# f = open('p_values.txt','w');

		# for a in self.atoms:
		# 	f.write(str(a.alpha[0]) + '\n')
		# f.close()

		P, flag = sc.sparse.linalg.minres(self.coeffMat,self.rhsMat)

		for i in range(self.N):
			atomi = self.atoms[i]
			atomi.px = P[3*i + 0]
			atomi.py = P[3*i + 1]
			atomi.pz = P[3*i + 2]

		return None
