import numpy as np
import scipy as sc


class Element(object):

	def __init__(self):

		self._R = 0.7
		self._alpha = np.ones((3,3))
		self._dia = 0
		self._den = 1
		self._mass = 1
		self._type = None

	@property
	def R(self):
		return self._R
	@property
	def alpha(self):
		return self._alpha
	@property
	def dia(self):
		return self._dia
	@property
	def den(self):
		return self._den
	@property
	def mass(self):
		return self._mass
	@property
	def type(self):
		return self._type
	
	
	@R.setter
	def R(self,data):
		self._R = data
		return
	@alpha.setter
	def alpha(self,data):
		self._alpha = data
		return
	@dia.setter
	def dia(self,data):
		self._dia = data
		return
	@den.setter
	def den(self,data):
		self._den = data
		return
	@type.setter
	def type(self,data):
		self._type = data
		return 