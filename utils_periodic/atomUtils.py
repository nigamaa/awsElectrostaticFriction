#!/usr/bin/env python

import numpy as np
#import matplotlib.pyplot as plt
from copy import deepcopy
import math
from operator import itemgetter
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy.special import erf

class Atom():
	def __init__(self):
		self._x, self._y, self._z = None, None, None
		self._id, self._type, self._q = None, None, 0.0
		self._v, self._mu = None, None
		self._px, self._py, self._pz = 0.0, 0.0, 0.0
		self._alpha = np.ones((3,3))
		self._ex, self._ey, self._ez = None, None, None
		self._diameter, self._density = 0, 1
		self._R = None
		self._element = None
		self._neighList = []
		self._neighListE = []
		self._fx, self._fy, self._fz = None, None, None
		self._f = None
		self._vx, self._vy, self._vz = None, None, None


	@property
	def x(self):
		return self._x
	@property
	def y(self):
		return self._y
	@property
	def z(self):
		return self._z
	@property
	def px(self):
		return self._px
	@property
	def py(self):
		return self._py
	@property
	def pz(self):
		return self._pz
	@property
	def alpha(self):
		return self._alpha
	@property
	def id(self):
		return self._id
	@property
	def type(self):
		return self._type
	@property
	def q(self):
		return self._q
	@property
	def v(self):
		return self._v
	@property
	def mu(self):
		return self._mu
	@property
	def ex(self):
		return self._ex
	@property
	def ey(self):
		return self._ey
	@property
	def ez(self):
		return self._ez
	@property
	def diameter(self):
		return self._diameter
	@property
	def density(self):
		return self._density
	@property
	def R(self):
		return self._R
	@property
	def neighList(self):
		return self._neighList
	@property
	def neighListE(self):
		return self._neighListE
	@property
	def fx(self):
		return self._fx
	@property
	def fy(self):
		return self._fy
	@property
	def fz(self):
		return self._fz
	@property
	def f(self):
		return self._f

	@property
	def vx(self):
		return self._vx
	@property
	def vy(self):
		return self._vy
	@property
	def vz(self):
		return self._vz














	# attribute setters
	@x.setter
	def x(self,data):
		self._x = data
		return
	@y.setter
	def y(self,data):
		self._y = data
		return
	@z.setter
	def z(self,data):
		self._z = data
		return
	@px.setter
	def px(self,data):
		self._px = data
		return
	@py.setter
	def py(self,data):
		self._py = data
		return
	@pz.setter
	def pz(self,data):
		self._pz = data
		return
	@alpha.setter
	def alpha(self,data):
		self._alpha = data
		return
	@id.setter
	def id(self,data):
		self._id = data
		return
	@type.setter
	def type(self,data):
		self._type = data
		return
	@q.setter
	def q(self,data):
		self._q = data
		return
	@v.setter
	def v(self,data):
		self._v = data
		return
	@mu.setter
	def mu(self,data):
		self._mu = data
		return
	@ex.setter
	def ex(self,data):
		self._ex = data
		return
	@ey.setter
	def ey(self,data):
		self._ey = data
		return
	@ez.setter
	def ez(self,data):
		self._ez = data
		return
	@diameter.setter
	def diameter(self,data):
		self._diameter = data
		return
	@density.setter
	def density(self,data):
		self._density = data
		return
	@R.setter
	def R(self,data):
		self._R = data
		return
	@neighList.setter
	def neighList(self,data):
		self._neighList = data
		return

	@neighListE.setter
	def neighListE(self,data):
		self._neighListE = data
		return

	@fx.setter
	def fx(self,data):
		self._fx = data
		return
	@fy.setter
	def fy(self,data):
		self._fy = data
		return
	@fz.setter
	def fz(self,data):
		self._fz = data
		return
	@f.setter
	def f(self,data):
		self._f = data
		return

	@vx.setter
	def vx(self,data):
		self._vx = data
		return
	@vy.setter
	def vy(self,data):
		self._vy = data
		return
	@vz.setter
	def vz(self,data):
		self._vz = data
		return


	def findDist(self,atom):
		rijx = atom.x - self.x
		rijy = atom.y - self.y
		rijz = atom.z - self.z

		rij = np.sqrt(rijx**2 + rijy**2 + rijz**2)

		return [rijx,rijy,rijz,rij]
