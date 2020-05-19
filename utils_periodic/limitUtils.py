import numpy as np
import math

class Limits():
	def __init__(self):
		self._xmin, self._xmax = math.inf, -math.inf
		self._ymin, self._ymax = math.inf, -math.inf
		self._zmin, self._zmax = math.inf, -math.inf 
		self._xcenter, self._ycenter = None, None


	@property
	def xmin(self):
		return self._xmin
	@property
	def xmax(self):
		return self._xmax
	@property
	def ymin(self):
		return self._ymin
	@property
	def ymax(self):
		return self._ymax
	@property
	def zmin(self):
		return self._zmin
	@property
	def zmax(self):
		return self._zmax
	@property
	def xcenter(self):
		return self._xcenter
	@property
	def ycenter(self):
		return self._ycenter
	
	


	@xmin.setter
	def xmin(self,data):
		self._xmin = data
		return
	@xmax.setter
	def xmax(self,data):
		self._xmax = data
		return
	@ymin.setter
	def ymin(self,data):
		self._ymin = data
		return
	@ymax.setter
	def ymax(self,data):
		self._ymax = data
		return
	@zmin.setter
	def zmin(self,data):
		self._zmin = data
		return
	@zmax.setter
	def zmax(self,data):
		self._zmax = data
		return 
	@xcenter.setter
	def xcenter(self,data):
		self._xcenter = data
		return
	@ycenter.setter
	def ycenter(self,data):
		self._ycenter = data
		return 


def compute_limits(atoms):
	limits = Limits()

	for atom in atoms:
		if atom.x < limits.xmin:
			limits.xmin = atom.x
		if atom.x > limits.xmax:
			limits.xmax = atom.x

		if atom.y < limits.ymin:
			limits.ymin = atom.y
		if atom.y > limits.ymax:
			limits.ymax = atom.y

		if atom.z < limits.zmin:
			limits.zmin = atom.z
		if atom.z > limits.zmax:
			limits.zmax = atom.z 

	limits.xcenter = 0.5*(limits.xmax + limits.xmin)
	limits.ycenter = 0.5*(limits.ymax + limits.ymin)


	return limits

	
def clip_atoms(init_atoms,lx,ly):

	limits = compute_limits(init_atoms)

	atoms = []
	for atom in init_atoms:
		cond1 = atom.x <= limits.xmin + lx
		cond2 = atom.y <= limits.ymin + ly 

		if cond1 and cond2:
			atoms.append(atom)


	return atoms


	

	
	
	
	
	
	