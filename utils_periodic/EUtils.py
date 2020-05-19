import numpy as np
import scipy as sc
import math
import copy
from copy import deepcopy

def computeE_j_to_i(atomi,atomj):
	const = 1/(4*np.pi*8.85)
	rijx = atomi.x - atomj.x
	rijy = atomi.y - atomj.y
	rijz = atomi.z - atomj.z

	rij = np.sqrt(rijx**2 + rijy**2 + rijz**2)
	EMag = (const*atomj.q)/(rij*rij)

	Ex = EMag*(rijx/rij)
	Ey = EMag*(rijy/rij)
	Ez = EMag*(rijz/rij)

	E = [Ex,Ey,Ez]

	return E

def computeE(source,atom):
	N = len(source)

	const = 1/(4*np.pi*8.85)

	Ex, Ey, Ez = 0,0,0

	for i in range(N):
		atomi = source[i]
		ex,ey,ez = computeE_j_to_i(atom,atomi)
		Ex += ex
		Ey += ey
		Ez += ez

	E = [Ex, Ey, Ez]

	return E




def computePeriodicAtoms(atom,atomi,cr,lmbdax,lmbday):

	neighAtoms = []

	nx = int(cr/lmbdax)
	ny = int(cr/lmbday)

	for r in range(-nx,nx+1):
		for c in range(-ny,ny+1):
			xshift = (r)*lmbdax
			yshift = (c)*lmbday

			nAtom = copy.copy(atomi)
			nAtom.x += xshift
			nAtom.y += yshift
			neighAtoms.append(nAtom)

	return neighAtoms






def computeEPeriodic(source,atom,cr,lmbdax,lmbday):

	N = len(source)
	const = 1/(4*np.pi*8.85)

	Ex,Ey,Ez = 0,0,0

	for i in range(N):
		atomi = source[i]
		periodicAtoms = computePeriodicAtoms(atom,atomi,cr,lmbdax,lmbday)
		for atomj in periodicAtoms:
			e = computeE_j_to_i(atom,atomj)
			Ex += e[0]
			Ey += e[1]
			Ez += e[2]

	E = [Ex,Ey,Ez]
	# print(Ex,Ey,Ez)
	return E
