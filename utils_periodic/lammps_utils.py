import numpy as np
import scipy as sc
import sys
import matplotlib.pyplot as plt
from charge_dipole import SourceConfig
from dipoleUtils import DipoleConfig
from elementUtils import Element
from assemblyUtils import Assembly
import math
from lammps import lammps
from ctypes import *
import ctypes
# loading MPI
from mpi4py import MPI
import numpy as np
comm = MPI.COMM_WORLD

def update_atoms(coords,charge,mu,source,tip,de=None):



	n_source_atoms = len(source.atoms)
	n_tip_atoms = len(tip.atoms)

	if de:
		n_de_atoms = len(de.atoms)
		n_total = n_source_atoms + n_tip_atoms + n_de_atoms

	else:
		n_total = n_source_atoms + n_tip_atoms

	if n_total != len(charge):
		print("*********************Error**********************")

		return None

	count = 0
	for i in range(n_source_atoms):
		source.atoms[i].x = coords[3*count]
		source.atoms[i].y = coords[3*count+1]
		source.atoms[i].z = coords[3*count+2]
		source.atoms[i].q = charge[count]
		source.atoms[i].px = mu[3*count]
		source.atoms[i].py = mu[3*count+1]
		source.atoms[i].pz = mu[3*count+2]



		count += 1

	for i in range(n_tip_atoms):
		tip.atoms[i].x = coords[3*count]
		tip.atoms[i].y = coords[3*count+1]
		tip.atoms[i].z = coords[3*count+2]
		tip.atoms[i].q = charge[count]
		tip.atoms[i].px = mu[3*count]
		tip.atoms[i].py = mu[3*count+1]
		tip.atoms[i].pz = mu[3*count+2]

		count += 1

	if de:
		for i in range(n_de_atoms):
			de.atoms[i].x = coords[3*count]
			de.atoms[i].y = coords[3*count+1]
			de.atoms[i].z = coords[3*count+2]
			de.atoms[i].q = charge[count]
			de.atoms[i].px = mu[3*count]
			de.atoms[i].py = mu[3*count+1]
			de.atoms[i].pz = mu[3*count+2]

			count += 1

		return source, tip, de

	else:
		return source, tip




def update_lammps_charge(charge,mu,source,tip,de=None):


	n_source_atoms = len(source.atoms)
	n_tip_atoms = len(tip.atoms)

	if de:
		n_de_atoms = len(de.atoms)
		n_total = n_source_atoms + n_tip_atoms + n_de_atoms
	else:
		n_total = n_source_atoms + n_tip_atoms

	if n_total != len(charge):
		print('--------------- lengths not matching----------------')
		return None

	count = 0

	for i in range(n_source_atoms):
		charge[count] = source.atoms[i].q
		count += 1

	for i in range(n_tip_atoms):
		mu[3*count]   = tip.atoms[i].px
		mu[3*count+1] = tip.atoms[i].py
		mu[3*count+2] = tip.atoms[i].pz

		count += 1

	if de:
		for i in range(n_de_atoms):
			mu[3*count] = de.atoms[i].px
			mu[3*count+1] = de.atoms[i].py
			mu[3*count+2] = de.atoms[i].pz
			count += 1

	return charge,mu


def update_neigh_list(rfname,source,sl=9):

	# reset the atoms neighbor list
	for atom in source.atoms:
		atom.neighList = []



	rfile = open(rfname,'r')
	lcount = 0
	n_source_atoms = len(source.atoms)
	nPairs,uniqueSet = 0, set()

	for l in rfile:
		if lcount >= sl:
			tmp = l.split()
			atom1,atom2 = int(tmp[0]), int(tmp[1])
			if atom1 <= n_source_atoms and atom2 <= n_source_atoms:
				source.atoms[atom1-1].neighList.append(atom2-1)
				source.atoms[atom2-1].neighList.append(atom1-1)
				uniqueSet.add(atom1)
				uniqueSet.add(atom2)
			nPairs += 1

		lcount += 1

	for atom in source.atoms:
		if len(atom.neighList) == 0:
			print('ERROR!!!!!! this atom does not have a neighbor list: ', atom.id)
			return None


	print('total number of pairs: ', nPairs, ' uniqSet Length: ', len(uniqueSet))

	return source


def isTip(atom,n_source_atoms):

	if atom > n_source_atoms:
		return True

	return False


def update_tip_neigh_list(rfname,tip,sl=9):

	n_source_atoms = len(tip.source)
	# print("number of n_source_atoms: ", n_source_atoms)
	n_tip_atoms = len(tip.atoms)

	rfile = open(rfname,'r')
	lcount = 0
	nPairs,uniqueSetE,  = 0, set()

	for atom in tip.atoms:
		atom.neighList = []
		atom.neighListE = []

	for l in rfile:
		if lcount >= sl:
			tmp = l.split()
			atom1,atom2 = int(tmp[0]),int(tmp[1])

			f1,f2 = isTip(atom1,n_source_atoms), isTip(atom2,n_source_atoms)
			tmp = [atom1,atom2,f1,f2]

			if (f1^f2):

				if f1:
					tip.atoms[atom1-n_source_atoms-1].neighListE.append(atom2-1)
					uniqueSetE.add(atom1-n_source_atoms)
					nPairs += 1
				if f2:
					tip.atoms[atom2 - n_source_atoms-1].neighListE.append(atom1-1)
					uniqueSetE.add(atom2-n_source_atoms)
					nPairs += 1

			if (f1 and f2):
				tip.atoms[atom1-n_source_atoms-1].neighList.append(atom2-n_source_atoms-1)
				tip.atoms[atom2-n_source_atoms-1].neighList.append(atom1-n_source_atoms-1)

		lcount += 1

	for atom in tip.atoms:
		if len(atom.neighList) == 0:
			print("this atom does not have neighbor list: ", atom.id)
			return None


	return tip



def plot_Q_figures(source,figname):
	Q = []
	for atom in source.atoms:
		Q.append(atom.q)

	plt.figure(figsize=(10,6))
	plt.plot(Q)
	plt.title('min = ' + str(min(Q)) + '  max = ' + str(max(Q)))
	plt.savefig(figname,dpi=300)
	plt.close()

	return

def plot_P_figures(atoms,figname):

	mux, muy, muz, mu = [],[],[],[]
	for atom in atoms:
		mux.append(atom.x)
		muy.append(atom.y)
		muz.append(atom.z)
		tmp = math.sqrt(atom.x**2 + atom.y**2 + atom.z**2)
		mu.append(tmp)

	plt.figure(figsize=(10,8))

	plt.subplot(221)
	plt.plot(mux)
	plt.title('mux')

	plt.subplot(222)
	plt.plot(muy)
	plt.title('muy')

	plt.subplot(223)
	plt.plot(muz)
	plt.title('muz')

	plt.subplot(224)
	plt.plot(mu)
	plt.title('mu')

	plt.savefig(figname,dpi=300)
	plt.close()
	return



def create_source_tip_configurations(folder,fname,v,tip_shift):
	C = Element()
	C.type, C.den = 1, 2.267

	Si = Element()
	Si.type, Si.den = 2, 2.238

	O = Element()
	O.type, O.den = 3, 2.4

	ElementMap = dict()
	ElementMap[1], ElementMap[2], ElementMap[3] = C, Si, O

	gateZ, R = 5, 0.7
	qtypes, dptypes = [1], [2,3]
	ntypes = 3

	rfname = folder + fname

	source = SourceConfig()
	source.ntypes, source.qtypes, source.dptypes = ntypes, qtypes, dptypes
	source.v, source.gateZ, source.R = v, gateZ, R
	source.read_data(rfname,sl=6,flag='hybrid')
	source.initialize(ElementMap)

	# initializing atoms in the tip
	print('------------------------ creating tip configuration-------------------------')
	tip = DipoleConfig()
	tip.ntypes, tip.qtypes, tip.dptypes = ntypes, qtypes, dptypes
	tip.read_data(rfname,sl=6,flag='hybrid')
	tip.initialize(ElementMap)
	tip.source = source.atoms

	return source,tip,ElementMap

def gather_atoms(lmp):
	coords = lmp.gather_atoms("x",1,3)
	charge = lmp.gather_atoms("q",1,1)
	mu = lmp.gather_atoms("mu",1,3)

	return coords, charge, mu
