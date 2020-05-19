import numpy as np
import scipy as sc
from atomUtils import Atom 

def read_data(fname,sl,types):

	rfile = open(fname,'r')

	count = 0
	atoms = []
	for l in rfile:
		if count >= sl:
			newAtom = Atom()
			tmp = l.split()
			newAtom.id = int(tmp[0])
			newAtom.type = int(tmp[1])
			newAtom.x = float(tmp[2])
			newAtom.y = float(tmp[3])
			newAtom.z = float(tmp[4])

			if newAtom.type in types:

				atoms.append(newAtom)

		count += 1

	return atoms

def read_data_charge(fname,sl,types):

	rfile = open(fname,'r')
	count = 0
	atoms = []

	for l in rfile:
		if count >= sl:
			newAtom = Atom()
			tmp = l.split()
			newAtom.id = int(tmp[0])
			newAtom.type = int(tmp[1])
			newAtom.q = float(tmp[2])
			newAtom.x = float(tmp[3])
			newAtom.y = float(tmp[4])
			newAtom.z = float(tmp[5])

			if newAtom.type in types:
				atoms.append(newAtom)

		count += 1

	return atoms


def read_data_dipole(fname,sl,types):

	rfile = open(fname,'r')
	count = 0
	atoms = []

	for l in rfile:
		if count >= sl:
			newAtom = Atom()
			tmp = l.split()
			newAtom.id = int(tmp[0])
			newAtom.type = int(tmp[1])
			newAtom.q = float(tmp[2])
			newAtom.x = float(tmp[3])
			newAtom.y = float(tmp[4])
			newAtom.z = float(tmp[5])
			newAtom.px = float(tmp[6])
			newAtom.py = float(tmp[7])
			newAtom.pz = float(tmp[8])

			if newAtom.type in types:
				atoms.append(newAtom)

		count += 1

	return atoms


def read_data_hybrid(fname,sl,types):

	rfile = open(fname,'r')
	count = 0
	atoms = []

	for l in rfile:
		if count >= sl:
			newAtom = Atom()
			tmp = l.split()
			newAtom.id = int(tmp[0])
			newAtom.type = int(tmp[1])
			newAtom.x = float(tmp[2])
			newAtom.y = float(tmp[3])
			newAtom.z = float(tmp[4])
			newAtom.q = float(tmp[5])
			newAtom.px = float(tmp[6])
			newAtom.py = float(tmp[7])
			newAtom.pz = float(tmp[8])
			newAtom.diameter = float(tmp[9])
			newAtom.density = float(tmp[10])

			if newAtom.type in types:

				atoms.append(newAtom)

		count += 1

	return atoms

def initialize_atoms(atoms,elementMap):

	for a in atoms:
		a.R = elementMap[a.type].R
		a.alpha = elementMap[a.type].alpha
		a.diameter = elementMap[a.type].dia 
		a.density = elementMap[a.type].den 
		

	return atoms