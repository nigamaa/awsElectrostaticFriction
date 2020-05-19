import numpy as np
from atomUtils import Atom 
from elementUtils import Element 

def write_data(wfname,atoms,N,ntypes):

	wfile = open(wfname,'w')
	wfile.write('# LAMMPS data file written by OVITO\n')
	wfile.write(str(N) + ' atoms\n')
	wfile.write(str(ntypes) + ' atom types\n\n')
	wfile.write('Atoms # atomic\n\n')
	count = 1
	for a in atoms:
		tmp = [count, a.type, a.x, a.y, a.z]
		wfile.write(' '.join(list(map(str,tmp))) + '\n')
		count += 1
	wfile.close()
	return None

def write_data_charge(wfname,atoms,N,ntypes):

	wfile = open(wfname,'w')
	wfile.write('# LAMMPS data file written by OVITO\n')
	wfile.write(str(N) + ' atoms\n')
	wfile.write(str(ntypes) + ' atom types\n\n')
	wfile.write('Atoms # charge\n\n')
	count = 1
	for a in atoms:
		tmp = [count, a.type, a.q, a.x, a.y, a.z]
		wfile.write(' '.join(list(map(str,tmp))) + '\n')
		count += 1
	wfile.close()
	return None

def write_data_dipole(wfname,atoms,N,ntypes):

	wfile = open(wfname,'w')
	wfile.write('# LAMMPS data file written by OVITO\n')
	wfile.write(str(N) + ' atoms\n')
	wfile.write(str(ntypes) + ' atom types\n\n')
	wfile.write('Atoms # dipole\n\n')
	count = 1
	for a in atoms:
		tmp = [count, a.type, a.q, a.x, a.y, a.z, a.px, a.py, a.pz]
		wfile.write(' '.join(list(map(str,tmp))) + '\n')
		count += 1
	wfile.close()
	return None


def write_data_hybrid(wfname,atoms,N,ntypes):

	wfile = open(wfname,'w')
	wfile.write('# LAMMPS data file written by OVITO\n')
	wfile.write(str(N) + ' atoms\n')
	wfile.write(str(ntypes) + ' atom types\n\n')
	wfile.write('Atoms # hybrid dipole sphere\n\n')
	count = 1
	for a in atoms:
		tmp = [count, a.type, a.x, a.y, a.z, a.q, a.px, a.py, a.pz, a.diameter, a.density]
		wfile.write(' '.join(list(map(str,tmp))) + '\n')

		count += 1
	wfile.close()
	return None