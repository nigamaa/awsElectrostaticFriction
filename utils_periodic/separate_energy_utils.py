import matplotlib.pyplot as plt
import sys
import numpy as np 
import scipy as sc
from atomUtils import Atom
from elementUtils import Element
from write_data_utils import write_data_hybrid
import os
import pickle

def read_dump(rfname,ElementMap,sl=9):
	rfile = open(rfname,'r')
	lcount = 0

	atoms = []
	for l in rfile:
		if lcount >= sl:
			tmp = l.split()
			newAtom = Atom()
			newAtom.id = int(tmp[0])
			newAtom.type = int(tmp[1])
			newAtom.x = float(tmp[2])
			newAtom.y = float(tmp[3])
			newAtom.z = float(tmp[4])
			newAtom.q = float(tmp[8])
			newAtom.px = float(tmp[9])
			newAtom.py = float(tmp[10])
			newAtom.pz = float(tmp[11])
			newAtom.diameter = 0
			newAtom.density = ElementMap[newAtom.type].den 

			atoms.append(newAtom)
		lcount += 1


	return atoms


def write_data_hybrid(wfname,atoms,N,ntypes):
	wfile = open(wfname,'w')
	wfile.write('# LAMMPS data file written by OVITO\n')
	wfile.write(str(N) + ' atoms\n')
	wfile.write(str(ntypes) + ' atom types\n\n')
	wfile.write('Atoms # hybrid dipole sphere\n\n')
	count = 1
	for a in atoms:
		tmp = [a.id, a.type, a.x, a.y, a.z, a.q, a.px, a.py, a.pz, a.diameter, a.density]
		wfile.write(' '.join(list(map(str,tmp))) + '\n')

		count += 1
	wfile.close()
	return None



def read_f_dump(rfname,sl=9,ntypes=[2,3]):


	# print('entered here')
	rfile = open(rfname,'r')
	lcount = 0

	zVec, fxVec, fyVec, fzVec = [],[],[],[]
	for l in rfile:
		if lcount == sl-1:
			tmp = l.split()
			zcol,fxcol = tmp.index('z'),tmp.index('fx')
			fycol, fzcol = tmp.index('fy'),tmp.index('fz')
			# print(zcol,fxcol,fycol,fzcol)
		if lcount >= sl:
			tmp = l.split()
			t = int(tmp[1])
			# print('type = ', t)
			if t ==2 or t == 3:

				z,fx = float(tmp[zcol-2]), float(tmp[fxcol-2]) 
				fy,fz = float(tmp[fycol-2]), float(tmp[fzcol-2])
				# print(z,fx,fy,fz)
				zVec.append(z)
				fxVec.append(fx)
				fyVec.append(fy)
				fzVec.append(fz)

		lcount += 1

	zmin = min(zVec)
	fx = sum(fxVec)
	fy = sum(fyVec)
	fz = sum(fzVec)


	return [fx,fy,fz,zmin]

def plot_sep_energies(main_directory,total_timesteps):

	# directory_total = os.fsencode(main_directory + 'config/')
	# n_files_total = 0
	# for file in os.listdir(directory_total):
	# 	filename = os.fsdecode(file)
	# 	if filename.startswith("avg_f."):
	# 		n_files_total += 1

	# directory_lj = os.fsencode(main_directory + 'config_lj/')
	# n_files_lj = 0
	# for file in os.listdir(directory_lj):
	# 	filename = os.fsdecode(file)
	# 	if filename.startswith("avg_f."):
	# 		n_files_lj += 1


	fzTotal, fzLj, fzQ, zVec = [],[],[],[]

	for i in total_timesteps:
		print(i)
		rfname = main_directory + 'config/avg_f.' + str(i) + '.dump'
		[fx,fy,fzT,zmin] = read_f_dump(rfname)
		fzTotal.append(fzT)
		zVec.append(zmin)

		rfname = main_directory + 'config_lj/avg_f.' + str(i) + '.dump'
		[fx,fy,fz,tmp] = read_f_dump(rfname)
		fzLj.append(fz)

		fzQ.append(fzT - fz)


	plt.figure(figsize=(7,5))
	plt.plot(zVec,fzTotal,label='Total')
	plt.plot(zVec,fzLj,label='LJ')
	plt.plot(zVec,fzQ,label='Q')
	plt.legend()
	plt.gca().invert_xaxis()
	plt.savefig(main_directory + 'txt_files/separate_energy.png',dpi=300)
	plt.close()


	forces = {'fzTotal' : fzTotal,
			  'fzLj' : fzLj,
			  'fzQ' : fzQ,
			  'zVec': zVec}

	# pickle.dump(forces,open('main_directory/forces.p','wb'))


	return forces






