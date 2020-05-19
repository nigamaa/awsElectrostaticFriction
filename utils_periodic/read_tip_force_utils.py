import numpy as np
import scipy as sc
import sys
sys.path.append('./../utils/')
sys.path.append('./../datafiles/Lammps/')
# sys.path.append('./../../storage/examples/HOPG-SiO2-dipole/4/txt_files/')
import matplotlib.pyplot as plt

def read_tip_force(rfname,sl,fields,sampleHeight):

	data = dict()

	for f in fields: 
		data[f] = []


	rfile = open(rfname,'r')

	for l in rfile:
		tmp = l.split()
		tmp2 = list(map(float,tmp))

		for i in range(len(tmp2)):
			data[fields[i]].append(tmp2[i])

	data['fnorm'] = []
	for i in range(len(data['fx'])):
		data['zmin'][i] = data['zmin'][i]  - sampleHeight
		tmp = np.sqrt(data['fx'][i]**2 + data['fy'][i]**2 + data['fz'][i]**2)
		data['fnorm'].append(tmp)
	return data

def plot_tip_force(fields,sampleHeight,rfname,figname,approach,retract,flag='CNT'):

	# rfname = '/home/nigamaa/storage/examples/' + flag + '-SiO2-dipole/' + str(folder) + '/txt_files/tip_force.txt'

	data = read_tip_force(rfname,0,fields,sampleHeight)
	print('min value of force: ', min(data['fz'][:]))
	
	plt.plot(data['zmin'][0:approach],1*np.array(data['fz'][0:approach]),label='approach')
	if retract == True :
		plt.plot(data['zmin'][approach+1:],1*np.array(data['fz'][approach+1:]), label='retract')
	
	# plt.plot(data['zmin'][:],1*np.array(data['fz'][:]))

	plt.plot(np.linspace(-2.5,15.5,100), np.zeros((100,)),'k--')
	# plt.ylim((-11,25))
	plt.xlim((-4,15.5))
	plt.xlabel('Tip-Sample seperation (Angstroms)')
	plt.ylabel('Force on the tip (nN)')
	plt.gca().invert_xaxis()
	plt.legend()
	# plt.savefig('/home/nigamaa/storage/examples/' + flag + '-SiO2-dipole/' + str(folder) + '/txt_files/tip_force.png')
	plt.savefig(figname)
	plt.show()


def plot_tip_force_overlay(fields,sampleHeight,rfnames,labels,figname,approach,retract):

	for i in range(len(rfnames)):
		# folder = folderVec[i]
		label = labels[i]
		# rfname = '/home/nigamaa/storage/examples/' + flag + '-SiO2-dipole/' + str(folder) + '/txt_files/tip_force.txt'
		rfname = rfnames[i]
		data = read_tip_force(rfname,0,fields,sampleHeight)
		print('min value of force: ', min(data['fz'][:]))
		stepsize = 1
		# plt.plot(data['zmin'][0:endframe],1*np.array(data['fz'][0:endframe]), label=label)
		plt.plot(data['zmin'][0:approach[i]],1*np.array(data['fz'][0:approach[i]]),label=label + ' approach')
		if retract[i] == True:
			plt.plot(data['zmin'][approach[i]:],1*np.array(data['fz'][approach[i]:]),label=label + ' retract')
	plt.plot(np.linspace(-2.5,15.5,100), np.zeros((100,)),'k.')
	plt.xlim((-4,15.5))
	# plt.ylim((-6,5))
	plt.xlabel('Tip-Sample seperation (Angstroms)')
	plt.ylabel('Force on the tip (nN)')
	plt.gca().invert_xaxis()
	plt.legend()
	# plt.savefig('/home/nigamaa/storage/examples/' + flag + '-SiO2-dipole/summary/tip_force.png')
	plt.savefig(figname)
	plt.show()


if __name__ == '__main__':

	fields = ['frame','c11','c12','c13','c_sum','sum','fx','fy','fz','zmin','numAtoms']
	sampleHeight = 57.0822
	folderVec = ['0K-50V']	
	
	for folderNum in folderVec:
		# rfname = '/home/nigamaa/storage/examples/' + flag + '-SiO2-dipole/' + str(folder) + '/txt_files/tip_force.txt'
		# figname = '/home/nigamaa/storage/examples/' + flag + '-SiO2-dipole/' + str(folder) + '/txt_files/tip_force.png'

		rfname = '/home/nigamaa/storage/examples/CNT-SiO2-dipole/' + str(folderNum) + '/speed1p0/txt_files/tip_force.txt'
		figname = '/home/nigamaa/storage/examples/CNT-SiO2-dipole/' + str(folderNum) + '/speed1p0/txt_files/tip_force.png'
		approach = 8500
		retract = False
		# plot_tip_force(fields,sampleHeight,rfname,figname,approach,retract,flag='CNT')


	rfnames, approach, retract = [],[],[]
	rfnames.append('/home/nigamaa/storage/examples/CNT-SiO2-dipole/0K-15V/speed0p2/txt_files/tip_force.txt')
	rfnames.append('/home/nigamaa/storage/examples/CNT-SiO2-dipole/0K-15V/speed2p0/txt_files/tip_force.txt')
	# rfnames.append('/home/nigamaa/storage/examples/CNT-SiO2-dipole/0K-10V/speed0p2/txt_files/tip_force.txt')
	# rfnames.append('/home/nigamaa/storage/examples/CNT-SiO2-dipole/0K-15V/speed0p2/txt_files/tip_force.txt')
	# rfnames.append('/home/nigamaa/storage/examples/CNT-SiO2-dipole/0K-25V/speed1p0/txt_files/tip_force.txt')
	# rfnames.append('/home/nigamaa/storage/examples/CNT-SiO2-dipole/0K-50V/speed1p0/txt_files/tip_force.txt')
	figname = '/home/nigamaa/storage/examples/CNT-SiO2-dipole/summary/tip_force_0k_15v_speed_0p2_2p0.png'
	approach.append(8500), retract.append(True)
	approach.append(750), retract.append(True)
	# approach.append(8500)
	# approach.append(8500)
	# approach.append(1500)
	# approach.append(1500)
	plot_tip_force_overlay(fields,sampleHeight,rfnames,['0K-15V-0P2','0K-15V-2P0'],figname,approach,retract)
	
	
