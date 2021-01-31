import os, sys, math
from numpy import *
import subprocess
import AtomInfo

class GaussianDFTRun:

	def __init__(self, functional, basis, nproc, value, sdf_file, error):

		self.sdf_file = sdf_file
		self.functional = functional.lower()
		self.basis = basis.lower()
		self.nproc = nproc
		self.value = value.lower()
		self.error = error

		self.decomposed_Energy = 0.0

	def Extract_values(self, infilename,nmr,uv,energy,gap,dipole,deen):

		ifile = open(infilename,'r') #open file for reading
		lines = ifile.readlines()
		ifile.close()

		output = {}
		
		if gap == 1:
			AlphaEigenVal = []
			BetaEigenVal = []
			for line in lines:
				if line.find(" basis functions, ") >=0:
					line_StateInfo = line.split()
			#		print (line_StateInfo[0])
					
					NumBasisFunc = int(line_StateInfo[0])

				if line.find(" alpha electrons ") >=0:
					line_StateInfo = line.split()
			#		print (line_StateInfo[0])
			#		print (line_StateInfo[3])
			
					NumAlphaElec = int(line_StateInfo[0])
					NumBetaElec = int(line_StateInfo[3])
				
				if line.find("Alpha  occ. eigenvalues --") >=0:

					if len(AlphaEigenVal) == NumBasisFunc:
						AlphaEigenVal=[]
					
					line_removed = line.replace("Alpha  occ. eigenvalues --", " ")
					line_StateInfo = line_removed.split()
					for i in range(len(line_StateInfo)):
						AlphaEigenVal.append(float(line_StateInfo[i]))

				if line.find("Alpha virt. eigenvalues --") >=0:
					line_removed = line.replace("Alpha virt. eigenvalues --", " ")
					line_StateInfo = line_removed.split()
					for i in range(len(line_StateInfo)):
						AlphaEigenVal.append(float(line_StateInfo[i]))


				if line.find("Beta  occ. eigenvalues --") >=0:

					if len(BetaEigenVal) == NumBasisFunc:
						BetaEigenVal=[]
					
					line_removed = line.replace("Beta  occ. eigenvalues --", " ")
					line_StateInfo = line_removed.split()
					for i in range(len(line_StateInfo)):
						BetaEigenVal.append(float(line_StateInfo[i]))

				if line.find("Beta virt. eigenvalues --") >=0:
					line_removed = line.replace("Beta virt. eigenvalues --", " ")
					line_StateInfo = line_removed.split()
					for i in range(len(line_StateInfo)):
						BetaEigenVal.append(float(line_StateInfo[i]))

			if BetaEigenVal == []:
				Alpha_gap = 27.211*(AlphaEigenVal[NumAlphaElec]-AlphaEigenVal[NumAlphaElec-1])
				output["gap"] = Alpha_gap
			#	return Alpha_gap
			else:
				Alpha_gap = 27.211*(AlphaEigenVal[NumAlphaElec]-AlphaEigenVal[NumAlphaElec-1])
				Beta_gap = 27.211*(BetaEigenVal[NumBetaElec]-BetaEigenVal[NumBetaElec-1])

				output["gap"] = [a.pha_gap, beta_gap]
			#	return Alpha_gap, Beta_gap


		if dipole == 1:
			Dipole_X = []
			Dipole_Y = []
			Dipole_Z = []
			Dipole_Total = []

			for line in lines:
				if line.find(" X= ") >=0:
					line_StateInfo = line.split()
					print (line_StateInfo[1])
					Dipole_X.append(float(line_StateInfo[1]))
					print (line_StateInfo[3])
					Dipole_Y.append(float(line_StateInfo[3]))
					print (line_StateInfo[5])
					Dipole_Z.append(float(line_StateInfo[5]))
					print (line_StateInfo[7])
					Dipole_Total.append(float(line_StateInfo[7]))
				
			output["dipole"] = [Dipole_X[-1], Dipole_Y[-1], Dipole_Z[-1], Dipole_Total[-1]]
		#	return Dipole_X[-1], Dipole_Y[-1], Dipole_Z[-1], Dipole_Total[-1]

		if energy == 1:
			Energy = []

			for line in lines:
        			if line.find("SCF Done:  ") >=0:
                			line_StateInfo = line.split()
                			print (line_StateInfo[4])
					Energy.append(float(line_StateInfo[4]))
                        #	return Energy
			output["Energy"] = Energy[-1]

		if uv == 1:
			WaveLength = []
			V_OS = []

			for line in lines:
        			if line.find("Excited State  ") >=0:
                			line_StateInfo = line.split()
                			print (line_StateInfo[6])
                			WaveLength.append(float(line_StateInfo[6]))
                			OS_info = line_StateInfo[8].split('=')
                			print(OS_info[1])
                			V_OS.append(float(OS_info[1]))

		#	return WaveLength, V_OS
			output["uv"] = [WaveLength, V_OS]

		if nmr == 1:
			Element = []
			ppm = []

			for line in lines:
				if line.find("Isotropic =  ") >=0:
					line_Info = line.split()
					print (line_Info[1])
					Element.append(line_Info[1])
					print(line_Info[4])
					ppm.append(float(line_Info[4]))

		#	return Element, ppm 
			output["nmr"] = [Element, ppm]


		#print (output)


		if deen == 1:
			try:
				Energy = Energy
			except NameError:
				Energy = []

				for line in lines:
        				if line.find("SCF Done:  ") >=0:
                				line_StateInfo = line.split()
                				print (line_StateInfo[4])
						Energy.append(float(line_StateInfo[4]))

		#	return deen 
			output["deen"] = Energy[-1] - (self.decomposed_Energy)


		return output


	def run_gaussian(self):
		
		infilename = self.sdf_file

		option_line = self.value	

		options = option_line.split()

		opt = 0
		nmr = 0
		uv = 0
		energy = 0
		gap = 0
		dipole = 0
		deen = 0

		for i in range(len(options)):
#			print(options[i])
			option = options[i]
			if option == 'opt':
				opt = 1
				print ('opt')
			elif option == 'nmr':
				nmr = 1
				print ('nmr')
			elif option == 'uv':
				uv = 1
				print ('uv')
			elif option == 'energy':
				energy = 1
				print ('energy')
			elif option == 'homolumo':
				gap = 1
				print ('HOMO/LUMO')
			elif option == 'dipole':
				dipole = 1
				print ('dipole')
			elif option == 'deen':
				deen = 1
				print ('Decomposition energy')
			else:
				print('invalid option: ', option)
		
		Mol_atom, X, Y, Z, TotalCharge, SpinMulti = self.read_sdf()

		PreGauInput = infilename.split('.')
	
		GauInputName = PreGauInput[0]+'.com'	

		line_chk = '%chk='+PreGauInput[0]

		line_method = '#'+self.functional+'/'+self.basis

		line_comment = infilename

		ofile = open(GauInputName ,'w')

		if self.nproc > 1 :
			line_proc = '%nproc='+str(self.nproc)
			ofile.write(line_proc)
			ofile.write('\n')
	
		ofile.write(line_chk)
		ofile.write('\n')
		ofile.write(line_method)
		ofile.write('\n')
		if opt == 1:
			ofile.write('Opt')
			ofile.write('\n')
		ofile.write('\n')

		ofile.write(line_comment)
		ofile.write('\n')
		ofile.write('\n')

		ofile.write('%5d %5d \n' % (TotalCharge, SpinMulti))

		for j in range(len(Mol_atom)):
       			ofile.write('%-4s % 10.5f  % 10.5f  % 10.5f \n'
			 % (Mol_atom[j],X[j], Y[j], Z[j]))

		ofile.write('\n')

		if nmr == 1:

			line_readMOGeom = 'Geom=AllCheck Guess=Read'

			ofile.write('--Link1--')
        		ofile.write('\n')

        		if self.nproc > 1 :
				line_proc = '%nproc='+str(self.nproc)
	                	ofile.write(line_proc)
				ofile.write('\n')

			ofile.write(line_chk)
			ofile.write('\n')
			ofile.write(line_method)
			ofile.write('\n')

			line_method_nmr = 'NMR'
			ofile.write(line_method_nmr)
			ofile.write('\n')

			ofile.write(line_readMOGeom)
			ofile.write('\n')
			ofile.write('\n')

		if uv == 1:
			line_readMOGeom = 'Geom=AllCheck Guess=Read'

			ofile.write('--Link1--')
        		ofile.write('\n')

        		if self.nproc > 1 :
				line_proc = '%nproc='+str(self.nproc)
	                	ofile.write(line_proc)
				ofile.write('\n')

			ofile.write(line_chk)
			ofile.write('\n')
			ofile.write(line_method)
			ofile.write('\n')

			line_method_TD = 'TD(Nstate=20, Singlet)'
			ofile.write(line_method_TD)
			ofile.write('\n')

			ofile.write(line_readMOGeom)
			ofile.write('\n')
			ofile.write('\n') 

		ofile.write('\n') 
		ofile.close()


############Run Gaussian##############################

                #subprocess.call(["mkdir", PreGauInput[0]])
		#subprocess.call(["mv", GauInputName,PreGauInput[0]])
		#os.chdir(PreGauInput[0])

		tmp_test=open(GauInputName,'r')

		route = tmp_test.readline()

		print (route)

		subprocess.call(["g16", PreGauInput[0]])

		logfile = PreGauInput[0]+'.log'

		output_dic = self.Extract_values(logfile,nmr,uv,energy,gap,dipole,deen)

		#os.chdir("..")


		return(output_dic)


#####################To get cartesian coordinates from sdf file###############
	def read_sdf(self):

		ifile = open(self.sdf_file, 'r')

		count = 0

		X = []
		Y = []
		Z = []
		element_symbol = []

		Bond_pair1 = []
		Bond_pair2 = []
		Bond_type = []

		TotalCharge = 0 
		CHG_atom = []
		CHG = []

		for line in ifile:
			if count == 0:
				Header1 = line

				count += 1
				continue

			if count == 1:
				Header2 = line

				count += 1
				continue

			if count == 2:
				Header3 = line

				count += 1
				continue

			if count == 3:
				a = line.split()
				N = int(a[0])
				N_Bond = int(a[1])

				count += 1
				continue

			if 3 < count <= N+4:
				i_atom = line.split()
				if len(i_atom) != 0:
					X.append(float(i_atom[0]))
					Y.append(float(i_atom[1]))
					Z.append(float(i_atom[2]))
					element_symbol.append(i_atom[3])
	
				count += 1
				continue

			if N+4 < count <= N+N_Bond+3 :
				bond_info = line.split()
				print (bond_info)
				bond_info = line.split()
				Bond_pair1.append(int(bond_info[0]))
				Bond_pair2.append(int(bond_info[1]))
				Bond_type.append(int(bond_info[2]))

				count +=1
				continue

			if count > N+N_Bond+3:
				mol_info = line.split()
				print (mol_info)
				if (mol_info[0] == "M"):
					if (mol_info[1] == "END"):
						break
					if (mol_info[1] == "CHG"):
						Num_CHGInfo = int(mol_info[2])
						for k in range(Num_CHGInfo):
							CHG_atom.append(int(mol_info[3+2*k]))
							CHG.append(int(mol_info[4+2*k]))
							TotalCharge += int(mol_info[4+2*k])
				else:
					print("The sdf file is invalid!")
					sys.exit()

				count +=1

#Copy to array of numpy###########################

		Mol_atom = []
		Mol_CartX = zeros(N)
		Mol_CartY = zeros(N)
		Mol_CartZ = zeros(N)

		CHG_atom = array(CHG_atom)
		CHG = array(CHG)

		for j in range(N):
			Mol_CartX[j] = X[j] 
			Mol_CartY[j] = Y[j] 
			Mol_CartZ[j] = Z[j] 
			Mol_atom.append(element_symbol[j])

	#del element_symbol[N:TotalStep]

		del element_symbol[:]
		del X[:]
		del Y[:]
		del Z[:]

##################################################

		print (Mol_atom)
		print (N)
		print (len(Mol_CartX))

		print('Reading the sdf file has finished')

###Calculating the total number of electrons#################
		TotalNum_electron = 0 

		for j in range(N):
			if (len(CHG_atom) != 0):
				Judge = CHG_atom-j
				if (any(Judge) == 0):
					TotalNum_electron += AtomInfo.AtomicNumElec(Mol_atom[j])-CHG[where(Judge == 0)]
				else:
					TotalNum_electron += AtomInfo.AtomicNumElec(Mol_atom[j])
			else:
				TotalNum_electron += AtomInfo.AtomicNumElec(Mol_atom[j])

		print('Total number of electron: %7d ' %  (TotalNum_electron))

		if (TotalNum_electron%2==0):
			print ("This system is a closed shell!")
			SpinMulti = 1
		else:
			print ("This system is a open shell!")
			SpinMulti = 2

#######Calculating Decomposed atoms total energy#######################

		for i in range(N):	
			print (Mol_atom[i], AtomInfo.One_Atom_Energy(Mol_atom[i], self.functional, self.basis))
			self.decomposed_Energy += AtomInfo.One_Atom_Energy(Mol_atom[i], self.functional, self.basis)
		print("Decomposed energy: ", self.decomposed_Energy)

#############################################################

		return Mol_atom, Mol_CartX, Mol_CartY, Mol_CartZ, TotalCharge, SpinMulti

		
