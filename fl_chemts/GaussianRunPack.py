import os, sys, math
import re
from numpy import *
import subprocess
import AtomInfo
import shutil

class GaussianDFTRun:

    def __init__(self, functional, basis, nproc, value, in_file, error):

        self.in_file = in_file
        self.functional = functional.lower()
        self.basis = basis.lower()
        self.nproc = nproc
        self.value = value.lower()
        self.error = error

        self.mem = ''

    def Extract_ExcitedState(self, lines):
        Egrd = 0.0
        Eext = 0.0
        Found = False
        WaveLength = []
        V_OS = []
        for line in lines:
            if line.find("SCF Done:  ") >=0:
                line_SCFEnergy = re.split("\s+", line)
                Egrd = float(line_SCFEnergy[5])
            if line.find("Total Energy, E(TD-HF/TD-DFT)") >=0:
                line_TotalEnergy = line.split('=')
                Eext = float(line_TotalEnergy[1])
            if line.find("Excitation energies and oscillator strengths:") >=0:
                WaveLength = []
                V_OS = []
            if line.find("Excited State  ") >=0:
                line_StateInfo = line.split()
                WaveLength.append(float(line_StateInfo[6]))
                OS_info = line_StateInfo[8].split('=')
                V_OS.append(float(OS_info[1]))
            if line.find("-- Stationary point found.") >=0:
                Found = True
        return Found, Egrd, Eext, WaveLength, V_OS


    def Extract_Coordinate(self, lines):

        print ("Start finding coordinates...")

        Atom_index = []
        NumElement = []
        AtomicType = []
        X = []
        Y = []
        Z = []

        count = 0
        count_Station = 0
        Total_NumStation =0
        Index_Station = 0

######Counting Stationary points##############
        for line in lines:
            if line.find("-- Stationary point found.") >= 0:
                count_Station += 1

##############################################
        Total_NumStation = count_Station

        if Total_NumStation > 0:
            count_Station = 0
            print ("Total number of stationary points: ", Total_NumStation)
            for line in lines:
                if line.find("-- Stationary point found.") >= 0:
                    count_Station += 1
                    Index_Station += 1
                    print ("A stationary point is found! #", Index_Station)
                    continue
                if line.find("Standard orientation:") >=0 and count_Station > 0:
#                    print ("Standard orientaion was found")
#                    print ("Start reading coordinate")
                    count += 1
                    continue
                if count == 1:
                    Border_1 = line
#                    print (Border_1)
                    count += 1
                    continue
                if count == 2:
                    Index_1 = line
#                    print (Index_1)
                    count += 1
                    continue
                if count == 3:
                    Index_2 = line
#                    print (Index_2)
                    count += 1
                    continue
                if count == 4:
                    Border_2 = line
#                    print (Border_2)
                    count += 1
                    continue
                if count >= 5:
                    i_atom = line.split()
#                    print (i_atom)
                    if len(i_atom) == 6:
                        Atom_index.append(int(i_atom[0]))
                        NumElement.append(int(i_atom[1]))
                        AtomicType.append(int(i_atom[2]))
                        X.append(float(i_atom[3]))
                        Y.append(float(i_atom[3]))
                        Z.append(float(i_atom[3]))
                        count += 1
                        continue
                    else :
                        print ("Reading atom coordinates is finished...")
                        N = count-5
                        print ("Number of atoms: ", N)
                        count = 0
                        count_Station = 0
                    continue

        else:
            for line in lines:
                if line.find("Standard orientation:") >=0:
                    print ("Standard orientaion was found")
                    print ("Start reading coordinate")
                    count += 1
                    continue
                if count == 1:
                    Border_1 = line
#                    print (Border_1)
                    count += 1
                    continue
                if count == 2:
                    Index_1 = line
#                    print (Index_1)
                    count += 1
                    continue
                if count == 3:
                    Index_2 = line
#                    print (Index_2)
                    count += 1
                    continue
                if count == 4:
                    Border_2 = line
#                    print (Border_2)
                    count += 1
                    continue
                if count >= 5:
                    i_atom = line.split()
#                    print (i_atom)
                    if len(i_atom) == 6:
                        Atom_index.append(int(i_atom[0]))
                        NumElement.append(int(i_atom[1]))
                        AtomicType.append(int(i_atom[2]))
                        X.append(float(i_atom[3]))
                        Y.append(float(i_atom[3]))
                        Z.append(float(i_atom[3]))
                        count += 1
                        continue
                    else :
                        print ("Reading atom coordinates is finished...")
                        N = count-5
                        print ("Number of atoms: ", N)
                        count = 0
                    continue


#Translating  atomic number to element symbol
        Mol_atom = []
#        Mol_CartX = zeros(N)
#        Mol_CartY = zeros(N)
#        Mol_CartZ = zeros(N)

        for i in range(N):
            Mol_atom.append(AtomInfo.AtomicNumElec(NumElement[i]))
#            print (Mol_atom[i])

#############
        del Atom_index[:]
        del NumElement[:]
        del AtomicType[:]
        del X[:]
        del Y[:]
        del Z[:]


        return Mol_atom


    def Extract_values(self, infilename,opt, nmr,uv,energy,gap,dipole,deen,fluor,tadf):

        ifile = open(infilename,'r') #open file for reading
        lines = ifile.readlines()
        ifile.close()

        output = {}
            
        for line in lines:
            if line.find("Error termination") >=0:
                print("Gausssian is stopped due to some errors")
                return output

        if gap == 1:
            AlphaEigenVal = []
            BetaEigenVal = []
            for line in lines:
                if line.find(" basis functions, ") >=0:
                    line_StateInfo = line.split()
            #        print (line_StateInfo[0])
                    
                    NumBasisFunc = int(line_StateInfo[0])
  
                if line.find(" alpha electrons ") >=0:
                    line_StateInfo = line.split()
            #        print (line_StateInfo[0])
            #        print (line_StateInfo[3])
            
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
            #    return Alpha_gap
            else:
                Alpha_gap = 27.211*(AlphaEigenVal[NumAlphaElec]-AlphaEigenVal[NumAlphaElec-1])
                Beta_gap = 27.211*(BetaEigenVal[NumBetaElec]-BetaEigenVal[NumBetaElec-1])
  
                output["gap"] = [a.pha_gap, beta_gap]
            #    return Alpha_gap, Beta_gap
  
  
        if dipole == 1:
            Dipole_X = []
            Dipole_Y = []
            Dipole_Z = []
            Dipole_Total = []
  
            for line in lines:
                if line.find(" X= ") >=0:
                    line_StateInfo = line.split()
                    #print (line_StateInfo[1])
                    Dipole_X.append(float(line_StateInfo[1]))
                    #print (line_StateInfo[3])
                    Dipole_Y.append(float(line_StateInfo[3]))
                    #print (line_StateInfo[5])
                    Dipole_Z.append(float(line_StateInfo[5]))
                    #print (line_StateInfo[7])
                    Dipole_Total.append(float(line_StateInfo[7]))
                
            output["dipole"] = [Dipole_X[-1], Dipole_Y[-1], Dipole_Z[-1], Dipole_Total[-1]]
        #    return Dipole_X[-1], Dipole_Y[-1], Dipole_Z[-1], Dipole_Total[-1]
  
        if energy == 1:
            Energy = []
  
            for line in lines:
                if line.find("SCF Done:  ") >=0:
                    line_StateInfo = line.split()
                    #print (line_StateInfo[4])
                    Energy.append(float(line_StateInfo[4]))
            #    return Energy
            output["Energy"] = min(Energy)
  
        ##########################################################################
        # Index of Links
        ##########################################################################
        # Index = 0         : (always blank)
        # Index = opt       : Optimization of S0             [if opt==1]
        # Index = opt+nmr   : NMR chemical shift of S0       [if nmr==1]
        # Index = opt+nmr+1 : Virtical excitation (S0 -> S1) [uv or fluor or tadf]
        # Index = opt+nmr+2 : Optimization of S1             [fluor or tadf] 
        # Index = opt+nmr+3 : Optimization of T1             [tadf]
        ##########################################################################
  
        Links = self.SplitLinks(infilename)
        n = len(Links)
  
        if nmr == 1:
            Element = []
            ppm = []
  
            for line in lines:
                if line.find("Isotropic =  ") >=0:
                    line_Info = line.split()
                    #print (line_Info[1])
                    Element.append(line_Info[1])
                    #print(line_Info[4])
                    ppm.append(float(line_Info[4]))

        # calculating chemical shift for H, C, or Si

            for i in range(len(Element)):
                if Element[i] =="H" or  Element[i] =="C" or  Element[i] =="Si":
                    ppm[i] = AtomInfo.One_TMS_refer(Element[i],  self.functional, self.basis)-ppm[i]
  
        #    return Element, ppm 
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
                        #print (line_StateInfo[4])
                        Energy.append(float(line_StateInfo[4]))
        
            Mol_atom = self.Extract_Coordinate(lines)
            
#######Calculating Decomposed atoms total energy#######################
            decomposed_Energy = 0
#
            for i in range(len(Mol_atom)):    
                #print (Mol_atom[i], AtomInfo.One_Atom_Energy(Mol_atom[i], self.functional, self.basis))
                decomposed_Energy += AtomInfo.One_Atom_Energy(Mol_atom[i], self.functional, self.basis)
            print("Decomposed energy: ", decomposed_Energy)
#
#############################################################
  
        #    return deen 
            output["deen"] = min(Energy) - (decomposed_Energy)

  
        if uv == 1 or fluor == 1 or tadf == 1:
            Index = opt+nmr+1
            lines = "" if Index >= n else Links[Index].splitlines()
            Found, Egrd, Eext, WaveLength, V_OS = self.Extract_ExcitedState(lines)
            output["uv"] = [WaveLength, V_OS]
  
        if fluor == 1 or tadf == 1:
            Index = opt+nmr+2
            lines = "" if Index >= n else Links[Index].splitlines()
            S1_Found, S1_Egrd, S1_Eext, S1_WaveLength, S1_V_OS = self.Extract_ExcitedState(lines)
            output["S1 Total Energy"] = S1_Eext
            output["S1 Wavelength and Oscillator strengths"] = [S1_WaveLength, S1_V_OS]
  
        if tadf == 1:
            Index = opt+nmr+3
            lines = "" if Index >= n else Links[Index].splitlines()
            T1_Found, T1_Egrd, T1_Eext, T1_WaveLength, T1_V_OS = self.Extract_ExcitedState(lines)
            output["T1 Total Energy"] = T1_Eext
            output["T1 Wavelength and Oscillator strengths"] = [T1_WaveLength, T1_V_OS]
            TADF_Eng = 0.0
            if S1_Found and T1_Found:
               TADF_Eng = S1_Eext - T1_Eext
            output["Energy difference (S1-T1)"] = TADF_Eng
  
       # Debug print for SplitLinks
       #      for i in range(0,n):
       #          f = open(infilename+"."+str(i), 'w')
       #          for s in Sections[i]:
       #              f.write(s)
       #          f.close() 
  
        return output

    def MakeLinkTD(self, line_chk, line_method, State, Opt):

        line_method_TD = 'TD(Nstate=10, '+ State + ')'
        line_readMOGeom = 'Geom=AllCheck Guess=Read'
        s = ''
        s = s + '--Link1--'
        s = s + '\n'
        if self.mem != '':
            line_mem = '%mem='+str(self.mem)
            s = s + line_mem      
            s = s + '\n'
        if self.nproc > 1 :
            line_proc = '%nproc='+str(self.nproc)
            s = s + line_proc
            s = s + '\n'
        s = s + line_chk
        s = s + '\n'
        s = s + line_method
        s = s + '\n'
        s = s + line_method_TD
        s = s + '\n'
        s = s + line_readMOGeom
        s = s + '\n'
        if Opt == True:
            s = s + " Opt"
            s = s + '\n'
        s = s + '\n' 
        return s

    def run_gaussian(self):
        
        infilename = self.in_file

        option_line = self.value    

        options = option_line.split()

        opt = 0
        nmr = 0
        uv = 0
        energy = 0
        gap = 0
        dipole = 0
        deen = 0
        fluor = 0
        tadf = 0

        for i in range(len(options)):
#            print(options[i])
            option = options[i]
            if option == 'opt':
                opt = 1
                #print ('opt')
            elif option == 'nmr':
                nmr = 1
                #print ('nmr')
            elif option == 'uv':
                uv = 1
                #print ('uv')
            elif option == 'energy':
                energy = 1
                #print ('energy')
            elif option == 'homolumo':
                gap = 1
                #print ('HOMO/LUMO')
            elif option == 'dipole':
                dipole = 1
                #print ('dipole')
            elif option == 'deen':
                deen = 1
                #print ('Decomposition energy')
            elif option == 'fluor':
                fluor = 1
                #print ('Fluorescence')
            elif option == 'tadf':
                tadf = 1
                #print ('Thermally Activated Delayed Fluorescence')
            else:
                print('invalid option: ', option)

        PreGauInput = infilename.split('.')
    
        GauInputName = PreGauInput[0]+'.com'    

####File type of input?######################
        ReadFromchk = 0 
        ReadFromsdf = 0 

        if PreGauInput[1] == "sdf":
            ReadFromsdf = 1 
            Mol_atom, X, Y, Z, TotalCharge, SpinMulti = self.read_sdf()
        elif PreGauInput[1] == "chk":
            ReadFromchk = 1 
        else:
            print ("Invalid input file")
############################################

        line_chk = '%chk='+PreGauInput[0]

        line_method = '#p '+self.functional+'/'+self.basis

        line_comment = infilename

        ofile = open(GauInputName ,'w')

#####For reading geomerty and MO from checkpoint file###
        line_readMOGeom = 'Geom=AllCheck Guess=Read'
        line_readGeom = 'Geom=Checkpoint'
########################################################

        if self.mem != '':
            line_mem = '%mem='+str(self.mem)
            ofile.write(line_mem+'\n')

        if self.nproc > 1 :
            line_proc = '%nproc='+str(self.nproc)
            ofile.write(line_proc+'\n')
    
        ofile.write(line_chk+'\n')
        ofile.write(line_method+'\n')
        if opt == 1:
            ofile.write('Opt\n')

#####Reading Geometry and MO from Checkpoint file
        if ReadFromchk == 1:
            ofile.write(line_readMOGeom+'\n')
#################################################

        ofile.write('\n')

#####Reading Geometry from sdf file#####

        if ReadFromsdf == 1:

            ofile.write(line_comment+'\n')
            ofile.write('\n')

            ofile.write('%5d %5d \n' % (TotalCharge, SpinMulti))

            for j in range(len(Mol_atom)):
                ofile.write('%-4s % 10.5f  % 10.5f  % 10.5f \n'
                 % (Mol_atom[j],X[j], Y[j], Z[j]))

            ofile.write('\n')
#######################################

        if nmr == 1:

            ofile.write('--Link1--\n')

            if self.mem != '':
                line_mem = '%mem='+str(self.mem)
                ofile.write(line_mem+'\n')


            if self.nproc > 1 :
                line_proc = '%nproc='+str(self.nproc)
                ofile.write(line_proc+'\n')

            ofile.write(line_chk+'\n')
            ofile.write(line_method+'\n')

            line_method_nmr = 'NMR'
            ofile.write(line_method_nmr+'\n')

            ofile.write(line_readMOGeom+'\n')
            ofile.write('\n')

        if uv == 1 or fluor==1 or tadf == 1:
            sTD = self.MakeLinkTD(line_chk, line_method, "Singlet", False)
            ofile.write(sTD) 

        if fluor == 1 or tadf == 1:
            sTD = self.MakeLinkTD(line_chk, line_method, "Singlet", True)
            ofile.write(sTD) 

        if tadf == 1:
            sTD = self.MakeLinkTD(line_chk, line_method, "Triplet", True)
            ofile.write(sTD) 

        ofile.write('\n') 
        ofile.close()


############Run Gaussian##############################

        if os.path.isdir(PreGauInput[0]):
            shutil.rmtree(PreGauInput[0])
        os.mkdir(PreGauInput[0])
        shutil.move(GauInputName, PreGauInput[0])

        if ReadFromchk == 1 :
            shutil.move(infilename, PreGauInput[0]) 

        os.chdir(PreGauInput[0])

        tmp_test=open(GauInputName,'r')
        #tmp_test.close()
        #route = tmp_test.readline()

        #print (route)
        
        subprocess.call(["g16", PreGauInput[0]])

        logfile = PreGauInput[0]+'.log'

        output_dic = self.Extract_values(logfile,opt,nmr,uv,energy,gap,dipole,deen,fluor,tadf)

        os.chdir("..")

        tmp_test.close()
        return(output_dic)

    def SplitLinks(self, logfile):
        f = open(logfile, 'r')
        lines = f.readlines()
        f.close()

        Links = []
        Link = ""

        for line in lines:
            if line.find("Initial command:")>0:
                Links.append(Link)
                Link = ""

            Link = Link + line

        Links.append(Link)
        return Links

#####################To get cartesian coordinates from sdf file###############
    def read_sdf(self):

        ifile = open(self.in_file, 'r')

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
                #print (bond_info)
                bond_info = line.split()
                Bond_pair1.append(int(bond_info[0]))
                Bond_pair2.append(int(bond_info[1]))
                Bond_type.append(int(bond_info[2]))

                count +=1
                continue

            if count > N+N_Bond+3:
                mol_info = line.split()
                #print (mol_info)
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

###For debug#######################################
#
#        print (Mol_atom)
#        print (N)
#        print (len(Mol_CartX))
#
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

        return Mol_atom, Mol_CartX, Mol_CartY, Mol_CartZ, TotalCharge, SpinMulti

        
