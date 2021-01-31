import os, sys, math
from numpy import *
import AtomInfo

usage ='Usage; %s infile' % sys.argv[0]

try:
    infile = sys.argv[1]
except:
    print (usage); sys.exit()

ifile = open(infile, 'r')
lines = ifile.readlines()
ifile.close()

Atom_index = []
NumElement = []
AtomicType = []
X = []
Y = []
Z = []

count = 0
NumStation = 0

for line in lines:
    if line.find("-- Stationary point found.") >= 0:
        print ("A stationary point is found!")
        NumStation += 1
        continue
    if line.find("Standard orientation:") >=0 and NumStation > 0:
        print ("Standard orientaion was found")
        print ("Start reading coordinate")
        count += 1
        continue
    if count == 1:
        Border_1 = line
        print (Border_1)
        count += 1
        continue
    if count == 2:
        Index_1 = line
        print (Index_1)
        count += 1
        continue
    if count == 3:
        Index_2 = line
        print (Index_2)
        count += 1
        continue
    if count == 4:
        Border_2 = line
        print (Border_2)
        count += 1
        continue
    if count >= 5:
        i_atom = line.split()
        print (i_atom)
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
            NumStation = 0
        continue

#Translating  atomic number to element symbol
Mol_atom = []
#        Mol_CartX = zeros(N)
#        Mol_CartY = zeros(N)
#        Mol_CartZ = zeros(N)

for i in range(N):
#    print(NumElement[i])
    Mol_atom.append(AtomInfo.AtomicNumElec(NumElement[i]))
    print (Mol_atom[i])




