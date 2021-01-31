import os, sys
from pprint import pprint
import GaussianRunPack
				
usage ='Usage; %s infile' % sys.argv[0]

try:
    infile = sys.argv[1]
except:
    print (usage); sys.exit()

ifile = open(infile, 'r')
lines = ifile.readlines()
ifile.close()

option = "opt deen UV "
test = GaussianRunPack.GaussianDFTRun('LC-BLYP', '3-21G*', 20, option, infile ,0)

test.Extract_Coordinate(lines)

print (test)

