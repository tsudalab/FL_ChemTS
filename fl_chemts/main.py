import os, sys
from pprint import pprint
import GaussianRunPack
				
usage ='Usage; %s infile' % sys.argv[0]

try:
    infilename = sys.argv[1]
except:
    print (usage); sys.exit()

option = "opt  UV nmr"
if len(sys.argv)==3: option = sys.argv[2]
#test_sdf = GaussianRunPack_py.GaussianDFTRun('B3LYP', '6-31G',12, 'OPT energy uv homolumo deen',infilename,0)
#test_sdf = GaussianRunPack_py.GaussianDFTRun('B3LYP', '6-31G',4, 'OPT fluor',infilename,0)
test_sdf = GaussianRunPack.GaussianDFTRun('LC-BLYP', '3-21G*', 12, option,infilename,0)

#print(test_sdf.read_sdf())

outdic = test_sdf.run_gaussian()
pprint (outdic)

#gap = test_sdf.Extract_values(infilename,1,0,1,1,0)
#print (gap)
#print(Energy[-1])


