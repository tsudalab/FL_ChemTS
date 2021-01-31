import os, sys
#import GaussianRunPack
from GaussianRunPack import GaussianDFTRun				

#test_sdf = GaussianRunPack.GaussianDFTRun('B3LYP', '3-21G*',12, 'OPT energy deen nmr uv homolumo',infilename,0)
#test_sdf = GaussianRunPack.GaussianDFTRun('B3LYP', '3-21G*',1, 'OPT energy deen nmr uv homolumo',infilename,0)   
#outdic = test_sdf.run_gaussian()
ind=1
calc_sdf = GaussianDFTRun('B3LYP', '3-21G*', 1, 'OPT energy deen nmr uv homolumo', 'CheckMolopt'+str(ind)+'.sdf', 0)
outdic = calc_sdf.run_gaussian()

#print(test_sdf.read_sdf())


print(outdic)

print('uv', outdic['uv'][0][0], 'deen',  outdic['deen'], 'uv', outdic['uv'][1][0])
#gap = test_sdf.Extract_values(infilename,1,0,1,1,0)
#print (gap)
#print(Energy[-1])


