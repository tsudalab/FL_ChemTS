from SDF2GauInput import GauTDDFT_ForDFT
from GaussianRunPack import GaussianDFTRun

ind = 1
calc_sdf = GaussianDFTRun('B3LYP', '3-21G*', 1, 'OPT energy deen nmr uv homolumo', 'CheckMolopt'+str(ind)+'.sdf', 0)
outdic = calc_sdf.run_gaussian()
print outdic
