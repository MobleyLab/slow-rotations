import torsions as tor
import matplotlib.pyplot as plt
import numpy as np

APPROX_0 = 0.000000000000000000000000000000000001


HiPen_Database = [
	[1, 61095, 'o', 'CCOc1ccc2nc(/N=C\c3ccccc3O)sc2c1',298.367, 0.884],
	[2, 77329, 'g', 'Cn1cc(Cl)c(/C=N/O)n1', 159.576, 0.718],
	[3, 79729, 'g', 'S=c1cc(-c2ccc(Cl)cc2)ss1', 244.793, 0.828],
	[4, 86442, 'g', 'CN1C(=O)C/C(=N\O)N(C)C1=O', 171.156, 0.734],
	[5, 87557, 'o', 'NNC(=O)[C@H]1C(c2ccccc2)[C@@H]1C(=O)NN', 234.259, 0.815],
	[6, 95858, 'o', 'CCO/C(O)=N/S(=O)(=O)c1ccccc1Cl', 263.702, 0.848], 
	[7, 107550, 'g', 'C/C(=N\O)c1oc(C)nc1C', 154.169, 0.709],    
	[8, 107778, 'r', 'O/N=C/C1=C(Cl)c2cc(Cl)ccc2OC1', 244.077, 0.827],
	[9, 123162, 'r', 'CC(=O)/C(=N/Nc1ccc(Cl)cc1)C(=O)c1ccccc1', 300.745, 0.886],
	[10, 133435, 'g', 'c1ccc(-c2nc3ccccc3nc2-c2ccccn2)nc1', 284.322, 0.870],
	[11, 138607, 'g', 'O=C(CC1=NO[C@H](c2ccccc2O)N1)N1CCCC1', 275.308, 0.861],
	[12, 140610, 'g', 'Cc1cc(C)c2c(=O)[nH]sc2n1', 180.232, 0.747],
	[13, 164361, 'g', 'CCON1C(=O)c2ccccc2C1=O', 191.186, 0.762],
	[14, 167648, 'g', 'Cc1ccc(COn2c(-c3ccccc3)nc3ccccc3c2=O)cc1', 342.398, 0.925],
	[15, 169358, 'g', 'CC1=Cn2c(=O)c3ccccc3c(=O)n2C1', 214.224, 0.792],
	[16, 1755198, 'o', 'CC(C)C(=O)NNC(=O)C(C)C', 172.228, 0.736],    
	[17, 1867000, 'g', 'c1ccc(-c2ccccc2-c2ccccc2)cc1', 230.31, 0.811],
	[18, 3127671, 'o', 'O=C(CSCC(=O)Nc1ccccc1)NNC(=O)c1ccccc1', 343.408, 0.926],
	[19, 4344392, 'o', 'CCOC(=O)NNC(=O)NCCCc1ccc2ccc3cccc4ccc1c2c34', 389.455, 0.966],
	[20, 4363792, 'r', 'Clc1cc(Cl)cc(/N=c2\ssnc2-c2ccccc2Cl)c1', 373.717, 0.953],
	[21, 6568023, 'g', 'O=C(NNC(=O)c1ccccc1)c1ccccc1', 240.262, 0.822],
	[22, 33381936, 'o', 'O=S(=O)(O/N=C1/CCc2ccccc21)c1ccc(Cl)cc1', 321.785, 0.907]
]

import time

start = time.time()



ligcode = "UNK"

for lig in HiPen_Database:

	ligname = lig[1]
	smiles = lig[3]
	topf_gas = f'/Users/megosato/Desktop/slow_rotations/TD_Project/Gas_Phase_Sim/NVT_Simulation_Data/nvt_{ligname}.gro'
	trajf_gas = f'/Users/megosato/Desktop/slow_rotations/TD_Project/Gas_Phase_Sim/NVT_Simulation_Data/nvt_{ligname}.xtc'

	ltf_gas = tor.LigandTorsionFinder(trajf_gas,topf_gas,ligcode,smiles)

	for t in ltf_gas.get_torsions():
		X,score,angle_min=ltf_gas.get_kde(t)
		freq_score = np.exp(score)
		for sw in np.arange(1,401,25):
			for pp in np.arange(0.001, 0.5, 0.01):

				num_peaks = ltf_gas.get_kde_num_peaks(score, smoothing_window=sw, peak_prominence=pp)

				min_max = ltf_gas.get_bounds_knn(X,num_peaks)
				min_max = sorted(min_max)
				#print(num_peaks)

				pdfs = {
					"x": X,
					"y": np.zeros(X.shape[0]),
				}

				for mm in min_max:
					gmm,x,pdf,pdf_individual,bounds = ltf_gas.get_individual_gmm(X, angle_min, mm)

					pdf_individual[1] = pdf_individual[1].flatten()
					pdfs['y'][bounds[0]:bounds[0]+len(pdf_individual[1])] = pdf_individual[1]

				# plt.plot(pdfs['x'], np.exp(score))
				# plt.plot(pdfs['x'],pdfs['y'])
				# plt.show()

				pdfs['y'][np.where(pdfs['y'] == 0)]= APPROX_0

				from scipy.special import kl_div

				print(num_peaks, sw, pp, sum(kl_div(freq_score, pdfs['y'])))


end = time.time()

print(start-end)

	# https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kl_div.html







