import torsions as tor
import rdkit_wrapper as rdw
import molconverter as mc
import mappings
import warnings

# gas is considered truth
topf_gas = f'/Users/megosato/Desktop/slow_rotations/TD_Project/Gas_Phase_Sim/NVT_Simulation_Data/nvt_61095.gro'
trajf_gas = f'/Users/megosato/Desktop/slow_rotations/TD_Project/Gas_Phase_Sim/NVT_Simulation_Data/nvt_61095.xtc'

# bound is inside the protein or in this case, restrained
topf_bnd = "/Users/megosato/Desktop/pseudo_binding/pseudo_trajectories/61095_6_10_20_16.gro"
trajf_bnd = "/Users/megosato/Desktop/pseudo_binding/pseudo_trajectories/61095_6_10_20_16.xtc"

smiles = "CCOc1ccc2nc(/N=C\c3ccccc3O)sc2c1"

ligtor_gas = tor.LigandTorsionFinder(trajf_gas,topf_gas,"UNK",smiles)
ligtor_bnd = tor.LigandTorsionFinder(trajf_bnd,topf_bnd,"UNK",smiles)

gas_rdmol = ligtor_gas.get_rdmol()
bnd_rdmol = ligtor_bnd.get_rdmol()

# mapping = dictionary mapping atom indices in the form {gas_idx: bnd_idx}
mapping = mappings.rd_map_mols(gas_rdmol, bnd_rdmol)

# we will determine the torsions of interest based off the 
# gas phase simulation
# bnd dihedrals will be based on the mapping conversion from gas_tor 
# to bnd_tor
gas_tor = ligtor_gas.get_torsions()
print(gas_tor)
bnd_tor = [ mappings.convert_dihedral(mapping, i) for i in gas_tor ]

for gt, bt in zip(gas_tor[3:], bnd_tor[3:]):
	# get the states based on the gas phase simulation
	gas_X, gas_scores, gas_angle_min = ligtor_gas.get_kde(gt)
	gas_num_peaks = ligtor_gas.get_kde_num_peaks(gas_scores)
	# gas_min_max contains the ranges for each of the states
	gas_min_max = ligtor_gas.get_bounds_knn(gas_X, gas_num_peaks)
	ligtor_gas.plot_dihedral_histogram(gt, angle_min=gas_angle_min, show=False, title="gas")
	ligtor_bnd.plot_dihedral_histogram(bt, angle_min=gas_angle_min, title="bnd")

	gas_angles = ligtor_gas.shift_torsion_angles(gt, angle_min=gas_angle_min)[1].flatten()
	bnd_angles = ligtor_bnd.shift_torsion_angles(bt,angle_min=gas_angle_min)[1].flatten()


