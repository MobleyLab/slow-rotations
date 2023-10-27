import torsions as tor
import rdkit_wrapper as rdw
import molconverter as mc
import mappings
import warnings
import compare

# gas is considered truth
topf_gas = f'/Users/megosato/Desktop/slow_rotations/TD_Project/Gas_Phase_Sim/NVT_Simulation_Data/nvt_61095.gro'
trajf_gas = f'/Users/megosato/Desktop/slow_rotations/TD_Project/Gas_Phase_Sim/NVT_Simulation_Data/nvt_61095.xtc'

# bound is inside the protein or in this case, restrained
topf_bnd = "/Users/megosato/Desktop/pseudo_binding/pseudo_trajectories/61095_6_10_20_16.gro"
trajf_bnd = "/Users/megosato/Desktop/pseudo_binding/pseudo_trajectories/61095_6_10_20_16.xtc"

smiles = "CCOc1ccc2nc(/N=C\c3ccccc3O)sc2c1"

ligcode = "UNK"

ligtor_gas = tor.LigandTorsionFinder(trajf_gas,topf_gas,ligcode,smiles)
ligtor_bnd = tor.LigandTorsionFinder(trajf_bnd,topf_bnd,ligcode,smiles)

lc = compare.LigandComparator(ligtor_gas, ligtor_bnd)

for g in lc.get_gas_torsions():
	lc.make_compare_torsions_img(g, show=True)


	