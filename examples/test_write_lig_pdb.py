import MDAnalysis as mda
import pdb_writer as pdbw
from pathlib import Path




def export_pdb_from_traj(mda_u, opdb, sel="all"):
    atms = mda_u.select_atoms(sel)
    atms.write(opdb, frames=mda_u.trajectory[[0,]])



base_path = Path('/Users/megosato/Desktop/mpro/md_traj/')

ligname = 'P2203'
system = f'{ligname}_0A'

# bound is inside the protein or in this case, restrained
topf_bnd = base_path / f'complex/{system}/npt.gro'
trajf_bnd = base_path / f'complex/{system}/prod.xtc'

u = mda.Universe(str(topf_bnd), str(trajf_bnd))


opdb = "/Users/megosato/Desktop/test_lig.pdb"

export_pdb_from_traj(u, opdb, sel="resname UNK")
opdb2 = "/Users/megosato/Desktop/test_lig2.pdb"

pdbw.rename_lig_pdb_atoms(opdb, opdb2)

