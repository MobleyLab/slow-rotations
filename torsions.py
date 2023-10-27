from pathlib import Path
import tempfile
import rdkit
import math
import warnings
from openeye import oechem
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD, GRO, XTC
from MDAnalysis.analysis.dihedrals import Dihedral
from MDAnalysis.analysis import distances

import utils
import mappings
import rdkit_wrapper as rdw
import molconverter as mc
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from scipy.signal import find_peaks
from sklearn.neighbors import KernelDensity
from sklearn.cluster import KMeans

class BadTorsionError():
    pass

class TorsionFinder():

    def __init__(self, trajf:str, topf: str):
        self.trajf = trajf
        self.topf = topf
        self.mda_universe = mda.Universe(self.topf, self.trajf)

    def get_torsions():
        raise utils.NotImplementedError

    @staticmethod
    def get_angle_shift_point(angles, num_bins=60):
        angles_of_bins = np.linspace(-180,180,num_bins+1)
        hist = np.histogram(angles, angles_of_bins)[0]
        
        max_streak = 0
        current_streak = 0
        first_zero_index = -1
        first_zero_index_of_max_streak = 0
        min_freq = min(hist)
        for i, n in enumerate(hist):
            if n == min_freq:
                if current_streak == 0:
                    first_zero_index = i
                current_streak += 1
            else:
                if current_streak > max_streak:
                    max_streak = current_streak
                    first_zero_index_of_max_streak = first_zero_index
                current_streak = 0
        if current_streak > max_streak:
            first_zero_index_of_max_streak = first_zero_index

        angle_min = angles_of_bins[first_zero_index_of_max_streak]
        
        return angle_min

    def get_torsion_angles(self, torsion:[int,int,int,int]):
        # ags = atom groups
        # mda_universe.atoms is indexed by a tuple of lists of atom
        # indices; each torsion must be a list
        ags = [mda.AtomGroup(self.mda_universe.atoms[(torsion)])]
        Run = Dihedral(ags).run()
        shape = (Run.results.angles.shape)
        angles = Run.results.angles
        return angles

    def shift_torsion_angles(self, torsion:[int,int,int,int], num_bins=60, angle_min=None):
        angles = np.array(self.get_torsion_angles(torsion))
        if not angle_min:
            angle_min = TorsionFinder.get_angle_shift_point(angles, num_bins)

        shifted_angles = np.array(angles)
        for i in range(len(shifted_angles)):
            if shifted_angles[i] < (angle_min + 360//num_bins):
                shifted_angles[i] = shifted_angles[i] + 360
        return angle_min, shifted_angles

    def get_kde(self, torsion, num_bins=60, angle_min=None):

        angle_min, X = self.shift_torsion_angles(torsion, num_bins, angle_min)
        X = np.sort(X.flatten())
        X = X.reshape(-1, 1)
        kde = KernelDensity(kernel='gaussian', bandwidth=6).fit(X)
        score = kde.score_samples(X)
        return X, score, angle_min

    def get_kde_num_peaks(self, scores, smoothing_window=160, peak_prominence=0.008):
        # moving average kde scores
        moving_avg = []
        for i in range(len(scores) - smoothing_window):
            moving_avg.append(np.mean(scores[i:i+smoothing_window]))

        # peaks from moving average
        peaks = find_peaks(moving_avg, prominence=peak_prominence)[0]
        num_peaks = peaks.shape[0]
        return num_peaks

    def get_bounds_knn(self, X, num_components, tolerance=10):
        data = np.column_stack((X[:,0], np.zeros(X.shape)))
        kmeans = KMeans(n_clusters=num_components, random_state=0)
        kmeans.fit(data)

        y_kmeans = kmeans.predict(data)

        min_max = []
        
        for c in range(kmeans.n_clusters):
            c_mask = np.where(y_kmeans == c)
            c_angles = X[:,0].transpose()[c_mask]            
            c_center = kmeans.cluster_centers_[c][0]
                        
            centroid_index = None
            min_diff = 999999999
            min_i = None
            for i,a in enumerate(c_angles):
                diff = abs(c_center - a)
                if diff < min_diff:
                    min_diff = diff
                    min_i = i
            centroid_index = min_i
            
            min_i = 0
            max_i = len(c_angles) - 1
            i = centroid_index
            while i > 0:
                if abs(c_angles[i] - c_angles[i-1]) > tolerance:
                    min_i = i
                    break
                i -= 1
            
            for i in range(centroid_index,len(c_angles) - 1):
                if abs(c_angles[i+1] - c_angles[i]) > tolerance:
                    max_i = i
                    break
                    
            min_max.append((c_angles[min_i],c_angles[max_i]))
        return min_max


    def get_individual_gmm(self, X, angle_min, min_max):
        num_components = 1
        min_bnd = min_max[0]
        max_bnd = min_max[1]
        min_max_range = max_bnd - min_bnd

        single_state_X = list()

        print(X)

        for point in X:
            if point[0] >= min_bnd and point[0] <= max_bnd:
                single_state_X.append(point)

        gmm, x, pdf, pdf_individual = self.get_gmm(single_state_X, num_components, angle_min)

        frac = len(single_state_X)/len(X)
        pdf_individual = pdf_individual*frac
        pdf_individual = np.array([x.flatten(), pdf_individual.flatten()])
        
        return gmm,x,pdf,pdf_individual


    def get_gmm(self, X, num_components, angle_min):
        flat_X = np.array(X).flatten()
        gmm = GaussianMixture(n_components=num_components).fit(np.array(X).reshape(-1,1))
        x = np.linspace(min(flat_X), max(flat_X), len(X))
        logprob = gmm.score_samples(x.reshape(-1, 1))
        responsibilities = gmm.predict_proba(x.reshape(-1, 1))
        pdf = np.exp(logprob)
        pdf_individual = responsibilities * pdf[:, np.newaxis]


        return gmm, x, pdf, pdf_individual

    def get_bounds_gmm():
        pass

    def transition_counter(self, angles, range_of_states):
        which_state = []
        for num in angles:
            for i in range(0, len(range_of_states)):
                if range_of_states[i][0] <= num <= range_of_states[i][1]:
                    which_state.append(i+1)

        # 2 entries for every state [transitions in, transitions out]
        transition_count = np.zeros((len(range_of_states),2),int)
        for i in range(1, len(which_state)):
            if which_state[i] != which_state[i-1]:
                
                # transition into state i
                transition_count[which_state[i]-1,0] += 1
                
                # transition out of state i-1
                transition_count[which_state[i-1]-1,1] += 1
        return transition_count

    def check_transitions(self, transition_matrix, min_transitions):
        if len(transition_matrix) == 1:
            return
        for i,in_out in enumerate(transition_matrix):
            if in_out[0] < min_transitions:
                warnings.warn(f"State {i} has fewer that {min_transitions} in")
            if in_out[1] < min_transitions:
                warnings.warn(f"State {i} has fewer that {min_transitions} out")


    ######################
    # PLOTTING FUNCTIONS #
    ######################

    def plot_dihedral_scatter(self, torsion, title=None, show=True, savepath=None):
        f, ax = plt.subplots()
        X = self.get_torsion_angles(torsion)
        ax.scatter(np.arange(len(X)), X)
        ax.set_ylabel("Dihedral Angle (˚)")
        ax.set_xlabel("Frame")
        ax.set_ylim([-180,180])
        if title: 
            ax.set_title(title)
        if show:
            plt.show()
        if savepath:
            plt.savefig(savepath)
        return f,ax

    def plot_dihedral_histogram(self, torsion, angles=[], angle_min=None, pdf_individual=[], ax=None, title=None, num_bins=60, alpha=0.5, show=True, color=None, save_path=None):
        if not ax:
            f, ax = plt.subplots()

        if angle_min == None:
            angle_min, angles = self.shift_torsion_angles(torsion, num_bins, angle_min=angle_min)

        X = np.array(angles).flatten()
        kwargs = dict()
        if color != None:
            kwargs['color']=color
        ax.hist(X, bins=num_bins, alpha=alpha, density=True, stacked=True, **kwargs)
        ax.set_xlim([angle_min, angle_min+360])
        ax.set_xlabel("Dihedral Angle (˚)")
        ax.set_ylabel("Frequency")

        pdf_individual_sort = sorted(pdf_individual, key=len)
        colors = ['red', 'orange', 'green', 'blue', 'purple', 'brown']

        for i,pdf in enumerate(pdf_individual_sort):
            ax.plot(pdf[0], pdf[1], color=colors[i])

        if title:
            ax.set_title(title)
        if show:
            plt.show()
        if save_path:
            plt.savefig(save_path, dpi=500)

    def plot_transition_matrix(self, transition_matrix, ax=None):

        if not ax:
            f,ax = plt.subplots()

        df = pd.DataFrame(data=transition_matrix, columns=['transitions in', 'transitions out'])
        table = ax.table(
            cellText=df.values, 
            colLabels=df.columns, 
            rowLabels=[f"s{i}" for i in range(len(df))],
            loc='center',
        )
        cells = table.properties()["celld"]
        for i in range(0, len(df)+1):
            cells[i, 0].set_text_props(ha="center")
            cells[i, 1].set_text_props(ha="center")

        ax.axis('off')


class ProteinTorsionFinder(TorsionFinder):

    def __init__(self, trajf: str, topf: str, ligcode: str):
        TorsionFinder.__init__(self, trajf, topf)
        self.ligcode = ligcode
        self.aa_only = self.mda_universe.select_atoms("protein and (name N or name CA or name C or name O or name CB)")
        self.ligand = self.mda_universe.select_atoms(f'resname {self.ligcode}')

    def get_binding_residues(self, A_cutoff: float):
        ''' a_cutoff: Angstrom cutoff from binding mode
        '''
        
        indices_list = []
        #for each ts (frame) in trajectory
        for ts in self.mda_universe.trajectory:
            #Calculate all possible distances between a reference set (ligand atoms) of postions and configureation set (aa atoms)
            dist_arr = distances.distance_array(self.ligand.positions,self.aa_only.positions)

            #Find where in distance array values are less than A_cutoff (angstroms) and extend to list
            indices = np.where(dist_arr < A_cutoff)
            indices_list.extend(indices[1].tolist())

        #only unique values
        unique_indices = list(set(indices_list))

        aa_interacting_dict = {}
        #For each number in list, find corresponding atom and store the number and name of its corresponding residue in dict
        for i in unique_indices:
            aa_of_atom = self.aa_only.atoms[i]
            aa_interacting_dict[int(aa_of_atom.resid)] = aa_of_atom.resname

        aa_interacting_list = sorted(list(aa_interacting_dict))
        return aa_interacting_list

    def get_chi1_torsions(self, resids):
        ''' 
        '''
        # MO TODO: Should this be consistent with get_torsions and return a 
        # list of a list of indices

        ags = [res.chi1_selection() for res in self.aa_only.residues[resids]]
        torsions = []

        for i,ag in enumerate(ags):
            if ag:
                torsions.append([ a.index for a in ag ])
            # else the residue is a "GLY" or "ALA"
            else:
                warnings.warn(f"skipping resid {resids[i]} (ALA or GLY); no chi1 torsion to select")

        return torsions

    # should I write a function to get all protein torsions?


    def get_intersection(self, list_of_sets: list):
        ''' returns a set that is an intersection of all the sets
            in the list
        '''
        if len(list_of_sets) == 0:
            return set()

        s1 = list_of_sets[0]

        for s in list_of_sets[1:]:
            s1 = s1.intersection(s)

        return s1


    def get_chi_x_residues(self, x, sel=None, a_cutoff=None):
        ''' returns a list of residues that have chi_x angles in
            the protein
        '''

        if x > 8:
            raise BadTorsionError

        if sel:
            residues_ags = self.mda_universe.select_atoms(sel)

        if a_cutoff:
            residues_ags = self.get_residues_ags(self.get_binding_residues(a_cutoff))

        residues = residues_ags.residues

        ags = self.get_chi_x_ags(x, residues)

        ags_res_set_list = []
        for ag in ags:
            ags_res_set_list.append(set([a.residue for a in ag]))

        return sorted(list(self.get_intersection(ags_res_set_list)))

    def get_residues_ags(self, resgrp):
        ''' 
        '''
        ags = self.mda_universe.select_atoms("")

        for r in resgrp:
            if type(r) == int:
                ags += self.mda_universe.select_atoms(f"resid {r}")
            elif type(r) == mda.Residue:
                ags += self.mda_universe.select_atoms(f"resid {r.resid}")
        return ags


    def get_chi_x_ags(self, x, resgrp):
        ''' returns a list of atom groups inolved in chi x
            based on the residue list
        '''

        x_atom_sel = [
            "name N", "name CA", 
            "name CB", "name CG CG1", 
            "name CD CD1 OD1 ND1 SD", "name NE OE1 CE", 
            "name CZ NZ", "name NH1"
        ]

        ags = list()

        x0base = x - 1
        for sel in x_atom_sel[x0base:x0base+4]:
            ags.append(resgrp.atoms.select_atoms(sel))

        return ags


    def get_chi_x_aid(self, x, resgrp):
        ''' returns a list of lists of atom indices (aid)
        '''
        ag_list = self.get_chi_x_ags(x, resgrp)

        chi_aid = list()

        for ag in ag_list:
            chi_aid.append(np.array([a.index for a in ag]))

        chi_aid_np = np.array(chi_aid)
        return chi_aid_np.transpose()


    def get_chi_x_torsions(self, x, sel=None, a_cutoff=None):
        ''' returns a list of chi x torsions atom indices in the 
            protein given the selection (sel)
        '''
        # list of the residues (not resid) that have a chiX torsion
        chi_x_reslst = self.get_chi_x_residues(x, sel, a_cutoff)
        chi_x_resgrp = mda.ResidueGroup(chi_x_reslst)

        resnames = list()
        for r in chi_x_resgrp:
            resnames.append(r.resname)

        return self.get_chi_x_aid(x, chi_x_resgrp)



    def save_traj_sel(self, sel, frames, save_path):
        ''' saves out the atoms selected in sel as a new file for the 
            purpose of this program meant to save out a single 
            amino acid as a pdb file

            sel:        selection string using MDAnalysis atom selection language
            frames:     2 membered tuple with the start and end (not inclusive)
                        frame to save out [start, end)
            save_path:  path to save 
        '''
        ag = self.mda_universe.select_atoms(sel)

        print(ag)

        start_incl, end_nincl = frames

        with mda.Writer(save_path, ag.n_atoms) as w:
            for ts in self.mda_universe.trajectory[start_incl:end_nincl]:
                w.write(ag)



    def highlight_dihedral(self, dihedral, save_path=None):
        ''' highlights the dihedral of a protein sidechain
            * exports the protein residue as a single AA pdb
            * loads in the single protein into rdkit with and
              without Hs
            * compares the 2 to get atom orderings via MCS
        '''
        sel_resid = self.get_residue_from_torsion(dihedral)
        min_res_aidx = min([a.index for a in sel_resid.atoms])

        with tempfile.NamedTemporaryFile(suffix='.pdb') as aapdb:

            print(f'resid {sel_resid.resid}')

            self.save_traj_sel(f'resid {sel_resid.resid}', (0,1), aapdb.name)

            rdmol_wH = rdw.load_rdmol_from_file(aapdb.name, removeHs=False)
            rdmol_woH = rdw.load_rdmol_from_file(aapdb.name, removeHs=True)

            mapping = mappings.map_mols(rdmol_wH, rdmol_woH)

            adj_dih = [idx - min_res_aidx for idx in dihedral]

            woH_adj_dih = mappings.convert_dihedral(mapping, adj_dih)

            rdw.highlight_dihedral(rdmol_woH, woH_adj_dih, save_path)


    def get_residue_from_torsion(self, torsion):
        sel_atm_in_dih = self.mda_universe.select_atoms(f"index {torsion[0]}")
        return sel_atm_in_dih[0].residue

    def determine_chi_x(self, torsion):
        raise utils.NotImplementedError


    def make_aa_torsions_img(self, torsion, show=False):
        
        d1,d2,d3,d4 = tuple(torsion)

        sel_a_in_dih = self.mda_universe.select_atoms(f"index {torsion[0]}")
        sel_resid = sel_a_in_dih[0].residue


        f,ax = plt.subplots(1, 3, figsize=(25, 6.25))
        sup_title = f"{sel_resid.resname} {sel_resid.resid} ({d1},{d2},{d3},{d4})"
        f.suptitle(sup_title,fontsize=60)
        f.tight_layout(pad=3.5)

        X, scores, angle_min = self.get_kde(torsion)

        

        num_peaks = self.get_kde_num_peaks(scores, smoothing_window=100, peak_prominence=0.008)
        
        min_max = self.get_bounds_knn(X, num_peaks)
        angles = self.shift_torsion_angles(torsion, angle_min=angle_min)[1].flatten()

        transition_matrix = self.transition_counter(angles, min_max)

        with tempfile.NamedTemporaryFile(suffix='.png') as highlightpng:

            self.highlight_dihedral(torsion, save_path=highlightpng)
            img = np.asarray(Image.open(highlightpng.name))
            ax[0].imshow(img)
            ax[0].axis('off')

        pdf_individual = []
        for mm in min_max:
            gmm,x,pdf,pdfi = self.get_individual_gmm(X, angle_min, mm)
            pdf_individual.append(pdfi)

        states_list = [f"s{i}" for i in range(num_peaks)]
        self.plot_dihedral_histogram(torsion, ax=ax[1], show=False, pdf_individual=pdf_individual)
        ax[1].legend(states_list)
        self.plot_transition_matrix(transition_matrix, ax=ax[2])







class LigandTorsionFinder(TorsionFinder):
    ''' This class requires openeye

        MO: problem... need to be able to distiguish if 2 dihedrals are the same
            in 2 representations of same molecule
    '''
    def __init__(self, trajf: str, topf: str, ligcode: str, smiles: str, complex=False, molfile=None):
        ''' molfile: with atom numbering that is equivalent to the trajectory mol atom numbering
        '''

        TorsionFinder.__init__(self, trajf, topf)
        warnings.warn("Ensure atom number in molecule structure file is the same as trajectory molecule atom numbering")

        self.smiles = smiles
        self.molfile = molfile
        self.ligcode = ligcode


        # this will need to be changed 
        with tempfile.NamedTemporaryFile(suffix='.pdb') as temppdb:
            selection = f"resname {ligcode}"
            self.export_pdb_from_traj(temppdb.name, sel=selection)
            self.rdmol = rdw.assign_bond_order_from_smiles(smiles, temppdb.name)
        
        self.rdmol = rdw.sanitize_rdmol(self.rdmol)
        self.oemol = mc.get_oemol_from_rdmol(self.rdmol)

        # using rdkit
        # load in smiles
        # load in mol
        # force the atom numbering from molfile onto smiles
        # use openff to convert rdkit to oemol

    def get_oemol(self):
        return self.oemol

    def get_rdmol(self):
        return self.rdmol

    def convert_ligidx_to_sysidx(self, idx):
        ''' when ligand is not the only molecule in the system
            converts the ligand indices starting at 0 to indices
            that are relevant to the system
        '''
        selection = f"resname {self.ligcode}"
        ags = self.mda_universe.select_atoms(selection)
        indices = [ a.index for a in ags ]
        return min(indices) + idx 


    def export_pdb_from_traj(self, opdb, sel="all"):
        atms = self.mda_universe.select_atoms(sel)
        atms.write(opdb, frames=self.mda_universe.trajectory[[0,]])


    def _get_torsion(self, bond):
        # credit to: Travis Dabbous
        # MO TODO: fix this 
        # this is NOT deterministic, same molecule can return DIFFERENT torsions each time
        pos_1 = []
        pos_2 = []
        pos_3 = []
        pos_4 = []
        
        bgn_atom = bond.GetBgn()
        pos_2 = bond.GetBgnIdx()
        end_atom = bond.GetEnd()
        pos_3 = bond.GetEndIdx()
        #find the neighbors of the atoms from their lists
        nbors_pos2 = []
        
        for atom in bgn_atom.GetAtoms():
            if atom.GetAtomicNum() != 1:
                nbors_pos2.append(atom.GetIdx())
                
        for atom in nbors_pos2:
            if atom != pos_3:
                pos_1 = atom 
                break
                
        nbors_pos3 = []
        
        for atom in end_atom.GetAtoms():
            if atom.GetAtomicNum() != 1:
                nbors_pos3.append(atom.GetIdx())
       
        for atom in nbors_pos3:
            if atom != pos_2:
                pos_4 = atom
                break

        # returned as a list to be able to properly index the atoms
        # of the mda universe
        return [pos_1, pos_2, pos_3, pos_4]


    def get_torsions(self):
        # MO: should this be an iterator instead
        # for bond in self.oemol.GetBonds(IsRotor()):
        #     yield self._get_torsion(bond)

        # only get torsions for bonds that are rotatable 
        # rotatable bonds cannot be terminal
        torsions = []
        for bond in self.oemol.GetBonds(oechem.IsRotor()):
            torsion = self._get_torsion(bond)
            torsions.append(torsion)
        return torsions