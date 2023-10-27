import torsions as tor
import mappings
import matplotlib.pyplot as plt
import tempfile
import rdkit_wrapper as rdw
from PIL import Image
import numpy as np
import pandas as pd
import os

class LigandComparator():
    def __init__(self, gas_ltf, bnd_ltf):
        ''' gas_ltf: torsions.LigandTorsionFinder for gas phase ligand "truth"
            bnd_ltf: torsions.LigandTorsionFinder for ligand bound to complex
        '''
        self.gas_ltf = gas_ltf
        self.bnd_ltf = bnd_ltf

        self.gas_rdmol = self.gas_ltf.get_rdmol()
        self.bnd_rdmol = self.bnd_ltf.get_rdmol()

        # gas_bnd_map index mapping 
        self.gas_bnd_map = mappings.rd_map_mols(self.gas_rdmol, self.bnd_rdmol)

        # gas_tor and bnd_tor: nested list of lists shape (n,4)
        # each nested list is 4 ints long, corresponding to one torsion
        self.gas_tor = self.gas_ltf.get_torsions()
        self.bnd_tor = [ mappings.convert_dihedral(self.gas_bnd_map, i) for i in self.gas_tor ]

    def get_gas_torsions(self):
        return self.gas_tor

    def get_bound_torsions(self):
        return self.bnd_tor

    def compare_torsions(self, gas_tor):
        '''
            Returns
            =======
            gas_angles:    shifted angles for the gas phase simulation
            gas_num_peaks: number of states in the gas phase
            gas_min_max:   boundaries of the gas phase
            bnd_angles:    shifted angles for the ligand bound complex 
                           simulation
        '''
        bnd_tor = mappings.convert_dihedral(self.gas_bnd_map, gas_tor)

        gas_X, gas_scores, gas_angle_min = self.gas_ltf.get_kde(gas_tor)

        print(gas_X)

        gas_num_peaks = self.gas_ltf.get_kde_num_peaks(gas_scores)

        gas_min_max = self.gas_ltf.get_bounds_knn(gas_X, gas_num_peaks)

        gas_angles = self.gas_ltf.shift_torsion_angles(gas_tor, angle_min=gas_angle_min)[1].flatten()
        bnd_angles = self.bnd_ltf.shift_torsion_angles(bnd_tor,angle_min=gas_angle_min)[1].flatten()
        return gas_X, gas_angles, gas_num_peaks, gas_min_max, gas_angle_min, bnd_angles

    def plot_compared_dihedrals_histogram(self, gas_angles, bnd_angles, angle_min, ax=None, num_bins=60, alpha=0.5, show=True, gascolor=None, bndcolor=None):
        self.gas_ltf.plot_dihedral_histogram(base_X, base_angle_min, ax=ax, color=gascolor)
        self.gas_ltf.plot_dihedral_histogram(base_X, base_angle_min, ax=ax, color=bndcolor)
        plt.title("Gas Phase Ligand vs. Bound Ligand")
        plt.legend(["Gas Phase", "Bound"])


    def make_compare_torsions_img(self, gas_tor, show=False):
        bnd_tor = mappings.convert_dihedral(self.gas_bnd_map, gas_tor)
        gas_X, gas_angles, gas_num_peaks, gas_min_max, gas_angle_min, bnd_angles = self.compare_torsions(gas_tor)

        gas_pdf_individual = []
        for min_max in gas_min_max:
        	gmm,x,pdf,pdf_individual = self.gas_ltf.get_individual_gmm(gas_X, gas_angle_min, min_max)
        	gas_pdf_individual.append(pdf_individual)

        g1,g2,g3,g4 = tuple(gas_tor)
        b1,b2,b3,b4 = tuple(mappings.convert_dihedral(self.gas_bnd_map, gas_tor))

        sup_title = f"Gas Phase Ligand ({g1},{g2},{g3},{g4}) vs. Bound Ligand ({b1},{b2},{b3},{b4})"
        gas_title = f"Gas Phase Ligand"
        bnd_title = f"Bound Ligand"

        f,ax = plt.subplots(2, 3, figsize=(25, 12.5))
        f.suptitle(sup_title,fontsize=60)
        f.tight_layout(pad=3.5)

        states_list = [f"s{i}" for i in range(gas_num_peaks)]

        highlightpng = tempfile.NamedTemporaryFile(suffix='.png', delete=False)


        self.gas_ltf.plot_dihedral_histogram(gas_tor, angles=gas_angles, angle_min=gas_angle_min, ax=ax[0,1], pdf_individual=gas_pdf_individual, show=False, title=gas_title,color='tab:blue')
        ax[0,1].legend(states_list)
        self.gas_ltf.plot_dihedral_histogram(bnd_tor, angles=bnd_angles, angle_min=gas_angle_min, ax=ax[1,1], pdf_individual=gas_pdf_individual, show=False, title=bnd_title,color='purple')
        ax[1,1].legend(states_list)
        self.gas_ltf.plot_dihedral_histogram(gas_tor, angles=gas_angles, angle_min=gas_angle_min, ax=ax[1,0], show=False, title=gas_title, color='tab:blue')
        self.gas_ltf.plot_dihedral_histogram(bnd_tor, angles=bnd_angles, angle_min=gas_angle_min, pdf_individual=gas_pdf_individual, ax=ax[1,0], show=False, title=bnd_title, color='purple')
        ax[1,0].set_title("Gas Phase Ligand vs. Bound Ligand")
        ax[1,0].legend(states_list + ["Gas", "Bound"])


        rdw.highlight_dihedral(self.gas_rdmol, gas_tor, highlightpng.name)
        img = np.asarray(Image.open(highlightpng.name))
        ax[0,0].imshow(img)
        ax[0,0].axis('off')

        gas_transition = self.gas_ltf.transition_counter(gas_angles, gas_min_max)
        gas_pd = pd.DataFrame(data=gas_transition, columns=['transitions in', 'transitions out'])
        table = ax[0,2].table(
            cellText=gas_pd.values, 
            colLabels=gas_pd.columns, 
            rowLabels=[f"s{i}" for i in range(len(gas_pd))],
            loc='center',
        )
        cells = table.properties()["celld"]
        for i in range(0, len(gas_pd)+1):
            cells[i, 0].set_text_props(ha="center")
            cells[i, 1].set_text_props(ha="center")


        bnd_transition = self.gas_ltf.transition_counter(bnd_angles, gas_min_max)
        bnd_pd = pd.DataFrame(data=bnd_transition, columns=['transitions in', 'transitions out'])
        table = ax[1,2].table(
            cellText=bnd_pd.values, 
            colLabels=bnd_pd.columns, 
            rowLabels=[f"s{i}" for i in range(len(bnd_pd))],
            loc='center',
        )
        cells = table.properties()["celld"]
        for i in range(0, len(bnd_pd)+1):
            cells[i, 0].set_text_props(ha="center")
            cells[i, 1].set_text_props(ha="center")

        ax[0,2].axis('off')
        ax[1,2].axis('off')

        if show:
            plt.show()
        os.unlink(highlightpng.name)












