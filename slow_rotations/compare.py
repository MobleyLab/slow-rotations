import torsions as tor
import mappings
import matplotlib.pyplot as plt
import tempfile
import rdkit_wrapper as rdw
from PIL import Image
import numpy as np
import pandas as pd
import os
from scipy.special import kl_div


class LigandComparisonResult():
    def __init__(
            self, 
            gas_tor, 
            gas_X, 
            gas_angles, 
            gas_num_peaks, 
            gas_min_max, 
            gas_angle_min, 
            gas_pdf_individual, 
            gas_transition, 
            bnd_tor, 
            bnd_angles, 
            bnd_num_peaks,
            bnd_min_max, 
            bnd_angle_min,
            bnd_pdf_individual,
            bnd_transition,
        ):
        self.gas_tor=gas_tor
        self.gas_X = gas_X
        self.gas_angles = gas_angles
        self.gas_num_peaks = gas_num_peaks
        self.gas_min_max = gas_min_max
        self.gas_angle_min = gas_angle_min
        self.gas_pdf_individual = gas_pdf_individual
        self.gas_transition = gas_transition
        self.bnd_tor = bnd_tor
        self.bnd_angles = bnd_angles
        self.bnd_angle_min = bnd_angle_min
        self.bnd_num_peaks = bnd_num_peaks
        self.bnd_min_max = bnd_min_max
        self.bnd_pdf_individual = bnd_pdf_individual
        self.bnd_transition = bnd_transition

    def to_dict(self):
        result = {
            "gas": {
                "indices": self.gas_tor,
                "num_peaks": self.gas_num_peaks,
                "transitions": self.gas_transition.to_dict(),
            },
            "bound": {
                "indices": self.bnd_tor,
                "num_peaks": self.bnd_num_peaks,
                "transitions": self.bnd_transition.to_dict()
            },
        }
        return result


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
        bnd_tor = [self.bnd_ltf.convert_ligidx_to_sysidx(i) for i in mappings.convert_dihedral(self.gas_bnd_map, gas_tor)]

        gas_X, gas_scores, gas_angle_min = self.gas_ltf.get_kde(gas_tor)

        gas_num_peaks, gas_peaks = self.gas_ltf.get_kde_num_peaks(gas_scores)

        gas_min_max = self.gas_ltf.get_bounds_mindist(gas_X, gas_num_peaks, gas_peaks)

        gas_angles = self.gas_ltf.shift_torsion_angles(gas_tor, angle_min=gas_angle_min)[1].flatten()
        

        gas_pdf_individual = []
        for min_max in gas_min_max:
            gmm,x,pdf,pdf_individual,bounds = self.gas_ltf.get_individual_gmm(gas_X, gas_angle_min, min_max)
            gas_pdf_individual.append(pdf_individual)

        gas_transition = self.gas_ltf.transition_counter(gas_angles, gas_min_max)



        bnd_X, bnd_scores, bnd_angle_min = self.bnd_ltf.get_kde(bnd_tor)

        bnd_angles = self.bnd_ltf.shift_torsion_angles(bnd_tor,angle_min=bnd_angle_min)[1].flatten()
        bnd_num_peaks, bnd_peaks = self.bnd_ltf.get_kde_num_peaks(bnd_scores)
        bnd_min_max = self.bnd_ltf.get_bounds_mindist(bnd_X, bnd_num_peaks, bnd_peaks)
        bnd_transition = self.bnd_ltf.transition_counter(bnd_angles, bnd_min_max)

        bnd_pdf_individual = []
        for min_max in bnd_min_max:
            gmm,x,pdf,pdf_individual,bounds = self.bnd_ltf.get_individual_gmm(bnd_X, bnd_angle_min, min_max)
            bnd_pdf_individual.append(pdf_individual)



        return LigandComparisonResult(
            gas_tor, gas_X, gas_angles, gas_num_peaks, gas_min_max, gas_angle_min, gas_pdf_individual, gas_transition, 
            bnd_tor, bnd_angles, bnd_num_peaks, bnd_min_max, bnd_angle_min, bnd_pdf_individual, bnd_transition
        )

    def get_individual_hist(self, angles, min_max):
        min_angle = min_max[0]
        max_angle = min_max[1]
        angles_np = np.array(angles)
        angles_np_1pk = angles_np[np.logical_and(angles_np >= min_angle, angles_np <= max_angle)]
        return angles_np_1pk


    def are_comparable(self, lig_comp_result: LigandComparisonResult, no_peak_threshold=0.05):

        kl_div_sums = list()

        gas_num_angles = len(lig_comp_result.gas_angles)
        bnd_num_angles = len(lig_comp_result.bnd_angles)

        for min_max in lig_comp_result.gas_min_max:
            gas_angles_1pk = self.get_individual_hist(lig_comp_result.gas_angles, min_max)

            bnd_angles_1pk = self.get_individual_hist(lig_comp_result.bnd_angles, min_max)
            if len(bnd_angles_1pk)/bnd_num_angles < no_peak_threshold:
                # condition that determines if there are enough samples in the region 
                # to even compare or if it should be labeled as the peak is not visited
                continue

            kl_div_sums.append(sum(kl_div(gas_angles_1pk, bnd_angles_1pk)))

        return kl_div_sums



    def plot_compared_dihedrals_histogram(self, gas_angles, bnd_angles, angle_min, ax=None, num_bins=40, alpha=0.5, show=True, gascolor=None, bndcolor=None):
        self.gas_ltf.plot_dihedral_histogram(base_X, base_angle_min, ax=ax, color=gascolor)
        self.gas_ltf.plot_dihedral_histogram(base_X, base_angle_min, ax=ax, color=bndcolor)
        plt.title("Solvent Ligand vs. Bound Ligand")
        plt.legend(["Solvent", "Bound"])



    def make_compare_torsions_img(self, lig_comparison_result, show=False, save_path=None):

        gas_tor=lig_comparison_result.gas_tor 
        gas_X=lig_comparison_result.gas_X
        gas_angles=lig_comparison_result.gas_angles
        gas_num_peaks=lig_comparison_result.gas_num_peaks
        gas_min_max=lig_comparison_result.gas_min_max
        gas_angle_min=lig_comparison_result.gas_angle_min
        bnd_tor=lig_comparison_result.bnd_tor
        bnd_angles=lig_comparison_result.bnd_angles
        gas_pdf_individual=lig_comparison_result.gas_pdf_individual
        gas_transition=lig_comparison_result.gas_transition
        bnd_angle_min=lig_comparison_result.bnd_angle_min
        bnd_min_max = lig_comparison_result.bnd_min_max
        bnd_pdf_individual = lig_comparison_result.bnd_pdf_individual
        bnd_transition=lig_comparison_result.bnd_transition

        g1,g2,g3,g4 = tuple(gas_tor)
        b1,b2,b3,b4 = tuple(bnd_tor)

        sup_title = f"Solvent Ligand vs. Bound Ligand"
        gas_title = f"Solventv Ligand ({g1},{g2},{g3},{g4}) "
        bnd_title = f"Bound Ligand ({b1},{b2},{b3},{b4})"

        f,ax = plt.subplots(2, 3, figsize=(18, 9))
        f.suptitle(sup_title,fontsize=60)
        f.tight_layout(pad=3.5)

        states_list = [f"s{i}" for i in range(gas_num_peaks)]

        highlightpng = tempfile.NamedTemporaryFile(suffix='.png', delete=False)

        self.gas_ltf.plot_dihedral_histogram(gas_tor, angles=gas_angles, angle_min=gas_angle_min, ax=ax[0,1], pdf_individual=gas_pdf_individual, show=False, title=gas_title,color='tab:blue')
        ax[0,1].legend(states_list)
        self.gas_ltf.plot_dihedral_histogram(bnd_tor, angles=bnd_angles, angle_min=bnd_angle_min, ax=ax[1,1], pdf_individual=bnd_pdf_individual, show=False, title=bnd_title,color='purple')
        ax[1,1].legend(states_list)
        self.gas_ltf.plot_dihedral_histogram(gas_tor, angles=gas_angles, angle_min=gas_angle_min, ax=ax[1,0], show=False, title=gas_title, color='tab:blue')
        pdf_colors = ['red', 'orange', 'green', 'blue', 'purple', 'brown']
        bnd_angles_as_gas = self.bnd_ltf.shift_torsion_angles(bnd_tor,angle_min=gas_angle_min)[1].flatten()
        self.gas_ltf.plot_dihedral_histogram(bnd_tor, angles=bnd_angles_as_gas, angle_min=gas_angle_min, pdf_individual=gas_pdf_individual, ax=ax[1,0], show=False, title=bnd_title, color='purple', pdf_colors=pdf_colors)
        ax[1,0].set_title("Solvent Ligand vs. Bound Ligand")
        ax[1,0].legend(states_list + ["Solvent", "Bound"])


        rdw.highlight_dihedral(self.gas_rdmol, gas_tor, highlightpng.name)
        img = np.asarray(Image.open(highlightpng.name))
        ax[0,0].imshow(img)
        ax[0,0].axis('off')

        self.gas_ltf.plot_transition_counts(gas_transition, ax=ax[0,2],colors=pdf_colors)

        self.gas_ltf.plot_transition_counts(bnd_transition, ax=ax[1,2], colors=pdf_colors)


    
        # ax[0,2].axis('off')
        # ax[1,2].axis('off')

        if show:
            plt.show()

        if save_path:
            plt.savefig(save_path)

        os.unlink(highlightpng.name)

        self.bnd_ltf.plot_kde(bnd_tor, save_path = str(save_path)+"kde")





