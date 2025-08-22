from  slow_rotations import torsions as tor
from slow_rotations import mappings
import matplotlib.pyplot as plt
import tempfile
from slow_rotations import rdkit_wrapper as rdw
from PIL import Image
import numpy as np
import pandas as pd
import os
from scipy.special import kl_div

class TorsionComparator():
	def __init__(self, tf_list: list):
		'''
		Torsion Comparator

		Args:
			tf_list: [TorsionFinder]; list of TorsionFinder to compare

		Returns:
			None
		'''
		self.tf_list = tf_list
		self.torsions = self.tf_list[0].get_torsions()


	def get_torsions(self):
		''' 
		Get all small molecule torsions using the indexing system first input LigandTorsionFinder

		Args:
			None

		Returns:
			[[int, int, int, int]]: list of torsions 
		'''
		return self.torsions

	def convert_torsion_indices(self, trajidx, torsion):	
		'''
		Converts torsion from first indexing system into indexing system of another trajectory

		Args: 
			trajidx: int; index of the trajectory we are converting to
			torsion: [int, int, int, int]; torsion in the indexing system of the first trjaectory, (torsions from get_torsions function)


		Returns
			[int, int, int, int]: torsion in the indexing system of the trajectory specified
		'''
		return torsion

 
	def get_all_angles(self, torsion):
		'''
		gathers all the angles across all timeseries for the torsion of interest
		Args:
			torsion: [int, int, int, int]; torsion in the indexing system of the first trjaectory, (torsions from get_torsions function)

		Returns: 
			[float]: timeseries of angles
		'''
		angles = []

		for idx, tf in enumerate(self.tf_list):
			angles.extend(tf.get_torsion_angles(self.convert_torsion_indices(idx, torsion)).flatten())

		return angles


	@staticmethod
	def get_angle_min(angles):
		'''
		gets the minimum angle across all angles provided to shift histogram by such that no peaks are split
		Args:
			angles: [float]; list of angless

		Returns: 
			int: minimum angle
		'''
		return tor.TorsionFinder.get_angle_shift_point(angles)


	def get_cumulative_angles(self, torsion):
		'''
		gets all angles of the torsion across all trajectors, shifted such that no peaks are split
		Args:
			torsion: [int, int, int, int]; torsion in the indexing system of the first trjaectory, (torsions from get_torsions function)

		Returns: 
			int: minimum angle
			[float]: shifted angles across all trajectories
		'''
		cum_angles = self.get_all_angles(torsion)
		angle_min = TorsionComparator.get_angle_min(cum_angles)

		angle_min, shifted_cum_angles = tor.TorsionFinder.shift_torsion_angles(cum_angles, angle_min=angle_min)
		return angle_min, shifted_cum_angles

	def plot_cumulative_distribution(self, torsion, ax=None, save_path=None, close=False):
		'''
		creates a matplotlib plot of the distribution across all the trajectories passed
		Args:
			torsion: [int, int, int, int]; torsion in the indexing system of the first trjaectory, (torsions from get_torsions function)
			ax: axis; matplot lib axis to put the plot in
			save_path: str; path to save 
			close: bool; close the figure when done

		Returns: 
			int: angle minimum
			[float]: angle location of peaks
			[[float,float]]: bounds of peaks
		'''
		angle_min, cum_angles = self.get_cumulative_angles(torsion)

		X,scores,angle_min = tor.TorsionFinder.get_kde(cum_angles, angle_min=angle_min)

		num_peaks, peaks = tor.TorsionFinder.get_kde_num_peaks(X,scores)

		min_max = tor.TorsionFinder.get_bounds_mindist(X, scores, num_peaks, peaks)

		pdf_individual = []
		for mm in min_max:
			gmm,x,pdf,pdfi,bounds = tor.TorsionFinder.get_individual_gmm(X, angle_min, mm)
			pdf_individual.append(pdfi)

		if not ax:
			f, ax = plt.subplots()

		tor.TorsionFinder.plot_dihedral_histogram(cum_angles, ax=ax, show=False,  angle_min=angle_min, pdf_individual=pdf_individual)
		ax.set_title("Cumulative Angle Distribution Across Repeats")

		if save_path:
			plt.savefig(save_path)
			if close:
				plt.close()

		return angle_min, peaks, min_max



class LigandTorsionComparator(TorsionComparator):
	def __init__(self, tf_list: list):
		'''
		Ligand Comparator

		Args:
			tf_list: [TorsionFinder]; list of TorsionFinder to compare

		Returns:
			None
		'''
		self.tf_list = tf_list
		
		self.rdmols = []
		for tf in self.tf_list:
			self.rdmols.append(tf.get_rdmol())

		self.lig_mappings = {}
		self.sys_mappings = {}

		for idx,rdmol in enumerate(self.rdmols):
			self.lig_mappings[(0,idx)] = mappings.rd_map_mols(self.rdmols[0], rdmol)

		self.torsions = self.tf_list[0].get_torsions()

	def convert_torsion_indices(self, trajidx, torsion):	
		'''
		Converts torsion from first indexing system into indexing system of another trajectory

		Args: 
			trajidx: int; index of the trajectory we are converting to
			torsion: [int, int, int, int]; torsion in the indexing system of the first trjaectory, (torsions from get_torsions function)


		Returns
			[int, int, int, int]: torsion in the indexing system of the trajectory specified
		'''
		converted_torsion = []
		for t_idx in torsion:
			converted_torsion.append(self.tf_list[trajidx].convert_idx_to_sysidx(self.lig_mappings[(0,trajidx)][t_idx]))
		return converted_torsion


	def plot_all_distributions(self, torsion, save_path=None, close=False):
		'''
		creates a matplotlib plot of the distribution across all the trajectories passed as well as one for each individual trajectory for comparison to each other and the cumulative distribution
		Args:
			torsion: [int, int, int, int]; torsion in the indexing system of the first trjaectory, (torsions from get_torsions function)
			save_path: str; path to save 
			close: bool; close the figure when done

		Returns: 
			dict: torsion information for futher analysis
		'''
		num_cols = 3
		fig_width = 8*num_cols
		
		num_rows = len(self.tf_list) + 1
		fig_height = 6.25*num_rows
		pdf_colors = ['red', 'orange', 'green', 'blue', 'purple', 'brown']
		f,ax = plt.subplots(num_rows, num_cols, figsize=(fig_width, fig_height))

		ax[0,0].axis('off')

		with tempfile.NamedTemporaryFile(suffix='.png') as highlightpng:
			self.tf_list[0].highlight_dihedral(torsion, save_path=highlightpng)
			img = np.asarray(Image.open(highlightpng.name))
			ax[0,2].imshow(img)
			ax[0,2].axis('off')

		angle_min, peaks, min_max_cum = self.plot_cumulative_distribution(torsion, ax=ax[0,1], save_path=None, close=False)
		
		results = {}


		for idx,tf in enumerate(self.tf_list):
			results[idx] = {}
			torsion_sys = self.convert_torsion_indices(idx, torsion)

			symmetry = False
			if mappings.check_symmetry(self.rdmols[idx], torsion):
				# if there is symmetry add a note on the image
				symmetry = True

			# Analyze torsion Data
			angles = tf.get_torsion_angles(torsion_sys).flatten()
			X,scores,angle_min = tor.TorsionFinder.get_kde(angles, angle_min=angle_min)
			min_max= tor.TorsionFinder.get_bounds_mindist(X, scores, len(peaks), peaks)
			for i,mm in enumerate(min_max):
				if mm == (None, None):
					min_max[i] = min_max_cum[i]
			angle_min, shifted_angles = tor.TorsionFinder.shift_torsion_angles(angles, angle_min=angle_min)
			
			pdf_individual = []
			missing_states_rpt = {}
			for si, mm in enumerate(min_max):
				try:
					gmm,x,pdf,pdfi,bounds = tor.TorsionFinder.get_individual_gmm(X, angle_min, mm)
					pdf_individual.append(pdfi)
					missing_states_rpt[si] = False
				
				except ValueError:
					pdf_individual.append(np.array([[],[]]))
					# missing a state
					missing_states_rpt[si] = True


			# plot histogram
			tor.TorsionFinder.plot_dihedral_histogram(angles, ax=ax[idx+1,1], show=False, angle_min=angle_min, pdf_individual=pdf_individual, pdf_colors=pdf_colors)

			# plot transitions
			transition_ctr = tor.TorsionFinder.transition_matrix(shifted_angles, min_max)
			tor.TorsionFinder.plot_transition_counts(transition_ctr, ax=ax[idx+1,2], colors=pdf_colors)
			transition_populations = tor.TorsionFinder.state_populations(angles, min_max)

			# plot scatter
			tor.TorsionFinder.plot_dihedral_scatter(shifted_angles, ax =ax[idx+1,0], angle_min=angle_min)

			label = f'Repeat {idx+1}'
			ax[idx+1,0].text(-0.2, 0.5, label, va='center', ha='center', rotation='vertical', fontsize=20, transform=ax[idx+1,0].transAxes)

			results[idx]['residue'] = self.tf_list[idx].get_residue_name_from_torsion(torsion_sys)
			results[idx]['torsion_idx'] = torsion
			results[idx]['torsion_system_idx'] = torsion_sys
			results[idx]['symmetry'] = symmetry
			results[idx]['missing_states'] = missing_states_rpt
			results[idx]['transitions'] = transition_ctr.to_dict()


		if symmetry:
			ax[0,1].text(0, -0.2, f'** Warning: torsion has symmetry, disregard transitions', fontsize=15, color='red', transform=ax[0,2].transAxes)

		if save_path:
			plt.savefig(save_path)

		return results



class ProteinTorsionComparator(TorsionComparator):
	def __init__(self, tf_list: list, a_cutoff: float):
		'''
		Protein Comparator

		Args:
			tf_list: [TorsionFinder]; list of TorsionFinder to compare
			a_cutoff: angstrom cutoff to consider the binding site

		Returns:
			None
		'''
		self.tf_list = tf_list

		self.a_cutoff = a_cutoff

		self.torsions = self.get_torsions()


	def get_torsions(self):
		'''
		gets all protein sidechain torsions in the system

		Args:
			None

		Returns: 
			[[int, int, int, int]]: list of torsions
		'''

		# torsions = None
		# for i in range(8): 
		# 	print(torsions)
		# 	chix_torsions = self.tf_list[0].get_chi_x_torsions(i, a_cutoff=self.a_cutoff)
		# 	print(chix_torsions)
		# 	if chix_torsions != []:
		# 		if not torsions:
		# 			torsions = chix_torsions

		# 		else:
		# 			np.concatenate(torsions, self.tf_list[0].get_chi_x_torsions(i, a_cutoff=self.a_cutoff))
		
		#i = self.tf_list[0].get_chi_x_torsions(3, a_cutoff=self.a_cutoff)
		#print(i)
		return self.tf_list[0].get_all_chi_x_torsions(a_cutoff=self.a_cutoff)




	def plot_all_distributions(self, torsion, save_path=None, close=False):
		'''
		creates a matplotlib plot of the distribution across all the trajectories passed as well as one for each individual trajectory for comparison to each other and the cumulative distribution
		Args:
			torsion: [int, int, int, int]; torsion in the indexing system of the first trjaectory, (torsions from get_torsions function)
			save_path: str; path to save 
			close: bool; close the figure when done

		Returns: 
			dict: torsion information for futher analysis
		'''
		num_cols = 3
		fig_width = 8*num_cols
		
		num_rows = len(self.tf_list) + 1
		fig_height = 6.25*num_rows
		pdf_colors = ['red', 'orange', 'green', 'blue', 'purple', 'brown']
		f,ax = plt.subplots(num_rows, num_cols, figsize=(fig_width, fig_height))

		ax[0,0].axis('off')

		with tempfile.NamedTemporaryFile(suffix='.png') as highlightpng:
			self.tf_list[0].highlight_dihedral(torsion, save_path=highlightpng)
			img = np.asarray(Image.open(highlightpng.name))
			ax[0,2].imshow(img)
			ax[0,2].axis('off')

		angle_min, peaks, min_max_cum = self.plot_cumulative_distribution(torsion, ax=ax[0,1], save_path=None, close=False)
		
		results = {}


		for idx,tf in enumerate(self.tf_list):
			results[idx] = {}
			torsion_sys = self.convert_torsion_indices(idx, torsion)

			# Analyze torsion Data
			angles = tf.get_torsion_angles(torsion_sys).flatten()
			X,scores,angle_min = tor.TorsionFinder.get_kde(angles, angle_min=angle_min)
			min_max= tor.TorsionFinder.get_bounds_mindist(X, scores, len(peaks), peaks)
			for i,mm in enumerate(min_max):
				if mm == (None, None):
					min_max[i] = min_max_cum[i]
			angle_min, shifted_angles = tor.TorsionFinder.shift_torsion_angles(angles, angle_min=angle_min)
			
			pdf_individual = []
			missing_states_rpt = {}
			for si, mm in enumerate(min_max):
				try:
					gmm,x,pdf,pdfi,bounds = tor.TorsionFinder.get_individual_gmm(X, angle_min, mm)
					pdf_individual.append(pdfi)
					missing_states_rpt[si] = False
				
				except ValueError:
					pdf_individual.append(np.array([[],[]]))
					# missing a state
					missing_states_rpt[si] = True


			# plot histogram
			tor.TorsionFinder.plot_dihedral_histogram(angles, ax=ax[idx+1,1], show=False, angle_min=angle_min, pdf_individual=pdf_individual, pdf_colors=pdf_colors)

			# plot transitions
			transition_ctr = tor.TorsionFinder.transition_matrix(shifted_angles, min_max)
			tor.TorsionFinder.plot_transition_counts(transition_ctr, ax=ax[idx+1,2], colors=pdf_colors)
			transition_populations = tor.TorsionFinder.state_populations(angles, min_max)

			# plot scatter
			tor.TorsionFinder.plot_dihedral_scatter(shifted_angles, ax =ax[idx+1,0], angle_min=angle_min)

			label = f'Repeat {idx+1}'
			ax[idx+1,0].text(-0.2, 0.5, label, va='center', ha='center', rotation='vertical', fontsize=20, transform=ax[idx+1,0].transAxes)

			results[idx]['residue'] = self.tf_list[idx].get_residue_name_from_torsion(torsion)
			results[idx]['chi'] = self.tf_list[idx].determine_chi_angle(torsion)
			results[idx]['torsion_idx'] = [int(t) for t in torsion]
			results[idx]['torsion_system_idx'] = [-1,-1,-1,-1]
			results[idx]['missing_states'] = missing_states_rpt
			results[idx]['symmetry'] = False
			results[idx]['transitions'] = transition_ctr.to_dict()

		if save_path:
			plt.savefig(save_path)

		return results


































 










		






