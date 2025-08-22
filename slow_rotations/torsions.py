from pathlib import Path
import tempfile
import math
import warnings
import os

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD, GRO, XTC
from MDAnalysis.analysis.dihedrals import Dihedral
from MDAnalysis.analysis import distances
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from PIL import Image
from sklearn.mixture import GaussianMixture
from scipy.signal import find_peaks
from sklearn.neighbors import KernelDensity
from sklearn.cluster import KMeans
from rdkit import Chem
from pymbar import timeseries

from rdkit.Chem.rdchem import AtomValenceException

from slow_rotations import utils
from slow_rotations import mappings
from slow_rotations import rdkit_wrapper as rdw
from slow_rotations import pdb_writer as pdbw
from slow_rotations import molconverter as mc
from slow_rotations import transitions


class BadTorsionError(Exception):
	'''
	Exception raised when there is a problem with the user input torsion
	'''
	pass

class TorsionFinder():
	'''
	Base class for performing analysis on all relevant torsions in a system
	'''

	def __init__(self, trajf:str, topf: str):
		'''
		Args:
			trajf: str; trajectory file
			topf: str; topology file

		Returns:
			float: KL divergence for each peak
		'''
		self.trajf = trajf
		self.topf = topf
		self.mda_universe = mda.Universe(self.topf, self.trajf)

		# solute = self.mda_universe.select_atoms("not resname HOH")

		# with tempfile.NamedTemporaryFile(suffix='.xtc') as no_solvent_xtc:
		# 	with tempfile.NamedTemporaryFile(suffix='.gro') as no_solvent_gro:

		# 		print("writing out set_xticks")
		# 		with mda.Writer(no_solvent_xtc.name, solute.n_atoms) as W:
		# 			for ts in self.mda_universe.trajectory:
		# 				W.write(solute)

		# 		print("writing out gro")
				
		# 		with mda.Writer(no_solvent_gro.name, solute.n_atoms) as W:
		# 			for ts in self.mda_universe.trajectory[0]:
		# 				W.write(solute)

		# 		print("finished writing out gro")

		# 		self.mda_universe = mda.Universe(no_solvent_gro.name, no_solvent_xtc.name)

	def get_torsions():
		raise utils.NotImplementedError

	@staticmethod
	def get_angle_shift_point(angles, num_bins=40):
		'''
		finds a new angle minimum for histogram so no peaks are split
		If no bins exist without histogram entries, shifts to the bin with the lowest number of entries.

		Args:
			angles: [float]; angles timeseries
			num_bins: int; number of bins to use for histogram of angles

		Returns:
			int: minimum angle for histogram 		
		'''
		angles_of_bins = np.linspace(-180,180,num_bins+1)
		hist = np.histogram(angles, angles_of_bins)[0]
		
		max_streak = 0
		current_streak = 0
		first_zero_index = -1
		first_zero_index_of_max_streak = 0
		last_zero_index_of_max_streak = -1
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
					last_zero_index_of_max_streak = i
				current_streak = 0
		if current_streak > max_streak:
			first_zero_index_of_max_streak = first_zero_index
			last_zero_index_of_max_streak = i

		angle_min = angles_of_bins[int((first_zero_index_of_max_streak + last_zero_index_of_max_streak)/2)]
		
		return angle_min

	def get_torsion_angles(self, torsion:[int,int,int,int]):
		'''
		Gets timeseries of torsion angles

		Args:
			torsion: [int, int, int, int]; list of torsion atom indices

		Returns:
			np.array: angle timeseries	
		'''

		ags = [mda.AtomGroup(self.mda_universe.atoms[(torsion)])]
		Run = Dihedral(ags).run()
		shape = (Run.results.angles.shape)
		angles = Run.results.angles
		return angles


	def convert_idx_to_sysidx(self, idx):
		'''
		converts the index to 0 based index within molecule to index of the entire system

		Args:
			idx: int; index to convert

		Returns:
			int: converted index
		'''		
		
		return idx


	@staticmethod
	def shift_torsion_angles(angles, num_bins=40, angle_min=None):
		'''
		shifts torsion angles in histogram to specified angle. if no angle is specified, shifts angle to 

		Args:
			angles: [float]; timeseries of angles
			num_bins: int; num bins in the historgram
			angle_min: int; force a minimum for the angles if specified otherwise 

		Returns:
			float: minimum angle of histogram
		'''	
		if angle_min == None:
			angle_min = TorsionFinder.get_angle_shift_point(angles, num_bins)

		angle_shift = abs(-180 - angle_min)

		angle_shift_rnd20 = math.ceil(angle_shift / 20) * 20

		angle_min = -180 + angle_shift_rnd20

		shifted_angles = np.array(angles)
		for i in range(len(shifted_angles)):
			if shifted_angles[i] < (angle_min + 360//num_bins):
				shifted_angles[i] = shifted_angles[i] + 360
		return angle_min, shifted_angles


	# def shift_torsion_angles(self, torsion:[int,int,int,int], num_bins=40, angle_min=None):
	# 	''' 
	# 	Args:
	# 		torsion: [int, int, int, int]; atom indices of atoms in torsion
	# 		num_bins: int; number of bins in the histogram so you know where the shift point is
	# 		angle_min: float; is the minimum angle of the histogram

	# 	Returns:
	# 		int: new angle minimum for histogram so no peaks are split
			
	# 	'''
	# 	angles = np.array(self.get_torsion_angles(torsion))
	# 	if angle_min == None:
	# 		angle_min = TorsionFinder.get_angle_shift_point(angles, num_bins)

	# 	shifted_angles = np.array(angles)
	# 	for i in range(len(shifted_angles)):
	# 		if shifted_angles[i] < (angle_min + 360//num_bins):
	# 			shifted_angles[i] = shifted_angles[i] + 360
	# 	return angle_min, shifted_angles


	@staticmethod
	def add_0_to_kde(X, score):
		'''
		shifts torsion angles in histogram to specified angle. if no angle is specified, shifts angle to 

		Args:
			angles: [float]; timeseries of angles
			num_bins: int; num bins in the historgram
			angle_min: int; force a minimum for the angles if specified otherwise 

		Returns:
			int: converted index
		'''	

		score = np.array(score) - np.min(np.array(score))
		new_X = [X[0]]
		new_score = [score[0]]

		prev = X[0]
		for i in range(1, len(X)):
			if X[i] - X[i-1] > 10:
				for newx in range(int(X[i-1]), int(X[i])):
					new_X.append(newx)
					new_score.append(0)

			
			new_X.append(X[i])
			new_score.append(score[i])

		return new_X, new_score


	@staticmethod
	def get_statistical_inefficiency(angles, num_bins=40, angle_min=None):
		'''
		calculates the statistical inefficiency of the angles timeseries 

		Args:
			angles: [float]; timeseries of angles
			num_bins: int; num bins in the historgram
			angle_min: int; force a minimum for the angles if specified otherwise 

		Returns:
			float: statisitcal inefficiency
		'''	

		angle_min, X = TorsionFinder.shift_torsion_angles(angles, num_bins=num_bins, angle_min=angle_min)

		g = timeseries.statistical_inefficiency(X.flatten())

		return g


	@staticmethod
	def get_kde(angles, num_bins=40, angle_min=None):
		'''
		calculates the kernel density estimator

		Args:
			angles: [float]; timeseries of angles
			num_bins: int; num bins in the historgram
			angle_min: int; force a minimum for the angles if specified otherwise 

		Returns:
			[[float]]: X angles
			[float]: kernel density estimator score
			int: minimum angle of the distribution
		'''	
		angle_min, X = TorsionFinder.shift_torsion_angles(angles, num_bins=num_bins, angle_min=angle_min)

		X = np.sort(X.flatten())
		X = X.reshape(-1, 1)
		kde = KernelDensity(kernel='gaussian', bandwidth=360/num_bins).fit(X)

		score = kde.score_samples(X)


		return X,score,angle_min

	@staticmethod
	def get_kde_num_peaks(X, scores, smoothing_window=30, peak_prominence=0.0001, peak_dist=30):
		'''
		counts the number of peaks in the Kernel Density estimator
		Args:
			X: [[float]] angles
			scores: [float] kernel density estimator score
			smoothing_window: int width of window to average over
			peak_prominence: float how prominent the peak must be in order to be counted
			peak_dist: int minimum distance apart peaks must fall
		Returns:
			int: num_peaks
			peaks: [float] angle location of the peaks 
		'''	
		# was 50 and 0.008
		# moving average kde scores
		peak_prominence=0.05 

		new_x, new_score = X, scores #self.add_0_to_kde(X, scores)

		moving_avg = []
		for i in range(len(new_score) - smoothing_window):
			moving_avg.append(np.mean(new_score[i:i+smoothing_window]))

		# peaks from moving average
		peak_info = find_peaks(new_score, prominence=peak_prominence, distance=peak_dist)
		peaks = peak_info[0]

		peaks = [new_x[p]+0.25*smoothing_window for p in peaks]

		num_peaks = len(peaks)

		return num_peaks, peaks


	def mindist(self, X, num_components, peaks):
		pass


	@staticmethod
	def get_closest_peak(angle, peaks):
		'''
		calculate the closest peak for that angle
		Args:
			angle: float; angle to get closest peak for
			peaks: [float]; angle locations of peaks
		Returns:
			int: index of peak in peaks that is the minimum
		'''	
		peaks = np.array(peaks)

		peaks = np.absolute(peaks - angle)

		return np.argmin(peaks)



	@staticmethod
	def get_bounds_mindist(X, scores, num_components, peaks, tolerance=10):
		'''
		get the bounds based of each peak by assigning each score value to its closest peak
		Args:
			X: [[float]]; angles
			scores: [float]; kernel density estimator score
			num_components: int; nubmer of peaks
			peaks: [float]; angle locations of peaks
			tolerance: int; the largest gap in angles we can have before it is considered not in that peak
		Returns:
			[[float, float]]: angle bounds of each peak
		'''	

		angles = X[:,0]

		angles, scores = TorsionFinder.add_0_to_kde(X, scores)

		cluster_labels = []

		n_clusters = len(peaks)

		new_angles = []
	
		# if the kde score is 0, we don't want stop the inclusion
		for i,a in enumerate(angles):
			if scores[i] == 0:
				continue

			else:
				cluster_labels.append(TorsionFinder.get_closest_peak(a, peaks))
				new_angles.append(a)

		min_max = []

		for c in range(n_clusters):
			found_c = False
			for i in cluster_labels:
				if i == c:
					found_c = True
			center = peaks[c]

			centroid_index = None
			centroid_min_diff = 999999999
			for i,a in enumerate(new_angles):
				diff = abs(center - a)
				if diff < centroid_min_diff:
					centroid_min_diff = diff
					centroid_index = i
			min_i = 0
			max_i = len(new_angles) - 1
			idx = centroid_index

			while idx > 0:
				if abs(new_angles[idx] - new_angles[idx-1]) > tolerance or cluster_labels[idx] != c:
					min_i = idx
					break
				idx -= 1

			for idx in range(centroid_index,len(new_angles) - 1):
				if abs(new_angles[idx+1] - new_angles[idx]) > tolerance or cluster_labels[idx] != c:
					max_i = idx - 1
					break
			if found_c:
				min_max.append((new_angles[min_i],new_angles[max_i]))
			else:
				min_max.append((None, None))
		return min_max


	@staticmethod
	def get_bounds_knn(X, num_components, peaks, tolerance=30):
		'''
		get the bounds based of each peak using k nearest neighbor algorithm
		Args:
			X: [[float]]; angles
			num_components: int; nubmer of peaks
			peaks: [float]; angle locations of peaks
			tolerance: int; the largest gap in angles we can have before it is considered not in that peak
		Returns:
			[[float, float]]: angle bounds of each peak
		'''	

		init_centers = [[p,0] for p in peaks]

		data = np.column_stack((X[:,0], np.zeros(X.shape)))
		kmeans = KMeans(n_clusters=num_components, init=init_centers, max_iter=1)
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


	@staticmethod
	def get_individual_gmm(X, angle_min, min_max):
		'''
		get gaussian mixture model representation of a single peak 
		Args:
			X: [[float]]; angles
			angle_min: int; minimum angle 
			min_max: [[float, float]] bounds of the peaks
		Returns:
			gmm: gaussian mixture model object
			x:
			pdf: probability distribution function
			pdf_individual: probability distribution function of the individual peak
			[float, float]: minimum and maximum used
		'''	
		num_components = 1
		min_bnd = min_max[0]
		max_bnd = min_max[1]
		min_max_range = max_bnd - min_bnd

		single_state_X = list()

		min_idx = None
		max_idx = None

		for i,point in enumerate(X):
			if point[0] >= min_bnd and point[0] < max_bnd:
				single_state_X.append(point)

			if point[0] >= min_bnd and not min_idx:
				min_idx = i

			if point[0] >= max_bnd and not max_idx:
				max_idx = i - 1

		if not max_idx:
			max_idx = len(X) - 1

		mm = [min_idx, max_idx]

		gmm, x, pdf, pdf_individual = TorsionFinder.get_gmm(single_state_X, num_components, angle_min)

		frac = len(single_state_X)/len(X)
		pdf_individual = pdf_individual*frac
		pdf_individual = np.array([x.flatten(), pdf_individual.flatten()])
		
		return gmm,x,pdf,pdf_individual,mm


	@staticmethod
	def get_gmm(X, num_components, angle_min):
		'''
		get gaussian mixture model representation the entire distribution
		Args:
			X: [[float]]; angles
			angle_min: int; minimum angle 
			min_max: [[float, float]] bounds of the peaks
		Returns:
			gmm: gaussian mixture model object
			x:
			pdf: probability distribution function
			pdf_individual: probability distribution function of the individual peakd
		'''	
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

	@staticmethod
	def transition_counter(angles, range_of_states):
		'''
		count the transitions in OR out of each state (less information than transition_matrix())
		Args:
			[float]: angles
			range_of_states: [[float, float]]; bounds of the peaks 

		Returns:
			TransitionCounter
		'''	
		num_states = len(range_of_states)

		which_state = []
		for num in angles:
			categorized = False
			for i in range(0, num_states):
				if range_of_states == (None, None):
					continue
				if range_of_states[i][0] <= num <= range_of_states[i][1]:
					which_state.append(i)
					categorized = True
			if not categorized:
				which_state.append(-1)


		transition_ctr = transitions.TransitionCounter(num_states)
		for i in range(1, len(which_state)):

			if which_state[i] != which_state[i-1]:
				# left which_state[i-1]
				# entered which_state[i]

				if which_state[i-1] != -1:
					transition_ctr.increment_out(which_state[i-1])

				if which_state[i] != -1:
					transition_ctr.increment_in(which_state[i])

		return transition_ctr


	@staticmethod
	def state_populations(angles, range_of_states):
		'''
		calculate the transitions between each state 
		Args:
			[float]: angles
			range_of_states: [[float, float]]; bounds of the peaks 

		Returns:
			[float]: populations of each state
		'''	
		num_states = len(range_of_states) 
		which_state = []
		for num in angles:
			categorized = False
			for i in range(0, num_states):
				if range_of_states[i][0] <= num <= range_of_states[i][1]:
					which_state.append(i)
					categorized = True
			if not categorized:
				which_state.append(-1)

		pops = {i: 0 for i in range(-1, num_states)}

		for s in which_state:
			pops[s] += 1

		for s in range(num_states):
			pops[s] /= len(angles)

		return pops




	@staticmethod
	def transition_matrix(angles, range_of_states):
		'''
		count the transitions in and out each state (matrix)
		Args:
			[float]: angles
			range_of_states: [[float, float]]; bounds of the peaks 
		Returns:
			TransitionMatrixCounter 
		'''

		num_states = len(range_of_states)

		which_state = []


		for num in angles:
			categorized = False
			for i in range(0, num_states):
				if range_of_states == (None, None):
					continue
				if range_of_states[i][0] <= num <= range_of_states[i][1]:
					which_state.append(i)
					categorized = True
			if not categorized:
				which_state.append(-1)

		transition_ctr = transitions.TransitionMatrixCounter(num_states)

		for i in range(1, len(which_state)):
			if which_state[i] != which_state[i-1]:
				# left which_state[i-1]
				# entered which_state[i]
				transition_ctr.increment_transition(which_state[i-1], which_state[i])
		return transition_ctr


	# def check_transitions(self, transition_counter, min_transitions):
	# 	raise utils.NotImplementedError
	# 	# below is broken
	# 	# if len(transition_matrix) == 1:
	# 	#	 return
	# 	# for i,in_out in enumerate(transition_matrix):
	# 	#	 if in_out[0] < min_transitions:
	# 	#		 warnings.warn(f"State {i} has fewer that {min_transitions} in")
	# 	#	 if in_out[1] < min_transitions:
	# 	#		 warnings.warn(f"State {i} has fewer that {min_transitions} out")


	######################
	# PLOTTING FUNCTIONS #
	######################

	# def highlight_dihedral(self, dihedral, save_path=None):
	# 	pass

	@staticmethod
	def plot_dihedral_scatter(angles, ax=None, angle_min=None, title=None, show=False, save_path=None):

		if not ax:
			f, ax = plt.subplots()
		ax.scatter(np.arange(len(angles)), angles)
		ax.set_ylabel("Dihedral Angle (˚) --")
		ax.set_xlabel("Frame")
		ax.set_ylim([angle_min,angle_min+360])
		if title: 
			ax.set_title(title)
		if show:
			plt.show()
		if save_path:
			plt.savefig(save_path)


	def plot_kde(self, torsion, smoothing_window=-1, angle_min=None, save_path=None):

		f, ax = plt.subplots()

		angles = self.get_torsion_angles(torsion)
		X, score, angle_min = TorsionFinder.get_kde(angles, num_bins=40, angle_min=None)

		new_X, new_score = X,score

		if smoothing_window != -1:
			moving_avg = []
			for i in range(len(new_score) - smoothing_window):
				moving_avg.append(np.mean(new_score[i:i+smoothing_window]))
			new_score = moving_avg
			new_X = np.arange(angle_min, angle_min+360, 360/len(new_score))

		ax.plot(new_X,new_score)
		ax.set_xlabel("Dihedral Angle (˚)")
		ax.set_ylabel("KDE Score")
		ax.set_xlim([angle_min,angle_min+360])
		if save_path:
			plt.savefig(save_path, dpi=500)

	@staticmethod
	def plot_dihedral_histogram(angles, angle_min=None, pdf_individual=[], ax=None, title=None, num_bins=40, alpha=0.5, show=True, color=None, pdf_colors=[], save_path=None):

		if not ax:
			f, ax = plt.subplots()

		angle_min, angles = TorsionFinder.shift_torsion_angles(angles, num_bins, angle_min=angle_min)

		X = np.array(angles).flatten()
		kwargs = dict()
		if color != None:
			kwargs['color']=color
		ax.hist(X, bins=num_bins, range=(angle_min,angle_min+360), alpha=alpha, density=True, stacked=True, **kwargs)
		ax.set_xticks(np.arange(round(angle_min,-1), round(angle_min, -1)+361, 30))
		ax.set_xlabel(f"Dihedral Angle (˚) shifted by {int(angle_min) - (-180)}˚ to right")
		ax.set_ylabel("Frequency")

		pdf_individual_sort = sorted(pdf_individual, key=len)
		if pdf_colors == []:
			pdf_colors = ['red', 'orange', 'green', 'blue', 'purple', 'brown']

		for i,pdf in enumerate(pdf_individual_sort):
			if len(pdf.flatten()) != 0:
				ax.plot(pdf[0], pdf[1], color=pdf_colors[i])

		if title:
			ax.set_title(title)
		if show:
			plt.show()
		if save_path:
			plt.savefig(save_path, dpi=500)

	@staticmethod
	def plot_transition_counts(transition_ctr, ax=None, colors=[]):

		if not ax:
			f,ax = plt.subplots()

		df = transition_ctr.to_pandas()

		# df = pd.DataFrame(data=transition_matrix, columns=['transitions in', 'transitions out'])

		table = ax.table(
			cellText=df.values, 
			colLabels=[f"s{d} in" for d in df.columns], 
			rowLabels=[f"s{d} out" for d in df.index],
			loc='center',
		)
		table.auto_set_font_size(False)
		table.scale(0.93, 3)

		for (row, col), cell in table.get_celld().items():
			if (row == 0) or (col == -1):
				cell.set_text_props(ha="center", weight='bold', fontsize=18)
			else:
				cell.set_text_props(ha="center", fontsize=18)

			alpha_value = 0.5
			if len(colors) >= len(df) and col == -1:
				cell.set_facecolor(colors[row-1])
				cell.set_alpha(alpha_value)


			if len(colors) >= len(df) and row == 0:
				cell.set_facecolor(colors[col])
				cell.set_alpha(alpha_value)

			if (row == len(df) and col == -1) or (col == len(df)-1 and row == 0):
			 	cell.set_facecolor('lightgrey')

			if row == col + 1:
			 	cell.set_facecolor('lightgrey')

		ax.axis('off')

	def get_torsion_name(self, torsion):
		'''
		count the transitions in and out each state (matrix)
		Args:
			[float]: angles
			range_of_states: [[float, float]]; bounds of the peaks 
		Returns:
			str: string name of torsion based on residue + residue number and indices
		'''

		d1,d2,d3,d4 = tuple(torsion)

		sel_a_in_dih = self.mda_universe.select_atoms(f"index {torsion[0]}")
		sel_resid = sel_a_in_dih[0].residue
		return f"{sel_resid.resname}  {sel_resid.resid}"


	def make_torsion_img(self, torsion, angle_min=None, save_path=None):
		'''
		create a image of all the torsion transitions
		Args:
			torsion: [int, int, int, int]; indices of the torsion of interest
			angle_min: int; minimum angle of the histogram  
			save_path: str; path to save the image
		Returns:
			int: minimum angle fo the histogram
		'''

		d1,d2,d3,d4 = tuple(torsion)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    

		sel_a_in_dih = self.mda_universe.select_atoms(f"index {torsion[0]}")
		sel_resid = sel_a_in_dih[0].residue


		f,ax = plt.subplots(1, 3, figsize=(25, 6.25))
		sup_title = f"{sel_resid.resname} {sel_resid.resid} ({d1},{d2},{d3},{d4})"
		f.suptitle(sup_title,fontsize=60)
		f.tight_layout(pad=3.5)


		angles = self.get_torsion_angles(torsion)
		X, scores, angle_min = TorsionFinder.get_kde(angles) 
		
		num_peaks, peaks = TorsionFinder.get_kde_num_peaks(X, scores, smoothing_window=100, peak_prominence=0.008)
		
		min_max = TorsionFinder.get_bounds_mindist(X, scores, num_peaks, peaks)

		angles = TorsionFinder.shift_torsion_angles(angles, angle_min=angle_min)[1].flatten()

		angles = self.get_torsion_angles(torsion)

		transition_ctr = TorsionFinder.transition_matrix(angles, min_max)

		pdf_individual = []
		for mm in min_max:
			gmm,x,pdf,pdfi,bounds = TorsionFinder.get_individual_gmm(X, angle_min, mm)
			pdf_individual.append(pdfi)

		with tempfile.NamedTemporaryFile(suffix='.png') as highlightpng:
			self.highlight_dihedral(torsion, save_path=highlightpng)
			img = np.asarray(Image.open(highlightpng.name))
			ax[0].imshow(img)
			ax[0].axis('off')

		states_list = [f"s{i}" for i in range(num_peaks)]

		TorsionFinder.plot_dihedral_histogram(angles, ax=ax[1], show=False, pdf_individual=pdf_individual, angle_min=angle_min)
		ax[1].legend(states_list)
		pdf_colors = ['red', 'orange', 'green', 'blue', 'purple', 'brown']
		TorsionFinder.plot_transition_counts(transition_ctr, ax=ax[2], colors=pdf_colors)

		if save_path:
			plt.savefig(save_path)

		plt.close()
		return angle_min




class ProteinTorsionFinder(TorsionFinder):
	''' Class for running analysis on protein sidechain torsions
	'''

	def __init__(self, trajf: str, topf: str, ligcode: str):
		TorsionFinder.__init__(self, trajf, topf)
		self.ligcode = ligcode
		self.aa_only = self.mda_universe.select_atoms("protein and (name N or name CA or name C or name O or name CB)")
		self.ligand = self.mda_universe.select_atoms(f'resname {self.ligcode}')

	def get_binding_residues(self, A_cutoff: float):
		''' 
		identify all residues within A_cutoff angstroms of the ligand
		Args:
			A_cutoff: float Angstrom cutoff from binding mode

		Returns: 
			[int]: list of interacting residue IDs
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
		identify the chi1 torsion atom indices of specified residue indices
		Args:
			resids: [int]; residue ids

		Returns: 
			[[int, int, int, int]]: atom indices of chi1 torsions
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


	@staticmethod
	def get_intersection(list_of_sets: list):
		''' 
			returns a set that is an intersection of all the sets in the list
			Args:
				list_of_sets: [set]; list of sets to intersect

			Returns: 
				set: set that represents the intersection

		'''
		if len(list_of_sets) == 0:
			return set()

		s1 = list_of_sets[0]

		for s in list_of_sets[1:]:
			s1 = s1.intersection(s)

		return s1


	def get_chi_x_residues(self, x, sel=None, a_cutoff=None):
		''' 
		identify the chi torsion atom indices of specified residue indices at the specified chi angle
		Args:
			x: int; which chi up to chi8
			sel: str; MDAnalysis selection language of residues of interest
			a_cutoff: float; Angstrom cutoff from ligand

		Returns: 
			[[int, int, int, int]]: atom indices of chi1 torsions
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

		return sorted(list(ProteinTorsionFinder.get_intersection(ags_res_set_list)))

	def get_residues_ags(self, resgrp):
		''' 
		given the residues in resgrp, create an atom group of all relevant atoms
		Args:
			resgrp: ResidueGroup; MDAnalysis residue Group
	
		Returns:
			AtomGroup: MDAnalysis atom group

		'''
		ags = self.mda_universe.select_atoms("")

		for r in resgrp:
			if type(r) == int:
				ags += self.mda_universe.select_atoms(f"resid {r}")
			elif type(r) == mda.Residue:
				ags += self.mda_universe.select_atoms(f"resid {r.resid}")
		return ags


	def get_chi_x_ags(self, x, resgrp):
		''' 
		returns a list of atom groups inolved in chi x anglebased on the residue group
		Args:
			resgrp: ResidueGroup; MDAnalysis residue Group
	
		Returns:
			[AtomGroup]: MDAnalysis atom group
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
		''' 
		returns a list of lists of atom indices (aid)
		
		Args:
			x: int; chi angle (ex. 2 for chi2)
			resgrp: ResidueGroup; residues of interest

		Returns
			[[int]]; atom indices for each chi X torsion of interest
		'''
		ag_list = self.get_chi_x_ags(x, resgrp)

		chi_aid = list()

		for ag in ag_list:
			chi_aid.append(np.array([a.index for a in ag]))

		chi_aid_np = np.array(chi_aid)
		return chi_aid_np.transpose()


	def get_chi_x_torsions(self, x, sel=None, a_cutoff=None):
		''' 
		returns a list of chi x torsions atom indices in the protein given the selection (sel)
		Args:
			x: int; chi angle (ex. 2 for chi2)
			sel: str; MDAnalysis selection string for residues of interest
			a_cutoff: float; Angstrom cutoff

		Returns
			[[int]]; atom indices for each chi X torsion of interest
		'''
		# list of the residues (not resid) that have a chiX torsion
		chi_x_reslst = self.get_chi_x_residues(x, sel, a_cutoff)

		if len(chi_x_reslst) == 0:
			return []
		chi_x_resgrp = mda.ResidueGroup(chi_x_reslst)

		resnames = list()
		for r in chi_x_resgrp:
			resnames.append(r.resname)

		return self.get_chi_x_aid(x, chi_x_resgrp)



	def save_traj_sel(self, sel, frames, save_path):
		''' 
		saves out the atoms selected in sel as a new file for the purpose of this program meant to save out a single amino acid as a pdb file

		Args:
			sel: str; selection string using MDAnalysis atom selection language
			frames:	(int,int); 2 membered tuple with the start and end (not inclusive) frame to save out [start, end)
			save_path: str; path to save 
		Returns:
			None
		'''
		ag = self.mda_universe.select_atoms(sel)

		start_incl, end_nincl = frames

		with mda.Writer(save_path, ag.n_atoms) as w:
			for ts in self.mda_universe.trajectory[start_incl:end_nincl]:
				w.write(ag)



	def highlight_dihedral(self, torsion, save_path=None):
		''' 
		highlights the dihedral of a protein sidechain
		* exports the protein residue as a single AA pdb
		* loads in the single protein into rdkit with and
		  without Hs
		* compares the 2 to get atom orderings via MCS

		Args:
			torsion: [int, int, int, int]; atom indices representing torsion of interest
			save_path: str; path to save image
		Returns:
			None
		'''
		sel_resid = self.get_residue_from_torsion(torsion)
		min_res_aidx = min([a.index for a in sel_resid.atoms])

		with tempfile.NamedTemporaryFile(suffix='.pdb') as aapdb:
			self.save_traj_sel(f'resid {sel_resid.resid}', (0,1), aapdb.name)

			rdmol_wH = rdw.load_rdmol_from_file(aapdb.name, removeHs=False)
			rdmol_woH = rdw.load_rdmol_from_file(aapdb.name, removeHs=True)

			mapping = mappings.map_mols(rdmol_wH, rdmol_woH)

			adj_dih = [idx - min_res_aidx for idx in dihedral]

			woH_adj_dih = mappings.convert_dihedral(mapping, adj_dih)

			rdw.highlight_dihedral(rdmol_woH, woH_adj_dih, save_path)


	def get_residue_from_torsion(self, torsion):
		'''
		get the residue the torsion atom indices belong to
		Args:
			torsion: [int, int, int, int] atom indices of atom involved in torsion

		Returns:
			Residue: MDAnalysis residue that torsion belongs to
		'''
		sel_atm_in_dih = self.mda_universe.select_atoms(f"index {torsion[0]}")
		return sel_atm_in_dih[0].residue


	def determine_chi_x(self, torsion):
		raise utils.NotImplementedError




class LigandTorsionFinder(TorsionFinder):
	''' 
	Analyze torsions in the ligand of a system
	'''
	def __init__(self, trajf: str, topf: str, ligcode: str, smiles: str):
		''' 
		Args
			trajf:  str; simulation trajectory file
			topf:	str; topology file 
			ligcode: str; 3 letter ligand code
			smiles: str; smiles string for the ligand of interest

		'''

		TorsionFinder.__init__(self, trajf, topf)

		self.smiles = smiles
		self.ligcode = ligcode
		self.trajectory_len = len(self.mda_universe.trajectory)




		pdb_incorrect_atype = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)
		pdb_fixed_atype = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)
		selection = f"resname {self.ligcode}"


		self.rdmol = None
		frame = 0
		while not self.rdmol:
			try:
				self.export_pdb_from_traj(frame, pdb_incorrect_atype.name, sel=selection)
				pdbw.rename_lig_pdb_atoms(pdb_incorrect_atype.name, pdb_fixed_atype.name)

				pdbw.rename_lig_pdb_atoms(pdb_incorrect_atype.name, "pdb_mol.pdb")

				# begin old:
				self.rdmol_unsanitized = rdw.assign_bond_order_from_smiles(smiles, pdb_fixed_atype.name)
				self.rdmol = rdw.sanitize_rdmol(Chem.Mol(self.rdmol_unsanitized))
				self.oemol = mc.get_oemol_from_rdmol(self.rdmol)

			except AtomValenceException:
				frame += 1
				continue

		os.unlink(pdb_fixed_atype.name)
		os.unlink(pdb_incorrect_atype.name)


	def _check_top_has_conect(self, topf):

		if not topf.endswith(".pdb"):
			return False

		seen_lig_conect = False

		selection = f"resname {self.ligcode}"

		ags = self.mda_universe.select_atoms(selection)
		indices = [ a.index for a in ags ]

		with open(topf, 'r') as t:
			for line in t:
				if line.strip() == "END":
					return False
				if line.startswith("CONECT"):
					if int(line[6:11]) in indices:
						return True
		return False

	def get_rdmol(self):
		'''
		return RDKit Mol for the small molecule 
		Args:
			None
		Returns:
			Mol
		'''
		return self.rdmol

	def get_rdmol_by_frame(self, frame):
		'''
		returns the conformation of the rdkit molecule in the specified frame
		Args:
			frame: int; frame of trajectory of interest
		Returns:
			Mol
		'''
		pdb_incorrect_atype = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)
		pdb_fixed_atype = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)
		selection = f"resname {self.ligcode}"
		self.export_pdb_from_traj(frame, pdb_incorrect_atype.name, sel=selection)
		pdbw.rename_lig_pdb_atoms(pdb_incorrect_atype.name, pdb_fixed_atype.name)
		rdmol_unsanitized = rdw.assign_bond_order_from_smiles(self.smiles, pdb_fixed_atype.name)
		rdmol = rdw.sanitize_rdmol(Chem.Mol(rdmol_unsanitized))

		os.unlink(pdb_fixed_atype.name)
		os.unlink(pdb_incorrect_atype.name)
		return rdmol


	def convert_idx_to_sysidx(self, idx):
		''' 
		when ligand is not the only molecule in the system converts the ligand indices starting at 0 to indices that are relevant to the system
		Args:
			idx: int; index to be converted
		Returns:
			int: index of the atom within the context of the system
		'''
		selection = f"resname {self.ligcode}"
		ags = self.mda_universe.select_atoms(selection)
		indices = [ a.index for a in ags ]
		return int(min(indices) + idx)


	def export_pdb_from_traj(self, frame, opdb, sel="all"):
		''' 
		exports the specified selection (sel) from the frame of interest as a pdb
		Args:
			frame: int; frame of trajectory
			opdb: str; name of the output pdb
			sel; str; MD Analysis selection language of what to output
		Returns:
			None
		'''

		atms = self.mda_universe.select_atoms(sel)
		atms.write(opdb, frames=self.mda_universe.trajectory[[frame,]], bonds='conect')
	
	def _get_torsion(self, bond):
		''' 
		gets the atom indices of the torsion from the rotatable bond of interest
		Args:
			frame: int; frame of trajectory
			opdb: str; name of the output pdb
			sel; str; MD Analysis selection language of what to output
		Returns:
			[int, int, int, int]: atom indices of torsion
		'''
		# credit to: Travis Dabbous

		pos_1 = []
		pos_2 = []
		pos_3 = []
		pos_4 = []
		
		bgn_atom = bond.GetBeginAtom()
		pos_2 = bgn_atom.GetIdx()
		end_atom = bond.GetEndAtom()
		pos_3 = end_atom.GetIdx()
		#find the neighbors of the atoms from their lists
		nbors_pos2 = []
		
		for atom in bgn_atom.GetNeighbors():
			if bgn_atom.GetAtomicNum() != 6:
				nbors_pos2.append(atom.GetIdx())
			elif atom.GetAtomicNum() != 1:
				nbors_pos2.append(atom.GetIdx())

		# sorting atom indices for nbors_pos2 so lowest index is always chosen
		nbors_pos2 = sorted(nbors_pos2)
				
		for atom in nbors_pos2:
			if atom != pos_3:
				pos_1 = atom 
				break
				
		nbors_pos3 = []
		
		for atom in end_atom.GetNeighbors():
			if end_atom.GetAtomicNum() != 6:
				nbors_pos3.append(atom.GetIdx())
			elif atom.GetAtomicNum() != 1:
				nbors_pos3.append(atom.GetIdx())

		# sorting atom indices for nbors_pos3 so lowest index is always chosen
		nbors_pos3 = sorted(nbors_pos3)
	   
		for atom in nbors_pos3:
			if atom != pos_2:
				pos_4 = atom
				break

		if any([pos_1, pos_2, pos_3, pos_4]) == []:
			raise BadTorsionError

		# returned as a list to be able to properly index the atoms
		# of the mda universe
		return [pos_1, pos_2, pos_3, pos_4]

	def get_torsions(self):
		''' 
		gets a list of all torsions
		Args:
			None
		Returns:
			[[int, int, int, int]]: list of torsions
		'''
		# only get torsions for bonds that are rotatable 
		# rotatable bonds cannot be terminal
		torsions = []
		for bond in rdw.get_rotatable_bonds(self.rdmol):
			try:
				torsion = self._get_torsion(bond)
				torsions.append(torsion)
			except BadTorsionError:
				pass
		return torsions

	def highlight_dihedral(self, dihedral, save_path=None):
		''' highlights the dihedral of a protein sidechain
			* exports the protein residue as a single AA pdb
			* loads in the single protein into rdkit with and
			  without Hs
			* compares the 2 to get atom orderings via MCS
		'''
		rdw.highlight_dihedral(self.rdmol_unsanitized, dihedral, save_path)


	def make_torsion_img_no_shift(self, torsion, save_path=None):
		''' 
		create a matplot lib figure where the histogram is not shifted (bounds are -180 to 180)
		Args:
			torsion: [int, int, int, int]: torsions atom indices
			save_path: str; path to save image to 
		Returns:
			None
		'''

		angle_min = -180
		# prevents shifting

		d1,d2,d3,d4 = tuple(torsion)
		torsion_sys = [self.convert_idx_to_sysidx(i) for i in torsion]
		sel_a_in_dih = self.mda_universe.select_atoms(f"index {torsion_sys[0]}")
		sel_resid = sel_a_in_dih[0].residue

		f,ax = plt.subplots(1, 2, figsize=(12, 6.25))

		sup_title = f"{sel_resid.resname} {sel_resid.resid} ({d1},{d2},{d3},{d4})"

		f.suptitle(sup_title,fontsize=60)

		f.tight_layout(pad=3.5)

		with tempfile.NamedTemporaryFile(suffix='.png') as highlightpng:
			self.highlight_dihedral(torsion, save_path=highlightpng)
			img = np.asarray(Image.open(highlightpng.name))
			ax[0].imshow(img)
			ax[0].axis('off')

		angles = self.get_torsion_angles(torsion_sys)

		angles = TorsionFinder.shift_torsion_angles(angles, angle_min=angle_min)[1].flatten()
		TorsionFinder.plot_dihedral_histogram(angles, ax=ax[1], show=False,  angle_min=angle_min)

		if save_path:
			plt.savefig(save_path)



	def plot_dihedral_scatter(self, torsion, ax=None, angle_min=None, title=None, save_path=None):
		''' 
		create a matplot lib figure of scatter plot
		Args:
			torsion: [int, int, int, int]: torsions atom indices
			ax: axis; matplotlib axis to use for plot
			angle_min: int; minimum angle of the histogram
			title: str; title of the figure
			save_path: str; path to save image to 
		Returns:
			None
		'''
		if not ax:
			f, ax = plt.subplots()

		d1,d2,d3,d4 = tuple(torsion)

		torsion_sys = [self.convert_idx_to_sysidx(i) for i in torsion]
		angles = self.get_torsion_angles(torsion)

		angle_min, angles = TorsionFinder.shift_torsion_angles(angles, angle_min=angle_min)
		angles = angles.flatten()

		ax.scatter(np.arange(len(angles)), angles)
		ax.set_ylabel("Dihedral Angle (˚) --")
		ax.set_xlabel("Frame")
		ax.set_ylim([angle_min,angle_min+360])
		if title: 
			ax.set_title(title)
		if save_path:
			plt.savefig(save_path)




	def make_torsion_img(self, torsion, angle_min=None,  save_path=None):
		''' 
		create a matplot lib figure of scatter plot
		Args:
			torsion: [int, int, int, int]: torsions atom indices
			angle_min: int; minimum angle of the histogram
			save_path: str; path to save image to 
		Returns:
			[[float, float]]: torsion state boundaries
			TransitionMatrixCounter: counts of transitions between states
			[str]: names of each state
			bool: True if torison has symmetry, False otherwise
			float: populations of each state

		'''

		d1,d2,d3,d4 = tuple(torsion)

		torsion_sys = [self.convert_idx_to_sysidx(i) for i in torsion]

		sel_a_in_dih = self.mda_universe.select_atoms(f"index {torsion_sys[0]}")
		sel_resid = sel_a_in_dih[0].residue

		f,ax = plt.subplots(1, 3, figsize=(25, 6.25))
		sup_title = f"{sel_resid.resname} {sel_resid.resid} ({d1},{d2},{d3},{d4})"
		f.suptitle(sup_title,fontsize=60)
		f.tight_layout(pad=3.5)

		angles = self.get_torsion_angles(torsion_sys)

		X, scores, angle_min = TorsionFinder.get_kde(angles, angle_min=angle_min)

		num_peaks, peaks = TorsionFinder.get_kde_num_peaks(X, scores, smoothing_window=30, peak_prominence=0.01)

		min_max = TorsionFinder.get_bounds_mindist(X, scores, num_peaks, peaks)

		angle_min, angles = TorsionFinder.shift_torsion_angles(angles, angle_min=angle_min)



		transition_ctr = TorsionFinder.transition_matrix(angles, min_max)

		pdf_individual = []
		for mm in min_max:
			gmm,x,pdf,pdfi,bounds = TorsionFinder.get_individual_gmm(X, angle_min, mm)
			pdf_individual.append(pdfi)

		with tempfile.NamedTemporaryFile(suffix='.png') as highlightpng:
			self.highlight_dihedral(torsion, save_path=highlightpng)
			img = np.asarray(Image.open(highlightpng.name))
			ax[0].imshow(img)
			ax[0].axis('off')

		states_list = [f"s{i}" for i in range(num_peaks)]

		pdf_colors = ['red', 'orange', 'green', 'blue', 'purple', 'brown']

		TorsionFinder.plot_dihedral_histogram(angles, ax=ax[1], show=False, angle_min=angle_min, pdf_individual=pdf_individual,pdf_colors = pdf_colors)

		ax[1].legend(states_list)
		g = self.get_statistical_inefficiency(angles)
		ax[1].set_title(f'statistical inefficiency, g = {round(g,2)}', fontsize=15, color='black')


		
		TorsionFinder.plot_transition_counts(transition_ctr, ax=ax[2], colors=pdf_colors)

		symmetry = False
		if mappings.check_symmetry(self.rdmol, torsion):
			# if there is symmetry add a note on the image
			symmetry = True
			ax[2].text(0.01, 0.01, f'** Warning: torsion has symmetry, disregard transitions', fontsize=15, color='red')

		populations = TorsionFinder.state_populations(angles, min_max)

		if save_path:
			plt.savefig(save_path)
		return min_max, transition_ctr, states_list, symmetry, populations




