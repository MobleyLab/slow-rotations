import pandas as pd

class TransitionCounter():
	def __init__(self, num_states):
		self.transition_dict = dict()
		for s in range(num_states):
			self.transition_dict[s] = {"in": 0, "out": 0}
		self.num_states = num_states

	def increment_out(self, state):
		self.transition_dict[state]["out"] += 1

	def increment_in(self, state):
		self.transition_dict[state]["in"] += 1

	def to_dict(self):
		return self.transition_dict

	def to_pandas(self):
		return pd.DataFrame.from_dict(self.transition_dict).transpose()

	def min_transitions(self):
		minimum = 9999999999999999
		for s,ct in self.transition_dict.items():
			if ct['in'] < minimum:
				minimum = ct['in']
		return minimum

	def max_transitions(self):
		maximum = -99999999999999
		for s,ct in self.transition_dict.items():
			if ct['in'] > maximum:
				maximum = ct['in']
		return maximum

	def total_transitions(self):
		total = 0
		for s,ct in self.transition_dict.items():
			total += ct['in']

		return total

	def num_states(self):
		return self.num_states

class TransitionMatrixCounter():
	def __init__(self, num_states):
		self.transition_dict = dict()
		for s in range(num_states):
			self.transition_dict[s] = {i: 0 for i in range(num_states)}
			self.transition_dict[s]['Ø'] = 0
		self.transition_dict['Ø'] = {i: 0 for i in range(num_states)}
		self.transition_dict['Ø']['Ø'] = 0
		self.num_states = num_states

	def increment_transition(self, state1, state2):
		# transitions from state1 into state2

		if state1 == -1:
			self.transition_dict['Ø'][state2] += 1
		elif state2 == -1:
			self.transition_dict[state1]['Ø'] += 1
		else:
			self.transition_dict[state1][state2] += 1

	def to_dict(self):
		return self.transition_dict

	def to_pandas(self):
		return pd.DataFrame.from_dict(self.transition_dict).transpose()

	def min_transitions(self):
		# note this will mark everything with 1 real state and 1 
		# catch all state as the default minimum
		minimum = 9999999999999
		for s1 in self.transition_dict.keys():
			for s2 in self.transition_dict[s1].keys():
				if s1 == s2: # or s1 == 'Ø' or s2 == 'Ø':
					continue
				if self.transition_dict[s1][s2] < minimum:
					minimum = self.transition_dict[s1][s2]

		return minimum


	def count_transitions_into_state(self, state):
		ct = 0
		for i in range(self.num_states): 
			ct += self.transition_dict[i][state]

		ct += self.transition_dict['Ø'][state]
		return ct


	def count_transitions_out_of_state(self, state):
		ct = 0
		for i in range(self.num_states): 
			ct += self.transition_dict[state][i]

		ct += self.transition_dict[state]['Ø']
		return ct


	def max_transitions(self):
		# note this will mark everything with 1 real state and 1 
		# catch all state as the default maximum

		maximum = -99999999999999
		for s1 in self.transition_dict.keys():
			for s2 in self.transition_dict[s1].keys():
				if s1 == s2 or s1 == 'Ø' or s2 == 'Ø':
					continue

				if self.transition_dict[s1][s2] > maximum:
					maximum = self.transition_dict[s1][s2]		
		return maximum


	def total_transitions(self):
		total = 0
		for s1 in self.transition_dict.keys():
			for s2 in self.transition_dict[s1].keys():
				total += self.transition_dict[s1][s2]		

		return total

	def num_states(self):
		return self.num_states

