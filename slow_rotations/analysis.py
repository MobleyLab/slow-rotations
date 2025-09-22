
DEFAULT_LOW_TRANSITION_THRESHOLD = 10

def check_transitions(torsionrepeat, min_transitions=DEFAULT_LOW_TRANSITION_THRESHOLD):
	'''
	Checks the torsion for missing states
	
	Args:
        torsionrepeat: TorsionRepeat
        min_transitions: int; minimum transitions that must be observed in and out of a state to be considered adequately sampled

    Returns:
        [int]: list of states that are missing
	'''

	transitions = torsionrepeat.transitions
	for s1,states in transitions.items():
		for s2,cts in states.items():
			if cts < min_transitions and 'Ø' not in s1 and 'Ø' not in s1:
				return True
	return False




def check_states(torsionrepeat):
	'''
	Checks the torsion for missing states
	
	Args:
        torsionrepeat: TorsionRepeat

    Returns:
        [int]: list of states that are missing
	'''
	missing = []
	for state,is_missing in torsionrepeat.missing_states.items():
		if is_missing:
			missing.append(state)
	return missing