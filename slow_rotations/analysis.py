

LOW_TRANSITION_THRESHOLD = 10

def check_transitions(torsionrepeat):
	'''
	Checks the torsion for missing states
	
	Args:
        torsionrepeat: TorsionRepeat

    Returns:
        [int]: list of states that are missing
	'''

	transitions = torsionrepeat.transitions
	for s1,states in transitions.items():
		for s2,cts in states.items():
			if cts < LOW_TRANSITION_THRESHOLD and 'Ø' not in s1 and 'Ø' not in s1:
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