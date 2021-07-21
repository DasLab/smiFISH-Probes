# PNAS 1: Filter out nucleotide compositions with A < 28%
def pnas_filter_1(probe_sequence):
	percent_a = probe_sequence.count('A')/len(probe_sequence)
	return 0.28 < percent_a

# PNAS 2: Remove anything with AAAA stacks
def pnas_filter_2(probe_sequence):
	A_stack = 'AAAA'
	if A_stack in probe_sequence:
		return False 
	else:
		return True

# PNAS 3: Enforce C composition of between 22-28%
def pnas_filter_3(probe_sequence):
	percent_c = probe_sequence.count('C')/len(probe_sequence)
	return 0.22 < percent_c and 0.28 > percent_c

# PNAS 4: No CCCC stacks in any 6 consecutive nucleotides in the first 12 positions
def pnas_filter_4(probe_sequence):
	C_stack = 'CCCC'
	sequence_subset = probe_sequence[0:12]
	if C_stack in sequence_subset:
		return False
	else:
		return True 

# PNAS 5: No 4 nonconsecutive Cs in any 6 consecutive nucleotides in the first 12 positions
def pnas_filter_5(probe_sequence):
	sequence_subset = probe_sequence[0:12]
	for p in range(6):
		if sequence_subset[p:p+6].count('C') >= 4:
			return False
	return True