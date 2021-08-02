"""
TODO Include info about the paper here
"""


# PNAS 1: Filter out nucleotide compositions with A < 28%
def pnas_filter_1(probe_sequence):
	probe_sequence = probe_sequence.upper()
	percent_a = probe_sequence.count('A')/len(probe_sequence)

	return 0.28 < percent_a


# PNAS 2: Remove anything with AAAA stacks
def pnas_filter_2(probe_sequence):
	probe_sequence = probe_sequence.upper()
	A_stack = 'AAAA'

	return A_stack not in probe_sequence


# PNAS 3: Enforce C composition of between 22-28%
def pnas_filter_3(probe_sequence):
	probe_sequence = probe_sequence.upper()
	percent_c = probe_sequence.count('C')/len(probe_sequence)

	return 0.22 < percent_c < 0.28


# PNAS 4: No CCCC stacks in any 6 consecutive nucleotides in the first 12 positions
def pnas_filter_4(probe_sequence):
	probe_sequence = probe_sequence.upper()
	C_stack = 'CCCC'
	sequence_subset = probe_sequence[0:12]

	return C_stack not in sequence_subset


# PNAS 5: No 4 nonconsecutive Cs in any 6 consecutive nucleotides in the first 12 positions
def pnas_filter_5(probe_sequence):
	probe_sequence = probe_sequence.upper()
	sequence_subset = probe_sequence[0:12]
	for p in range(6):
		if sequence_subset[p:p+6].count('C') >= 4:
			return False
	return True
