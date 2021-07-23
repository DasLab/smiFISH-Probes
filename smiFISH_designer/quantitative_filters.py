import pandas as pd
import numpy as np

# Calculate deltaG using Nearest-Neighbor Parameters
def deltaG1(probe_sequence):
	table = {"AA":0.2, "AC":-1.4, "AG":-0.4, "AT":-0.4,
	"CA":-1.6, "CC":-2.3, "CG":-1.4, "CT":-1.3,
	"GA":-1.4, "GC":-2.0, "GG":-1.7, "GT":-1.5, 
	"TA":-0.5, "TC":-1.5, "TG":-1.2, "TT":-0.7}

	for p in range(len(probe_sequence)):
		# Add deltaG of pairing, then next pairing including the last base-pair (ex. ACGA = AC, CG, GA)
		list_of_pairs = [probe_sequence[p:2 + p] for p in range(0, len(probe_sequence))]
		# Delete last element of list, not a paired base
		list_of_pairs.pop()
		if len(list_of_pairs) > 0:
			# Make a series using the dictonary of calculated deltaG's
			C = pd.Series(list_of_pairs).map(table)
			D = list(C)
			return sum(D)
		else:
			return 'Error'

# Calculate GC percentage
def GC_percent(probe_sequence):
	G = probe_sequence.count('G')
	C = probe_sequence.count('C')
	percent_GC = (G + C)/len(probe_sequence)
	return percent_GC

# Filter for GC content between 40% and 60%
def GC_filter(probe_sequence):
    return 0.40 <= GC_percent(probe_sequence) <= 0.60

# Find start index of probe within complementary strand
def start_index(probe_sequence, comp):
	if probe_sequence in comp:
		start = comp.find(probe_sequence)
		return start

	return 'Error'

# Find end index of probe within complementary strand
def end_index(probe_sequence, comp):
	if probe_sequence in comp:
		end = comp.find(probe_sequence)
		end += len(probe_sequence) - 1
		return end

	return 'Error'

# Well position for ordering .xlsx
def well_position_list(n):
    well_y = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    newrow = [[i+str(j) for j in range(1,13)] for i in well_y]
    well_name = [j for i in newrow for j in i]
    well_name_array = np.array(well_name)
    
    well_names_output = []
    for p in range(n):
        well_names_output.append(well_name_array[p % len(well_name_array)])

    return well_names_output