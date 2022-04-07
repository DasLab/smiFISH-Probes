"""
TODO: Explain what this file does overall
"""

import pandas as pd

# Takes in two intervals, and returns boolean if they overlap or not
# if there is desired spacing between intervals, return true only if they are at least n away from each other
# e.g. a spacing of 2 means that there need to be at least 2 nt between intervals
def overlapping_intervals(a, b, spacing=0):
    if a[0] < b[0]:
        first = a
        second = b
    elif a[0] == b[0]:
        return True
    else:
        first = b
        second = a

    return first[1] + spacing > second[0]


# Take a probe and see if it does or does not overlap with probes in a list of probes
# Note that probes are supplied as dictionaries with 'stert' and 'stop' coordinates
def check_probe_fits(probe_candidate, probe_list, spacing=0):
	for probe in probe_list:
		if overlapping_intervals([probe['start'], probe['end']], [probe_candidate['start'], probe_candidate['end']]):
			return False
	return True


# Given dataframe and list of filters, return filtered dataframes
# Filters are names (strings) that correspond to boolean columns in the pandas DataFrame
def filter_df(df, filters):
    return df[df[filters].all(True)]


def find_nonoverlapping_probes_around_dG_median(df, dG_max_distance=0.1):
    """
    We want our probes to all have similar dG of annealing at 37 C for all-or-nothing signal

    Here, we start at the probe with medium dG and iteratively accumulate the next closest non-overlapping
    probes until we exhaust all available probes.

    """
    # Calculate median delta G at 37 of all possible probes
    median = df['deltaG'].median()
    
    # Order rows in probes dataframe by distance to median with Pandas-function
    for probe in range(len(df['deltaG'])):
        result_index = df['deltaG'].sub(median).abs().argsort()[0:].tolist()
        ordered_rows = df.iloc[result_index]

    # Seed our probe list with the probe with median dG
    final_probe_list = [df.iloc[result_index[0]]]

    # March down the list of probes ordered by dG clostest to median
    # If the probe doesn't overlap with a probe already in our set, add it to the set
    for index in result_index[1:]:
        probe_candidate = df.iloc[index]
        if check_probe_fits(probe_candidate, final_probe_list):
            if median*(1 - dG_max_distance) > probe_candidate['deltaG'] > median*(1 + dG_max_distance):
                final_probe_list.append(probe_candidate)
    
    return pd.DataFrame(final_probe_list)



def iteratively_find_probe_set(df, filters, gc_filter, min_probe_count=20):
    """
    We want a set of probes that passes as many of our filters as possible

    Here, we first try to make a probe set with all filters on. If the resulting probe set is too small,
    we relax the least stringent filter and try again.

    TODO: RE-decide on if we want to remove the most or least stringent filters first.

    """

    # First, we count how many probes are filtered out by each of the filters
    filtered_out_counts = {}

    for f in filters:
        filtered_out_count = len(df) - df[f].sum()
        filtered_out_counts[f] = filtered_out_count

    # Order the filters by how many probes they each filter out
    ordered_filters = sorted(filtered_out_counts, key=filtered_out_counts.get, reverse=False)

    gc_filter_on = df[df[gc_filter].all(True)]
    
    median_and_filters = find_nonoverlapping_probes_around_dG_median(filter_df(gc_filter_on, filters))

    while len(median_and_filters) < min_probe_count:
        if len(median_and_filters) >= min_probe_count:
            break
        
        if len(ordered_filters) > 0:
            removed_filter = ordered_filters.pop()
            print('Removing filter ' + removed_filter)
            filtered_filters_df = filter_df(gc_filter_on, ordered_filters)
            median_and_filters = find_nonoverlapping_probes_around_dG_median(filter_df(gc_filter_on, ordered_filters))
        else:
            print('No more PNAS filters to remove...')
            quit()


        

    
    if len(median_and_filters) >= min_probe_count:
        print('Found {} probes using filters {}'.format(len(median_and_filters), ' '.join(ordered_filters)))
        return median_and_filters
    else:
        raise ValueError('Not enough probes even with all filters off')
