import pandas as pd

# Takes in two intervals, and returns boolean if they overlap or not
def overlapping_intervals(a, b):
    if a[0] < b[0]:
        first = a
        second = b
    elif a[0] == b[0]:
        return True
    else:
        first = b
        second = a
        
    return first[1] > second[0]

# Takes start and end range of probe and compares it to start and end of other, lower-ranking probes
def check_probe_fits(probe_candidate, probe_list):
	for probe in probe_list:
		if overlapping_intervals([probe['start'], probe['end']], [probe_candidate['start'], probe_candidate['end']]):
			return False
	return True

# Find probes that are nonoverlapping based on their ranking of distance away from the median
def find_nonoverlapping_probes_around_median(df):
    
    # Calculate median delta G
    median = df['deltaG'].median()
    
    # Order rows in probes dataframe by distance to median with Pandas-function
    for probe in range(len(df['deltaG'])):
        result_index = df['deltaG'].sub(median).abs().argsort()[0:].tolist()
        ordered_rows = df.iloc[result_index]

    # Iteratively build list of non-overlapping probes around the median
    final_probe_list = [df.iloc[result_index[0]]]
    # Probe vs. all other probe_candidates based on ranking
    for index in result_index[1:]:
        probe_candidate = df.iloc[index]
        if check_probe_fits(probe_candidate, final_probe_list):
            final_probe_list.append(probe_candidate)
    
    return pd.DataFrame(final_probe_list)



# Given dataframe and filters, return filtered dataframes
def filter_df(df, filters):
    return df[df[filters].all(True)]

# NOTE: Requiers 
def iteratively_find_probe_set(df, filters, gc_filter, min_probes):
    
    filtered_out_counts = {}

    for f in filters:
        filtered_out_count = len(df) - df[f].sum()
        filtered_out_counts[f] = filtered_out_count

    ordered_filters = sorted(filtered_out_counts, key=filtered_out_counts.get, reverse=False)
    
    gc_filter_on = df[df[gc_filter].all(True)]
    
    median_and_filters = find_nonoverlapping_probes_around_median(filter_df(gc_filter_on, filters))

    while len(median_and_filters) < min_probes:
        ordered_filters.pop(0)
        filtered_filters_df = filter_df(gc_filter_on, ordered_filters)
        median_and_filters = find_nonoverlapping_probes_around_median(filter_df(gc_filter_on, ordered_filters))
        
        if len(median_and_filters) >= min_probes:
            break
    
    if len(median_and_filters) >= min_probes:
        print('Found {} probes using filters {}'.format(len(median_and_filters), ' '.join(ordered_filters)))
        return median_and_filters
    else:
        raise ValueError('Not enough probes even with all filters off')