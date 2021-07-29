import os

import pandas as pd
import click

from filters import pnas_filters, overlap_filters, quantitative_filters
from utils import read_fasta, write_probes_fasta_from_dataframe
 
@click.command()
@click.option('--fasta', type=str, required=True, help='Path to fasta file of target RNA sequence')
@click.option('--output', type=str, required=True, help='Path to output folder')
# @click.STRING('--FLAP_name', required=True, help='Name of FLAP sequence')
# @click.STRING('--FLAP_sequence', required=True, help='Name of FLAP sequence')

@click.option('--acceptable_max', type=int, required=False, default=60, help='Maximum length of acceptable probe length. Based on the company you are ordering probes with (with FLAP sequence)')
@click.option('--probe_min', type=int, required=False, default=26, help='Minimum length of desired probe (without FLAP sequence)')
@click.option('--probe_max', type=int, required=False, default=32, help='Maximum length of desired probes (without FLAP sequence')


def main(fasta, output, probe_min, acceptable_max, probe_max):

	# Read in target fasta file
	sequence_name, target_sequence = read_fasta(fasta)

	# List of GC content filter
	gc_filter = ['GC_filter']

	# List of 5 PNAS rules for optimum probe performance
	filters = ['pnas1', 'pnas2', 'pnas3', 'pnas4', 'pnas5']

	# Flip the inputed target RNA and replace bases to create complement strand
	reverse = target_sequence[::-1]
	comp = target_sequence.lower()[::-1].replace("g", "C").replace("c", "G").replace("t", "A").replace("a", "T").replace("u", "A")

	# Split complement strand (comp) into all possible probes between range (ex. 20-30 bases long)
	split_probes = []
	for n in range(probe_min, probe_max): 
		for i in range(0, len(comp), n):
			probe = comp[i:i + n]

			# If length of probe is shorter than minimum range, remove from list
			if len(probe) > 20:
				split_probes.append(probe)


	# Main DataFrame
	df = pd.DataFrame()

	# Inputed names/sequences for FLAP and target RNA
	target_name = sequence_name
	FLAP_name = 'X'
	FLAP_sequence = 'CCTCCTAAGTTTCGAGCTGGACTCAGTG'
	final_sequence = [probe + FLAP_sequence for probe in split_probes]

	# Adding columns to main DataFrame (boolean --> integers)
	df['probe_sequence'] = split_probes
	df['probe_length'] = df['probe_sequence'].apply(len)
	df['start'] = df['probe_sequence'].apply(quantitative_filters.start_index, comp=comp)
	df['end'] = df['probe_sequence'].apply(quantitative_filters.end_index, comp=comp)
	df['pnas1'] = df['probe_sequence'].apply(pnas_filters.pnas_filter_1).astype(int)
	df['pnas2'] = df['probe_sequence'].apply(pnas_filters.pnas_filter_2).astype(int)
	df['pnas3'] = df['probe_sequence'].apply(pnas_filters.pnas_filter_3).astype(int)
	df['pnas4'] = df['probe_sequence'].apply(pnas_filters.pnas_filter_4).astype(int)
	df['pnas5'] = df['probe_sequence'].apply(pnas_filters.pnas_filter_5).astype(int)
	df['deltaG'] = df['probe_sequence'].apply(quantitative_filters.deltaG1)
	df['GC'] = df['probe_sequence'].apply(quantitative_filters.GC_percent)
	df['GC_filter'] = df['probe_sequence'].apply(quantitative_filters.GC_filter).astype(int)
	df['target_name'] = target_name
	df['FLAP_name'] = FLAP_name
	df['FLAP_sequence'] = FLAP_sequence
	df['final_sequence'] = final_sequence

	# Create DataFrame with GC filter and using filters found in overlap_filters file
	fdf = overlap_filters.iteratively_find_probe_set(df, filters, gc_filter)

	# Final sequence name (target name + FLAP name + ranking number)
	fdf_index_reset = fdf.reset_index()
	fdf_index = fdf_index_reset.index.map(str)
	final_sequence_name = (target_name + '-' + FLAP_name + '-' + fdf_index)

	# Output of final, filtered DataFrame
	fdf['final_sequence_name'] = final_sequence_name
	fdf['well_position'] = quantitative_filters.well_position_list(len(fdf))
	print(fdf)

	output_file_path = os.path.join(output, sequence_name + '_smiFISH_probes')
	ordering_output_file_path = os.path.join(output, sequence_name + '_ordering_smiFISH_probes')

	'''
	Output excel with two sheets: 
	1. all columns of final DataFrame
	2. Ordering form
	'''
	fdf.to_excel(output_file_path + '.xlsx', sheet_name = 'All Data', index=False)
	fdf.to_excel(ordering_output_file_path + '.xlsx', columns=['well_position', 'final_sequence_name', 'final_sequence'], sheet_name = "Ordering Form", index=False)

	fdf.to_csv(output_file_path + '.csv', index=False)
	fdf.to_csv(ordering_output_file_path + '.csv', columns=['well_position', 'final_sequence_name', 'final_sequence'], index=False)

	"""
	Output a fasta file of probes for easy BLASTing
	"""
	write_probes_fasta_from_dataframe(fdf, os.path.join(output, sequence_name + '_filtered_probes.fa'))



if __name__ == "__main__":
	main()
