import os

import pandas as pd
import click
from configparser import ConfigParser

from filters import pnas_filters, overlap_filters, quantitative_filters
from quantitative_filters import max_probe
from utils import read_fasta, write_probes_fasta_from_dataframe
 
@click.command()
# Must input either Fasta file or Config File
@click.option('--fasta', '-fa', type=str, required=False, help='Path to target sequence fasta file')
@click.option('--output', '-o', type=str, required=True, help='Path to output folder')
@click.option('--config', '-cfg' , type=str, help='Path to optional config file with parameters')

def main(fasta, config, output):
	"""
	There are two main ways to 
	"""
	print(fasta, config, output)
	if config is not None:
		parser = ConfigParser()
		parser.read(config)

		target_name = parser['REQUIRED']['TargetName']
		target_sequence = parser['REQUIRED']['TargetSequence']

		min_probe_number = int(parser['OPTIONAL']['MinProbeNumber'])
		probe_min = parser['OPTIONAL']['ProbeMin']
		probe_max = parser['OPTIONAL']['ProbeMax']
		AcceptableProbeMax = parser['OPTIONAL']['AcceptableProbeMax']
		FLAP_name = parser['OPTIONAL']['FLAPName']
		FLAP_sequence = parser['OPTIONAL']['FLAPSequence']

		
		if fasta is not None:
			print('Detected config file and inputted fasta. Ignoring inputted fasta file for config-defined sequence')

	else:
		print('Sit tight, will implement other functionality soon...')


	"""Parameters that are optional to change:
		probe minimum and maximum length
		acceptable probe length for ordering company
		FLAP name and sequence
	"""

	




	# if input_file.endswith('.fa'):
	# 	fasta = input_file
	# 	sequence_name, target_sequence = read_fasta(fasta)
	# elif input_file.endswith('.cfg'):
	# 	config_required = input_file
	# 	parser.read(config_required)
	# 	sequence_name = parser['REQUIRED']['TargetName']
	# 	target_sequence = parser['REQUIRED']['TargetSequence']
	# else:
	# 	print('Use Fasta or Config file as input!')


	# List of GC content filter
	gc_filter = ['GC_filter']

	# List of 5 PNAS rules for optimum probe performance
	filters = ['pnas1', 'pnas2', 'pnas3', 'pnas4', 'pnas5']

	# Flip the inputed target RNA and replace bases to create complement strand
	reverse = target_sequence[::-1]
	comp = target_sequence.lower()[::-1].replace("g", "C").replace("c", "G").replace("t", "A").replace("a", "T").replace("u", "A")

	# Split complement strand (comp) into all possible probes between range (ex. 20-30 bases long)
	split_probes = []
	for n in range(int(probe_min), int(probe_max)): 
		for i in range(0, len(comp), n):
			probe = comp[i:i + n]

			# If length of probe is shorter than minimum range, remove from list
			if len(probe) > 20:
				split_probes.append(probe)


	# Main DataFrame
	df = pd.DataFrame()

	# Inputed names/sequences for FLAP and target RNA
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
	fdf = overlap_filters.iteratively_find_probe_set(df, filters, gc_filter, min_probe_number)

	# Final sequence name (target name + FLAP name + ranking number)
	fdf_index_reset = fdf.reset_index()
	fdf_index = fdf_index_reset.index.map(str)
	final_sequence_name = (target_name + '-' + FLAP_name + '-' + fdf_index)

	# Output of final, filtered DataFrame
	fdf['final_sequence_name'] = final_sequence_name
	fdf['well_position'] = quantitative_filters.well_position_list(len(fdf))
	print(fdf)

	output_file_path = os.path.join(output, target_name + '_smiFISH_probes')
	ordering_output_file_path = os.path.join(output, target_name + '_ordering_smiFISH_probes')

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
	write_probes_fasta_from_dataframe(fdf, os.path.join(output, target_name + '_filtered_probes.fa'))



if __name__ == "__main__":
	main()
