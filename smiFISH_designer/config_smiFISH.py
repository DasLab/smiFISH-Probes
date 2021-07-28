from configparser import ConfigParser

config = ConfigParser()

# Section required for running smiFISH probe designer
config['Required'] = {
	# Name of target RNA
	'Target RNA Name' : 'MBP'

	# Sequence of target RNA
	'Target RNA Sequence' : ''

	# Name of FLAP sequence
	'FLAP Name' : 'X'

	#Sequence of FLAP
	'FLAP Sequence' : 'CCTCCTAAGTTTCGAGCTGGACTCAGTG'

	# Desired path for output files
	'Output' : '../output'
}

# Section optional for running smiFISH probe designer
config['Optional'] = {
	# Minimumm length of desired probes (without FLAP sequence)
	'Probe Min' : '20'

	# Maximum length of desired probes (without FLAP sequence)
	'Probe Max' : '30'

	# Maximum length of acceptable probe length. Based on the company you are ordering probes with (with FLAP sequence)
	'Acceptable Probe Max' : '60'
}