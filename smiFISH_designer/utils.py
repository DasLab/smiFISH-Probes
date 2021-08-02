def read_fasta(fasta_path):

    with open(fasta_path) as f:
        name = f.readline().strip().replace(' ', '').split('>')[-1]
        sequence = f.readline().strip()

    return [name, sequence]

def write_probes_fasta_from_dataframe(df, fasta_path):

    with open(fasta_path, 'w') as f:
        for index, row in df.iterrows():
            f.write('> {} {}-{}\n'.format(row['final_sequence_name'], row['start'], row['end']))
            f.write(row['probe_sequence'])
            f.write('\n')
