import os

list_of_sequences = ['P59594', 'P0DTC2', 'K9N5Q8', 'P36334', 'Q0ZME7', 'P15423', 'Q6Q1S2', 'Q5MQD0', 'Q14EB0']


def alignment(seq_list):
    seq_string = ''
    for i in seq_list:
        seq_string = seq_string + 'sp:' + i + ','
    seq_string = seq_string[:-1]
    #os.system(
    #    'python msa_algos/clustalo.py --email erald.bb@gmail.com --sequence ' + seq_string + ' --outfile msa_results/clustal_omega/clustal_spike  --outfmt fa')
    #os.system(
    #    'python msa_algos/kalign.py --email erald.bb@gmail.com --stype protein --sequence ' + seq_string + ' --outfile msa_results/kalign/kalign_spike --format fasta')
    os.system(
        'python msa_algos/mafft.py --email erald.bb@gmail.com --stype protein --sequence ' + seq_string + ' --outfile msa_results/mafft/mafft_spike --format fasta')
    os.system(
        'python msa_algos/muscle.py --email erald.bb@gmail.com --sequence ' + seq_string + ' --outfile msa_results/muscle/muscle_spike --format fasta')
    #os.system(
    #    'python msa_algos/tcoffee.py --email erald.bb@gmail.com --stype protein --sequence ' + seq_string + ' --outfile msa_results/tcoffee/tcoffee_spike --format fasta_aln')


alignment(list_of_sequences)
