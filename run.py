from msa import alignment
from msa import get_conserved_sequences
from prediction import predict_all
from analyse import make_inputs_for_analysis
from analyse import analyse_all
from analyse import make_csv_from_results
from analyse import read_prediction_results



example_seq_dict = {'P0DTC2': ['SRSARSAIEDLLFDKVTIADPGYMQGYDDC','WSYTGSSFYAPEPITSLNTKY'],
'P36334': ['WMYTGSGYYYPEPITENNVVV','FKEELDQWFKNQTSVAPDL']}
list_of_swissprot_ids = ['P59594', 'P0DTC2', 'K9N5Q8', 'P36334', 'Q0ZME7', 'P15423', 'Q6Q1S2', 'Q5MQD0', 'Q14EB0']
path_to_mafft_alignment = '/home/erald/Desktop/ScrapyEpitope/msa_results/mafft/mafft.aln-fasta.fasta'
path_to_muscle_alignment = '/home/erald/Desktop/ScrapyEpitope/msa_results/muscle/muscle.aln-fasta.fasta'
mhci_alleles = ['HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*03:01']
mhci_lengths = [8, 9, 10]
mhcii_alleles = ['HLA-DRB1*01:01', 'HLA-DRB1*07:01', 'HLA-DRB1*03:01', 'HLA-DRB1*15:01', 'HLA-DRB3*01:01', 'HLA-DRB3*02:02']
mhcii_lengths = [11, 12, 13]
list_of_pdb_ids = ['6vxx']

# alignment(list_of_swissprot_ids, matrix='bl62', gapopen=1.53, gapext=0.123, order='aligned', nbtree=2, treeout='true', maxiterate=2, ffts='none')
# conserved_sequences_mafft = get_conserved_sequences(path_to_mafft_alignment, min_seq_conserved_pos='default', min_seq_flank_pos='default', max_contigous_nonconserved_pos = 8, min_length_block= 10, allowed_gap_pos='None')
predict_all(example_seq_dict, mhci_alleles, mhci_lengths, mhcii_alleles, mhcii_lengths, list_of_pdb_ids)
prediction_results = read_prediction_results()
analysis_input = make_inputs_for_analysis(prediction_results, list_of_swissprot_ids)
print("Analysing sequences...")
analysis_results = analyse_all(analysis_input)
make_csv_from_results(prediction_results, analysis_results)
print("Analysis finished. Sequences saved in 'results' folder in csv format")



def run_pipeline():

    """This function will run the pipeline. It first requires the user to input the swissprot and 
    pdb ids, then the mhc class I and II alleles and lengths. Then it runs the alignment algorithms 
    on the swissprot ids, asks the user for the path to the alingment file, thus the user can choose
    either the muscle or the mafft alignment, gets the conserved sequences, predicts the epitopes and
    then analyses all the epitopes and saves the results in csv files in folder 'results'."""

    swissprot_ids_str = input("Enter your Swissprot IDs separated by commas: ")
    pdb_ids_str = input("Enter your PDB IDs separated by commas: ")
    mhci_alleles_str = input("Enter your MHC class I alleles separated by commas: ")
    mhci_lengths_str = input("Enter your MHC class I lengths separated by commas (one length for each allele): ")
    mhcii_alleles_str = input("Enter your MHC class II alleles separated by commas: ")
    mhcii_lengths_str = input("Enter your MHC class II lengths separated by commas (one length for multiple alleles): ")

    list_of_swissprot_ids = [x.strip() for x in swissprot_ids_str.split(',')]
    list_of_pdb_ids = [x.strip() for x in pdb_ids_str.split(',')]
    mhci_alleles = [x.strip() for x in mhci_alleles_str.split(',')]
    mhci_lengths = [x.strip() for x in mhci_lengths_str.split(',')]
    mhcii_alleles = [x.strip() for x in mhcii_alleles_str.split(',')]
    mhcii_lengths = [x.strip() for x in mhcii_lengths_str.split(',')]

    alignment(list_of_swissprot_ids, matrix='bl62', gapopen=1.53, gapext=0.123, order='aligned', nbtree=2, treeout='true', maxiterate=2, ffts='none')
    path_to_alignment = input("Enter the path to your alignment file: ")
    conserved_sequences_mafft = get_conserved_sequences(path_to_alignment, min_seq_conserved_pos='default', min_seq_flank_pos='default', max_contigous_nonconserved_pos = 8, min_length_block= 10, allowed_gap_pos='None')
    predict_all(conserved_sequences_mafft, mhci_alleles, mhci_lengths, mhcii_alleles, mhcii_lengths, list_of_pdb_ids)
    prediction_results = read_prediction_results()
    analysis_input = make_inputs_for_analysis(prediction_results, list_of_swissprot_ids)
    print("Analysing sequences...")
    analysis_results = analyse_all(analysis_input)
    make_csv_from_results(prediction_results, analysis_results)
    print("Analysis finished. Sequences saved in 'results' folder in csv format")

#run_pipeline()