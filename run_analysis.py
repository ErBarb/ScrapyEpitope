from msa import prediction_choice
from prediction import predict_all
from prediction import get_pdb_from_swissprot
from prediction import analysis_choice
from prediction import swissprotIDSequenceLength
from prediction import epitope_distribution_plots
from prediction import dssp_analysis
from analyse import make_inputs_for_analysis
from analyse import analyse_all
from analyse import make_csv_from_results
from analyse import read_prediction_results



# path_to_mafft_alignment = '/home/erald/Desktop/ScrapyEpitope/msa_results/mafft/mafft.aln-fasta.fasta'
# path_to_muscle_alignment = '/home/erald/Desktop/ScrapyEpitope/msa_results/muscle/muscle.aln-fasta.fasta'
# list_of_swissprot_ids = ['P59594', 'P0DTC2', 'K9N5Q8', 'P36334', 'Q0ZME7', 'P15423', 'Q6Q1S2', 'Q5MQD0', 'Q14EB0']

list_of_swissprot_ids = ['P0DTC2', 'P59594', 'P36334', 'K9N5Q8', 'U3N9S7', 'P15423', 'Q6Q1S2', 
                        'P0DTC9', 'P59595', 'P33469', 'K9N4V7', 'Q5MQC6', 'P15130', 'Q6Q1R8',
                        'P0DTC5', 
                        'P0DTD1', 'P0C6X7', 'P0C6U8', 'P0C6X6', 'K9N7C7', 'P0C6X2', 'P0C6X1', 'P0C6X5', 
                        'P0DTC3', 'P59632']

mhci_alleles = ['HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*02:06', 'HLA-A*03:01', 'HLA-A*11:01', 'HLA-A*23:01', 'HLA-A*24:02', 'HLA-A*25:01', 'HLA-A*26:01', 'HLA-A*29:02', 'HLA-A*30:01',
                'HLA-A*30:02', 'HLA-A*31:01', 'HLA-A*32:01', 'HLA-A*33:03', 'HLA-A*68:01', 'HLA-A*68:02', 'HLA-A*74:01', 'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*13:01', 'HLA-B*13:02',
                'HLA-B*14:02', 'HLA-B*15:01', 'HLA-B*15:02', 'HLA-B*15:25', 'HLA-B*18:01', 'HLA-B*27:02', 'HLA-B*27:05', 'HLA-B*35:01', 'HLA-B*35:03', 'HLA-B*37:01', 'HLA-B*38:01',
                'HLA-B*39:01', 'HLA-B*40:01', 'HLA-B*40:02', 'HLA-B*44:02', 'HLA-B*44:03', 'HLA-B*46:01', 'HLA-B*48:01', 'HLA-B*49:01', 'HLA-B*50:01', 'HLA-B*51:01', 'HLA-B*52:01',
                'HLA-B*53:01', 'HLA-B*55:01', 'HLA-B*56:01', 'HLA-B*57:01', 'HLA-B*58:01', 'HLA-B*58:02', 'HLA-C*01:02', 'HLA-C*02:02', 'HLA-C*02:09', 'HLA-C*03:02', 'HLA-C*03:03',
                'HLA-C*03:04', 'HLA-C*04:01', 'HLA-C*05:01', 'HLA-C*06:02', 'HLA-C*07:01', 'HLA-C*07:02', 'HLA-C*07:04', 'HLA-C*08:01', 'HLA-C*08:02', 'HLA-C*12:02', 'HLA-C*12:03',
                'HLA-C*14:02', 'HLA-C*15:02', 'HLA-C*16:01', 'HLA-C*17:01', 'HLA-E*01:01', 'HLA-E*01:03', 'HLA-G*01:01', 'HLA-G*01:02', 'HLA-G*01:03', 'HLA-G*01:04', 'HLA-G*01:06']
mhci_lengths = [8, 9, 10, 11, 12, 13, 14]

mhcii_alleles = ['HLA-DRB1*01:01', 'HLA-DRB1*07:01', 'HLA-DRB1*03:01', 'HLA-DRB1*04:01', 'HLA-DRB1*04:05', 'HLA-DRB1*08:02', 'HLA-DRB1*09:01', 'HLA-DRB1*11:01',
                'HLA-DRB1*12:01', 'HLA-DRB1*13:02', 'HLA-DRB1*15:01', 'HLA-DRB3*01:01', 'HLA-DRB3*02:02', 'HLA-DRB4*01:01', 'HLA-DRB5*01:01', 'HLA-DQA1*05:01/DQB1*02:01',
                'HLA-DQA1*05:01/DQB1*03:01', 'HLA-DQA1*03:01/DQB1*03:02', 'HLA-DQA1*04:01/DQB1*04:02', 'HLA-DQA1*01:01/DQB1*05:01', 'HLA-DQA1*01:02/DQB1*06:02', 'HLA-DPA1*02:01/DPB1*01:01',
                'HLA-DPA1*01:03/DPB1*02:01', 'HLA-DPA1*01:03/DPB1*04:01', 'HLA-DPA1*03:01/DPB1*04:02', 'HLA-DPA1*02:01/DPB1*05:01', 'HLA-DPA1*02:01/DPB1*14:01']
mhcii_lengths = [12, 13, 14, 15, 16, 17, 18]

list_of_pdb_ids = ['6VSB', '5WRG', '7M51', '4L72', '7CYD', '7KIP', 
                    '6WJI', '2CJR', '4J3K', '6KL2', '7LGT', '5EPW', 
                    '5R7Y', '1Q2W', '7NH7', '4WUR', '4RS4', '5NH0', 
                    '6XDC']


# sequences = prediction_choice(list_of_swissprot_ids)
# print("Protein sequences collected")
# print("Trying to get PDB IDs...")
# list_of_pdb_ids = get_pdb_from_swissprot(list_of_swissprot_ids)
# print("PDB IDs collected")
# predict_all(mhci_alleles, mhci_lengths, mhcii_alleles, mhcii_lengths, list_of_pdb_ids)
# print("Epitope prediction done")
# analysis_choice(list_of_swissprot_ids, list_of_pdb_ids)

prediction_results = read_prediction_results()
analysis_input = make_inputs_for_analysis(prediction_results, list_of_swissprot_ids)
print("Analysing sequences...")
analysis_results = analyse_all(analysis_input)
#make_csv_from_results(prediction_results, analysis_results)
print("Analysis finished. Sequences saved in 'results' folder in csv format")




def run_pipeline():

    """This function will run the pipeline. It first requires the user to input the swissprot and pdb ids, then the mhc class 
    I and II alleles and lengths. Then it runs the alignment algorithms on the swissprot ids, asks the user for the path to the 
    alingment file, thus the user can choose either the muscle or the mafft alignment, gets the conserved sequences, predicts the 
    epitopes and then analyses all the epitopes and saves the results in csv files in folder 'results'."""

    swissprot_ids_str = input("Enter your Swissprot IDs separated by commas: ")
    mhci_alleles_str = input("Enter your MHC class I alleles separated by commas: ")
    mhci_lengths_str = input("Enter your MHC class I lengths separated by commas (one length for each allele): ")
    mhcii_alleles_str = input("Enter your MHC class II alleles separated by commas: ")
    mhcii_lengths_str = input("Enter your MHC class II lengths separated by commas (one length for multiple alleles): ")

    list_of_swissprot_ids = [x.strip() for x in swissprot_ids_str.split(',')]
    mhci_alleles = [x.strip() for x in mhci_alleles_str.split(',')]
    mhci_lengths = [x.strip() for x in mhci_lengths_str.split(',')]
    mhcii_alleles = [x.strip() for x in mhcii_alleles_str.split(',')]
    mhcii_lengths = [x.strip() for x in mhcii_lengths_str.split(',')]

    sequences = prediction_choice()
    list_of_pdb_ids = get_pdb_from_swissprot(list_of_swissprot_ids)
    predict_all(sequences, mhci_alleles, mhci_lengths, mhcii_alleles, mhcii_lengths, list_of_pdb_ids)
    analysis_choice(list_of_swissprot_ids, list_of_pdb_ids)
    
    prediction_results = read_prediction_results()
    analysis_input = make_inputs_for_analysis(prediction_results, list_of_swissprot_ids)
    print("Analysing sequences...")
    analysis_results = analyse_all(analysis_input)
    make_csv_from_results(prediction_results, analysis_results)
    print("Analysis finished. Sequences saved in 'results' folder in csv format")

#run_pipeline()