from msa import alignment
from msa import get_conserved_sequences
from prediction import predict_all

def run_pipeline():
    swissprot_ids_str = input("Enter your Swissprot IDs separated by commas: ")
    pdb_ids_str = input("Enter your PDB IDs separated by commas: ")
    list_of_pdb_ids = [x.strip() for x in pdb_ids_str.split(',')]
    path_to_alignment = input("Enter the path to your alignment file: ")
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
    conserved_sequences_mafft = get_conserved_sequences(path_to_alignment, min_seq_conserved_pos='default', min_seq_flank_pos='default', max_contigous_nonconserved_pos = 8, min_length_block= 10, allowed_gap_pos='None')
    predicted_epitopes = predict_all(conserved_sequences_mafft, mhci_alleles, mhci_lengths, mhcii_alleles, mhcii_lengths, list_of_pdb_ids)