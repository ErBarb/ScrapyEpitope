from prediction import predict_all
from msa import alignment
from msa import get_conserved_sequences

example_mhci_results = [['protein_id', 'conserved_sequence', 'allele', 'seq_num', 'start', 'end', 'length', 'peptide', 'core', 'icore', 'score', 'percentile_rank'],
['P0DTC2', 'SRSARSAIEDLLFDKVTIADPGYMQGYDDC', 'HLA-A*02:01', '1', '11', '19', '9', 'LLFDKVTIA', 'LLFDKVTIA', 'LLFDKVTIA', '0.88', 0.04], 
['P0DTC2', 'SRSARSAIEDLLFDKVTIADPGYMQGYDDC', 'HLA-A*03:01', '1', '16', '25', '10', 'VTIADPGYMQ', 'VTIDPGYMQ', 'VTIADPGYMQ', '0.0015', 9.9],
['P0DTC2', 'WSYTGSSFYAPEPITSLNTKY', 'HLA-A*02:01', '1', '9', '17', '9', 'YAPEPITSL', 'YAPEPITSL', 'YAPEPITSL', '0.586', 0.21],
['P36334', 'WMYTGSGYYYPEPITENNVVV', 'HLA-A*01:01', '1', '3', '10', '8', 'YTGSGYYY', 'YTG-SGYYY', 'YTGSGYYY', '0.559', 0.16]]
example_mhci_proc_results = [['protein_id', 'conserved_sequence', 'allele', 'start', 'end', 'peptide_length', 'peptide', 'proteasome_score', 'tap_score', 'mhci_score', 'processing_score', 'total_score', 'mhci_ic50'],
['P36334', 'WMYTGSGYYYPEPITENNVVV', 'HLA-A*03:01', '1', '1', '10', '10', 'WMYTGSGYYY', '1.1608', '1.3381', '-2.0228', '2.4988', 0.476, '105.4']]
example_mhcii_results = [['protein_id', 'conserved_sequence', 'allele', 'seq_num', 'start', 'end', 'length', 'core_peptide', 'peptide', 'ic50', 'rank', 'adjusted_rank'],
['P0DTC2', 'SRSARSAIEDLLFDKVTIADPGYMQGYDDC', 'HLA-DRB1*03:01', '1', '7', '20', '14', 'Consensus (smm/nn/sturniolo)', 'AIEDLLFDKVTIAD', '7.60', 8.18, '-', '-', '-', '-', 'LLFDKVTIA', '562.00', '7.60', '8.18', 'LLFDKVTIA', '221.60', '8.20', '8.83', '-', '-', '-', '-', 'LLFDKVTIA', '5.20', '1.30', '1.40'],
['P0DTC2', 'SRSARSAIEDLLFDKVTIADPGYMQGYDDC', 'HLA-DRB1*03:01', '1', '9', '21', '13', 'Consensus (smm/nn/sturniolo)', 'EDLLFDKVTIADP', '6.40', 9.98, '-', '-', '-', '-', 'LLFDKVTIA', '610.00', '6.40', '9.98', 'LLFDKVTIA', '376.30', '9.70', '15.13', '-', '-', '-', '-', 'LLFDKVTIA', '5.20', '1.10', '1.72'],
['P0DTC2', 'SRSARSAIEDLLFDKVTIADPGYMQGYDDC', 'HLA-DRB1*03:01', '1', '10', '25', '16', 'Consensus (smm/nn/sturniolo)', 'DLLFDKVTIADPGYMQ', '8.30', 9.58, '-', '-', '-', '-', 'LLFDKVTIA', '1625.00', '13.00', '15.01', 'LLFDKVTIA', '273.60', '8.30', '9.58', '-', '-', '-', '-', 'LLFDKVTIA', '5.20', '1.60', '1.85']]
example_bepipred2_results = [['protein_id', 'conserved_sequence', 'predicted_epitope', 'start_position', 'end_position'],
['P0DTC2', 'WSYTGSSFYAPEPITSLNTKY', 'GSSFYAPEPITSLN', 4, 17],
['P36334', 'WMYTGSGYYYPEPITENNVVV', 'YYYPEPITEN', 7, 16],
['P36334', 'FKEELDQWFKNQTSVAPDL', 'LDQWFKNQTSVA', 4, 15]]
example_bepipred_results = [['protein_id', 'conserved_sequence', 'predicted_epitope', 'start_position', 'end_position'],
['P0DTC2', 'SRSARSAIEDLLFDKVTIADPGYMQGYDDC', 'IADPGYMQGYDDC', 18, 30],
['P0DTC2', 'WSYTGSSFYAPEPITSLNTKY', 'SSFYAPEPITSL', 6, 17],
['P36334', 'WMYTGSGYYYPEPITENNVVV', 'SGYYYPEPITE', 6, 16],
['P36334', 'FKEELDQWFKNQTSVAPDL', 'QTSVAPDL', 12, 19]]
example_emini_results = [['protein_id', 'conserved_sequence', 'predicted_epitope', 'start_position', 'end_position'],
['P0DTC2', 'SRSARSAIEDLLFDKVTIADPGYMQGYDDC', 'PGYMQGY', 21, 27],
['P36334', 'WMYTGSGYYYPEPITENNVVV', 'YYYPEPIT', 8, 15]]
example_choufasman_results = [['protein_id', 'conserved_sequence', 'predicted_epitope', 'start_position', 'end_position'],
['P0DTC2', 'SRSARSAIEDLLFDKVTIADPGYMQGYDDC', 'ADPGYMQGY', 19, 27],
['P36334', 'WMYTGSGYYYPEPITENNVVV', 'TGSGYYY', 4, 10],
['P36334', 'FKEELDQWFKNQTSVAPDL', 'WFKNQTSVA', 8, 16]]
example_karplusschulz_results = [['protein_id', 'conserved_sequence', 'predicted_epitope', 'start_position', 'end_position']]
example_kolaskar_results = [['protein_id', 'conserved_sequence', 'predicted_epitope', 'start_position', 'end_position'],
['P0DTC2', 'SRSARSAIEDLLFDKVTIADPGYMQGYDDC', 'FDKVTIA', 13, 19],
['P36334', 'WMYTGSGYYYPEPITENNVVV', 'SGYYYPE', 6, 12]]
example_parker_results = [['protein_id', 'conserved_sequence', 'predicted_epitope', 'start_position', 'end_position']]
example_ellipro_results = ([['pdb_id', 'chain', 'start', 'end', 'peptide', 'nr_of_residues', 'score'], 
['6vxx', 'A', '1071', '1147', 'QEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDS', '77', 0.88], 
['6vxx', 'A', '92', '192', 'FASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVNCTFEYVSFKNLREF', '67', 0.788], 
['6vxx', 'A', '328', '364', 'RFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVAD', '37', 0.767], 
['6vxx', 'A', '433', '537', 'VIAWNSNNLDSKGNYNYLYRKPFERDIYFPLQSYGFQPTNVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNK', '75', 0.765], 
['6vxx', 'A', '236', '267', 'TRFQTLLALHAAYYV', '15', 0.757], ['6vxx', 'A', '62', '86', 'VTWFHAIHDNPVLPF', '15', 0.735], 
['6vxx', 'A', '391', '430', 'CFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFT', '40', 0.72], ['6vxx', 'A', '553', '564', 'TESNKKFLPFQQ', '12', 0.721], 
['6vxx', 'A', '205', '219', 'SKHTPINLVRDLPQG', '15', 0.703], ['6crz', 'A', '1049', '1120', 'YVPSQERNFTTAPAICHEGKAYFPREGVFVFNGTSWFITQRNFFSPQIITTDNTFVSGNCDVVIGIINNTVY', '72', 0.879], 
['6crz', 'A', '420', '511', 'VLAWNTRNIDATSTGNYNYKYRYLRHGKLRPFERDISNVPFSPDGKPCTPPALNCYWPLNDYGFYTTTGIGYQPYRVVVLSFETVC', '86', 0.853], 
['6crz', 'A', '232', '252', 'RAILTAFSPIWGTSAAAY', '18', 0.802], ['6crz', 'A', '18', '31', 'RCTTFDDVQAPNYT', '14', 0.801], 
['6crz', 'A', '105', '183', 'STMNNKSQSVIIINNSTNVVIRACNFELCDNPFFAVSKPMGTQTHTMIFDNAFNCTFEYISDAFSLDVSEKSGNFKHLR', '79', 0.758], 
['6crz', 'A', '68', '80', 'GFHTINHTFGNPV', '13', 0.733], ['6crz', 'A', '684', '701', 'DSSIAYSNNTIAIPTNFS', '18', 0.713], 
['5x5f', 'A', '1141', '1206', 'YYPSNHIEVVSAYGLCDAANPTNCIAPVNGYFIKTNNTRIVDEWSYTGSSFYAPEPITSLNTKYVA', '66', 0.873], 
['5x5f', 'A', '481', '572', 'LATVPHNLTTITKPLKYSYINKCSRLLSDDRTEVPQLVNANQYSPCVSIVPSTVWEDGDYYRKQLSPLEGGGWLVASGSTVAMTEQLQMGFG', '92', 0.869], 
['5x5f', 'A', '86', '104', 'VYSAGHATGTTPQKLFVAN', '19', 0.804], ['5x5f', 'A', '296', '312', 'IPHSIRSIQSDRKAWAA', '17', 0.752], 
['5x5f', 'A', '18', '45', 'YVDVGPDSVKSACIEVDIQQTFFDKTWP', '28', 0.735], ['5x5f', 'A', '765', '787', 'NHPIQVDQLNSSYFKLSIPTNFS', '23', 0.733], 
['5x5f', 'A', '122', '144', 'AAANSTGTVIISPSTSATIRKIY', '23', 0.719], ['5x5f', 'A', '376', '410', 'EQVECDFSPLLSGTPPQVYNFKRLVFTNCNYN', '32', 0.716], 
['5x5f', 'A', '172', '250', 'LPDGCGTLLRAFYCILEPRSGNHCPAGNSYTSFATYHTPATDCSDGNYNRNASLNSFKEYFNLRNCTFMYTYNITEDEI', '79', 0.703]],
[['pdb_id', 'chain', 'start', 'end', 'peptide', 'nr_of_residues', 'score'], 
['6vxx', 'A', '27', '267', 'AYTVTWFHAIHDNPVLPNFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVNCTFEYVSFKNLREFGSKHTPINLVRDLPSLLVLPIGIITRFQTLLALHAAYYV', '124', '0.741'], 
['6vxx', 'A', '707', '1147', 'YSNNSIAIPTNFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRPLLEMQYSAAGITSGWTFGAGAALQIPFAMQMAYFNGIGVTQNVLYENQKLIANLGQKQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDS', '180', '0.732'], 
['6vxx', 'A', '328', '586', 'RFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADVLSASFSTYNLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTVIAWNSNNLDSKGNYNYLYRKPFERDIFPLQSYGFQPTNVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKNTESNKFLPFQQVRDPQTLEILD', '185', '0.729'], 
['6crz', 'A', '18', '252', 'RCTTFDDVQAPNTGFHTINHTFGNPVATEKSNVVRGWVSTMNNKSQSVIIINNSTNVVIRACNFELCDNPFFAVSKPMGTQTHTMIFDNAFNCTFEYISDAFSLDVSEKSGNFKHLRQPIDVVRDLPLIRAILTAFSPIWGTSAAAY', '147', '0.746'], 
['6crz', 'A', '315', '571', 'RFPNITNLCPFGEVFNATKFPSVYAWERKKISNCVADYVLNSTFFSTFYLCFSNVYADSFVVKGDDVRQIAPGQTGVIADYNYKLPDDFMGVLAWNTRNIDATSTGNYNYKYRYLRHGKLRPERDISNVPFSPDGKPCTPPALNCYWPLNDYGFYTTTGIGYQPYRVVVLSFETVCGKLSTDLIKNNPSSKFQPFQQRDPKTSEIL', '206', '0.734'], 
['6crz', 'A', '684', '1120', 'DSSIAYSNNTIAIPTNFSAQVKQMYKTPTLKYFGGFNFSQILPDPLKPTKPLLDMAYAAVSGTATAGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKQIANQNCVLGQSGYLYPSQERNFTTAPAICHEGKAYFPGVFVFNGTSWFITQRNFFSPQIITTDNTFVSGNCDVVIGIINNTVY', '187', '0.729'], 
['5x5f', 'A', '376', '649', 'EQVECDFSPLLSGTPPQVYNFKRLVFTNCNNKLSLFSSSLILDYFSYPLSMKSDSSSAGPISQFNYKQLATVPHNLTTITKPLKYSYINKCSRLLSDDRTEVPQLVNANQYSPCVSIVPSTVWEDGDYYKQLSPLEGGGWLVASGSTVAMTEQLQMGFGTDTNSKLLYVFQNCTAVGVQQRFVGYYSDDGNYY', '193', '0.75'], 
['5x5f', 'A', '765', '1206', 'NHPIQVDQLNSSYFKLSIPTNFSFDRNFASVKSSQSSPIIPGFGGDFNLTLLEPVARPLMDVNMAAYTSSLLGSIAGVGWTAGLSSFAAIPFAQSIFYRLNGVGITQQVLSENQKLIANKFNNKAQSKRGTIVYYPSNHIEVVSAYGLCDAANPTNCIAPVNGYFIKTNNTRIVDEWSYTGSSFYAPEPITSLNTKYVA', '199', '0.737'], 
['5x5f', 'A', '18', '312', 'YVDVGPDSVKSACIEVDIQQTFFDKTWVYSAGHATGTTPQKLFVANQGAAANSTGTVIISPSTSATIRKIYNFSDGKRFNLPDGCGTLRFCILEPRSGNHCPAGNSYTSFATYHTPATDCSDGNYNRNASLNSFKEYFNLRNCTMTNITEDEIIPHSIRSIQSDRKAWAA', '170', '0.727']])
tuple_of_lists = (example_mhci_results,example_mhci_proc_results,example_mhcii_results,example_bepipred2_results,example_bepipred_results,example_emini_results,example_choufasman_results,example_karplusschulz_results,example_kolaskar_results,example_parker_results,example_ellipro_results)


def analyse_antigenicity_vaxijen2():
    pass
def analyse_toxicity_and_other_toxinpred():
    pass
def analyse_allergenicity_algpred2():
    pass
def analyse_protparam():
    pass
def analyse_immunogenicity_both_mhc_classes():
    pass
def analyse_population_coverage():
    pass
def cluster_analysis():
    pass
def conservancy_analysis():
    pass
def pymol_visualization():
    pass