import os

def analyse(seq_list):

    # """This function uses REST API services from https://www.ebi.ac.uk/ to output alignment files using the multiple
    # sequence alignment methods "MAFFT" and "MUSCLE". The argument is a list of swissprot protein IDs. The files are
    # saved in their respective folders in msa_results"""

    for seq in seq_list:
        os.system(
            'python embosspepstats.py --email erald.bb@gmail.com --sequence ' + seq + ' --outfile pepstats_results')
        with open('pepstats_results.out.txt') as f:
             lines = f.readlines()
             aa_properties = lines[-11:-1]
             exp_inclusion_bodies = float(lines[7].split()[-1])
             print(aa_properties)
             print(exp_inclusion_bodies)
        #os.remove(os.getcwd()+"/pepstats_results.out.txt")
        #os.remove(os.getcwd()+"/pepstats_results.sequence.txt")

sequence_list = ['WMYTGSGYYYPEPITENNVVV','FKEELDQWFKNQTSVAPDL','SRSARSAIEDLLFDKVTIADPGYMQGYDDC','WSYTGSSFYAPEPITSLNTKY']
analyse(sequence_list)