import os

def pepstats(seq_list):

    # """This function uses REST API services from https://www.ebi.ac.uk/ to output alignment files using the multiple
    # sequence alignment methods "MAFFT" and "MUSCLE". The argument is a list of swissprot protein IDs. The files are
    # saved in their respective folders in msa_results"""

    for seq in seq_list:
        os.system(
            'python embosspepstats.py --email erald.bb@gmail.com --sequence ' + seq + ' --outfile pepstats_results --quiet')
        results_row = []
        with open('pepstats_results.out.txt') as f:
            lines = f.readlines()
            aa_properties = lines[-11:-1]
            exp_inclusion_bodies = float(lines[7].split()[-1])

            tiny_aa = float(aa_properties[1].split('\t')[-1][:-1])
            small_aa = float(aa_properties[2].split('\t')[-1][:-1])
            aliphatic_aa = float(aa_properties[3].split('\t')[-1][:-1])
            aromatic_aa = float(aa_properties[4].split('\t')[-1][:-1])
            non_polar_aa = float(aa_properties[5].split('\t')[-1][:-1])
            polar_aa = float(aa_properties[6].split('\t')[-1][:-1])
            charged_aa = float(aa_properties[7].split('\t')[-1][:-1])
            basic_aa = float(aa_properties[8].split('\t')[-1][:-1])
            acidic_aa = float(aa_properties[9].split('\t')[-1][:-1])
             
            results_row.append(tiny_aa)
            results_row.append(small_aa)
            results_row.append(aliphatic_aa)
            results_row.append(aromatic_aa)
            results_row.append(non_polar_aa)
            results_row.append(polar_aa)
            results_row.append(charged_aa)
            results_row.append(basic_aa)
            results_row.append(acidic_aa)
            results_row.append(exp_inclusion_bodies)

            print(aa_properties)
            print(exp_inclusion_bodies)
        print(results_row)

        os.remove(os.getcwd()+"/pepstats_results.out.txt")
        os.remove(os.getcwd()+"/pepstats_results.sequence.txt")

sequence_list = ['WMYTGSGYYYPEPITENNVVV','FKEELDQWFKNQTSVAPDL','SRSARSAIEDLLFDKVTIADPGYMQGYDDC','WSYTGSSFYAPEPITSLNTKY']
pepstats(sequence_list)