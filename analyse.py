import re
import os
import requests
import time
import csv
from csv import reader
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from urllib import request
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as ec
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.keys import Keys



def read_prediction_results():

    """This function will read all the results from the csv files created by the prediction methods in prediction.py, then return them
        in a tuple of lists"""

    with open('epitope_prediction_results/mhci_epitopes.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        mhci_prediction_list = list(csv_reader)
    with open('epitope_prediction_results/mhci_proc_epitopes.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        mhci_proc_prediction_list = list(csv_reader)
    with open('epitope_prediction_results/mhcii_epitopes.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        mhcii_prediction_list = list(csv_reader)
    with open('epitope_prediction_results/bepipred2.0_epitopes.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        bepipred2_prediction_list = list(csv_reader)
    with open('epitope_prediction_results/bepipred1.0_epitopes.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        bepipred_prediction_list = list(csv_reader)
    with open('epitope_prediction_results/emini_epitopes.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        emini_prediction_list = list(csv_reader)
    with open('epitope_prediction_results/choufasman_epitopes.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        choufasman_prediction_list = list(csv_reader)
    with open('epitope_prediction_results/karplusschulz_epitopes.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        karplusschulz_prediction_list = list(csv_reader)
    with open('epitope_prediction_results/kolaskartongaonkar_epitopes.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        kolaskartongaonkar_prediction_list = list(csv_reader)
    with open('epitope_prediction_results/parker_epitopes.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        parker_prediction_list = list(csv_reader)
    with open('epitope_prediction_results/ellipro_linear_epitopes.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        ellipro_linear_prediction_list = list(csv_reader)
    with open('epitope_prediction_results/ellipro_discontinous_epitopes.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        ellipro_discontinous_prediction_list = list(csv_reader)
    with open('epitope_prediction_results/discotope_epitopes.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        discotope_prediction_list = list(csv_reader)

    return mhci_prediction_list, mhci_proc_prediction_list, mhcii_prediction_list, bepipred2_prediction_list, \
        bepipred_prediction_list, emini_prediction_list, choufasman_prediction_list, karplusschulz_prediction_list, \
            kolaskartongaonkar_prediction_list, parker_prediction_list, ellipro_linear_prediction_list, ellipro_discontinous_prediction_list, \
                discotope_prediction_list



def make_inputs_for_analysis(results_from_prediction, list_of_swissprot_ids):

    """This function will get the prediction results and swissprot ids and make the inputs needed for all of the analysis tools and
    return them"""

    pop_cov_input = ''
    for i in range(len(results_from_prediction)):
        if i == 0 or i == 1:
            for x in results_from_prediction[i]:
                peptide = x[7]
                allele = x[2]
                if peptide != 'peptide' or allele != 'allele':
                    pop_cov_input = pop_cov_input + peptide + '\t' + allele + '\n'
                x.remove(x[7])
                x.insert(0, peptide)
        elif i == 2:
           for x in results_from_prediction[i]:
                peptide = x[8]
                allele = x[2]
                if peptide != 'peptide' or allele != 'allele':
                    pop_cov_input = pop_cov_input + peptide + '\t' + allele + '\n'
                x.remove(x[8])
                x.insert(0, peptide)
        elif i == 3 or i == 4 or i == 5 or i == 6 or i == 7 or i == 8 or i == 9:
            for x in results_from_prediction[i]:
                peptide = x[2]
                x.remove(x[2])
                x.insert(0, peptide)
        elif i == 10:
            for x in results_from_prediction[i]:
                peptide = x[4]
                x.remove(x[4])
                x.insert(0, peptide)
        elif i == 11 or i == 12:
            for x in results_from_prediction[i]:
                peptide = x[2]
                x.remove(x[2])
                x.insert(0, peptide)


    list_of_all_linear_epitopes = []
    immunogenicity_input = ''
    immunogenicity_indexes = len(results_from_prediction[0][1:]) + len(results_from_prediction[1][1:])
    for i in range(len(results_from_prediction)-2):
        for x in results_from_prediction[i][1:]:
            list_of_all_linear_epitopes.append(x[0])
        if i == 0 or i == 1:
            for x in results_from_prediction[i][1:]:
                immunogenicity_input = immunogenicity_input + x[0] + '\n'

    input_string = ''
    for i in range(len(list_of_all_linear_epitopes)):
        input_string = input_string + '>seq' + str(i+1) + '\n' + list_of_all_linear_epitopes[i] + '\n'

    toxinpred_chunks = []
    toxinpred_epitopes = []
    toxinpred_excluded_indexes = []

    for i in range(len(list_of_all_linear_epitopes)):
        if len(list_of_all_linear_epitopes[i]) <= 50:
            toxinpred_epitopes.append([list_of_all_linear_epitopes[i]])
        else:
            toxinpred_excluded_indexes.append(i)

    for i in range(0, len(toxinpred_epitopes), 400):
        toxinpred_400e_chunk = toxinpred_epitopes[i:i + 400]
        toxinpred_input = ''
        for n in range(len(toxinpred_400e_chunk)):
            toxinpred_input = toxinpred_input + '>seq' + str(n + 1) + '\n' + toxinpred_400e_chunk[n][0] + '\n'
        if toxinpred_input != '':
            toxinpred_chunks.append(toxinpred_input)

    algpred_chunks = []
    for i in range(0, len(list_of_all_linear_epitopes), 400):
        algpred_400e_chunk = list_of_all_linear_epitopes[i:i + 400]
        algpred_input = ''
        for n in range(len(algpred_400e_chunk)):
            algpred_input = algpred_input + '>seq' + str(n + 1) + '\n' + algpred_400e_chunk[n][0] + '\n'
        if algpred_input != '':
            algpred_chunks.append(algpred_input)

    all_protein_fasta = ''
    for id in list_of_swissprot_ids:
        baseUrl="http://www.uniprot.org/uniprot/"
        currentUrl=baseUrl+id+".fasta"
        response = requests.post(currentUrl)
        cData=''.join(response.text)
        all_protein_fasta = all_protein_fasta + cData

    seq_file = open('seq_file.txt', 'a')
    seq_file.write(input_string[:-1])
    seq_file.close()

    pop_cov_file = open('pop_cov_file.txt', 'a')
    pop_cov_file.write(pop_cov_input[:-1])
    pop_cov_file.close()

    immunogenicity_file = open('immunogenicity_file.txt', 'a')
    immunogenicity_file.write(immunogenicity_input)
    immunogenicity_file.close()

    cluster_file = open('cluster_file.txt', 'a')
    cluster_file.write(input_string)
    cluster_file.close()

    conservancy_file = open('conservancy_seq_file.txt', 'a')
    conservancy_file.write(input_string)
    conservancy_file.close()
    conservancy_file = open('conservancy_protein_file.txt', 'a')
    conservancy_file.write(all_protein_fasta)
    conservancy_file.close()

    current_directory = os.getcwd()
    final_directory = os.path.join(current_directory, r'results')
    if not os.path.exists(final_directory):
        os.makedirs(final_directory)

    return list_of_all_linear_epitopes, toxinpred_chunks, toxinpred_epitopes, toxinpred_excluded_indexes, immunogenicity_indexes, algpred_chunks

def analyse_all(tuple_inputs):

    """This function uses Selenium to access the following tools: Toxinpred, Algpred2, Vaxijen and IEDB tools such Immunogenicity for
    MHCI, Population Coverage, Cluster and Conservancy analysis. It also uses Protparam module from Biopython to analyse each sequence.
    Population coverage, cluster and conservancy results are stored in files. The other ones are returned as list of lists."""

    options = Options()
    options.headless = True

    list_of_linear_epitopes = tuple_inputs[0]
    toxinpred_chunks = tuple_inputs[1]
    toxinpred_epitopes = tuple_inputs[2]
    toxinpred_excluded_indexes = tuple_inputs[3]
    immunogenicity_indexes = tuple_inputs[4]
    algpred_chunks = tuple_inputs[5]

    toxinpred_url = 'https://webs.iiitd.edu.in/raghava/toxinpred/multi_submit.php'
    algpred2_url = 'https://webs.iiitd.edu.in/raghava/algpred2/batch.html'
    vaxijen_url = 'http://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html'
    immunogenicity_mhci_url = 'http://tools.iedb.org/immunogenicity/'
    # immunogenicity_mhcii_url = ''
    population_coverage_url = 'http://tools.iedb.org/population/'
    cluster_analysis_url = 'http://tools.iedb.org/cluster/'
    conservancy_analysis_url = 'http://tools.iedb.org/conservancy/'

    analysis_results = []
    for index, value in enumerate(list_of_linear_epitopes):
        row = []
        row.insert(0,list_of_linear_epitopes[index])
        analysis_results.append(row)

    # Protparam
    for index, epitope in enumerate(list_of_linear_epitopes):
        row = []
        x = ProteinAnalysis(epitope)
        #aa_percentage = x.get_amino_acids_percent()
        #aa_count = x.count_amino_acids()
        mol_weight = x.molecular_weight()
        mol_weight = round(mol_weight, 3)
        isoel_point = x.isoelectric_point()
        isoel_point = round(isoel_point, 3)
        aromaticity = x.aromaticity()
        aromaticity = round(aromaticity, 3)
        insta_index = x.instability_index()
        insta_index = round(insta_index, 3)
        sec_struc = x.secondary_structure_fraction()
        sec_struc = [round(num, 3) for num in sec_struc]
        helix_2_struc = sec_struc[0]
        turn_2_struc = sec_struc[1]
        sheet_2_struc = sec_struc[2]
        epsilon_prot = x.molar_extinction_coefficient()
        reduCys = epsilon_prot[0]
        disulfBridge = epsilon_prot[1]
        hydropathicity = x.gravy()
        hydropathicity = round(hydropathicity, 3)
        #flex_list = [str(round(num, 2)) for num in x.flexibility()]
        #flexibility = ': '.join(flex_list)
        chpH = x.charge_at_pH(7)
        chpH = round(chpH, 3)
        #row.extend((mol_weight, isoel_point, aromaticity, insta_index, helix_2_struc, turn_2_struc, sheet_2_struc, reduCys, disulfBridge, hydropathicity, flexibility, chpH))
        row.extend((mol_weight, isoel_point, aromaticity, insta_index, helix_2_struc, turn_2_struc, sheet_2_struc, reduCys, disulfBridge, hydropathicity, chpH))
        for parameter in row:
            analysis_results[index].append(parameter)


    #options=options,

    immunogenicity_mhci = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
    immunogenicity_mhci.get(immunogenicity_mhci_url)
    immunogenicity_mhci.find_element(By.NAME, "sequence_file").send_keys(os.getcwd()+"/immunogenicity_file.txt")
    immunogenicity_mhci.find_element(By.NAME, "submit").click()

    vaxijen = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
    vaxijen.get(vaxijen_url)
    vaxijen.find_element(By.XPATH, "/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[1]/td[2]/p/input").send_keys(os.getcwd()+"/seq_file.txt")
    vaxijen.find_element(By.XPATH, "/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[2]/td[2]/p/select/option[2]").click()
    vaxijen.find_element(By.XPATH, "/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[2]/td[3]/input").send_keys('0.5')
    vaxijen.find_element(By.XPATH, "/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[3]/td[2]/input[1]").click()

    if len(list_of_linear_epitopes) < 3000:
        cluster = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
        cluster.maximize_window()
        cluster.get(cluster_analysis_url)
        cluster.find_element(By.NAME, "sequence_file").send_keys(os.getcwd()+"/cluster_file.txt")
        cluster.find_element(By.NAME, "submit").click()
        wait_cluster = WebDriverWait(cluster, 6000)
        wait_cluster.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div[3]/form/table/tbody/tr[2]/th")))
        cluster.find_element(By.NAME, "submit").click()
        wait_cluster.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div[3]/form/table/tbody")))
        cluster.find_element(By.NAME, "submit").click()

    conservancy = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
    conservancy.maximize_window()
    conservancy.get(conservancy_analysis_url)
    conservancy.find_element(By.NAME, "epitope_file").send_keys(os.getcwd()+"/conservancy_seq_file.txt")
    conservancy.find_element(By.NAME, "protein_file").send_keys(os.getcwd()+"/conservancy_protein_file.txt")
    conservancy.find_element(By.NAME, "submit").click()

    population_coverage = webdriver.Firefox(executable_path = '../ScrapyEpitope/geckodriver')
    population_coverage.maximize_window()
    population_coverage.get(population_coverage_url)
    population_coverage.find_element(By.ID, "id_epitope_allele_file").send_keys(os.getcwd()+"/pop_cov_file.txt")
    population_coverage.find_element(By.XPATH, "/html/body/div[3]/form/table[2]/tbody/tr[2]/td[1]/select/option[1]").click()
    population_coverage.find_element(By.XPATH, "/html/body/div[3]/form/table[2]/tbody/tr[4]/td/input").click()
    time.sleep(3)
    population_coverage.find_element(By.NAME, "submit").click()


    wait_vaxijen = WebDriverWait(vaxijen, 6000)
    wait_vaxijen.until(ec.visibility_of_element_located((By.CLASS_NAME, "boilerplate")))
    vaxijen_results_body = vaxijen.find_element(By.XPATH, '/html/body/div/table/tbody/tr[4]/td[3]/table/tbody').text.splitlines()
    vaxijen.close()
    print("Vaxijen done")
    vaxijen_results = []
    for value in vaxijen_results_body:
        if value.startswith('Overall') == True:
            result = []
            antigenicity_val = float(re.search(r'[-+]?\d*\.*\d+', value).group())
            result.append(antigenicity_val)
            if antigenicity_val >= 0.5:
                result.append('Probable Antigen')
            else:
                result.append('Probable Non-Antigen')
            vaxijen_results.append(result)
    for i in range(len(vaxijen_results)):
       analysis_results[i].append(vaxijen_results[i][0])
       analysis_results[i].append(vaxijen_results[i][1])



    wait_conservancy = WebDriverWait(conservancy, 6000)
    wait_conservancy.until(ec.visibility_of_element_located((By.ID, "result_table")))
    conservancy_results_columns = conservancy.find_element(By.XPATH, '/html/body/div[3]/table/thead/tr').text.splitlines()
    conservancy_results_table = conservancy.find_element(By.XPATH, '/html/body/div[3]/table/tbody').text.splitlines()
    conservancy.close()
    print("Conservancy done")
    conservancy_list_of_rows = []
    conservancy_list_of_rows.append(conservancy_results_columns)
    for row in conservancy_results_table:
        conservancy_row = row.split(' ')
        conservancy_list_of_rows.append(conservancy_row)
    with open('results/conservancy_analysis.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerows(conservancy_list_of_rows)



    if len(list_of_linear_epitopes) < 3000:
        wait_cluster.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div[3]/div/div[3]/div[2]")))
        cluster_results_table = cluster.find_element(By.XPATH, '/html/body/div[3]/div/div[3]/div[2]/div[1]/table/tbody').text.splitlines()
        cluster.close()
        print("Cluster done")
        cluster_list_of_rows = []
        cluster_columns = ['Cluster.Sub-Cluster Number','Peptide Number','Alignment','Position','Description','Peptide']
        cluster_list_of_rows.append(cluster_columns)
        for i in range(len(cluster_results_table)):
            if i == 0:
                continue
            else:
                cluster_row = cluster_results_table[i].split(' ')
                if len(cluster_row) == 6:
                    cluster_list_of_rows.append(cluster_row)
                else:
                    indexes_to_remove = len(cluster_row) - 6
                    seq = cluster_row.pop(4)
                    for i in range(indexes_to_remove):
                        removed_duplicate = cluster_row.pop(4)
                        seq = seq + removed_duplicate
                    cluster_row.insert(4, seq)
                    cluster_list_of_rows.append(cluster_row)
        with open('results/cluster_analysis.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerows(cluster_list_of_rows)



    wait_population_coverage = WebDriverWait(population_coverage, 6000)
    wait_population_coverage.until(ec.visibility_of_element_located((By.CLASS_NAME, "popcov")))
    population_coverage_results_table = population_coverage.find_element(By.XPATH, '/html/body/div[3]/table[1]/tbody').text.splitlines()
    with open('results/population_coverage_graph.png', 'wb') as file:
        file.write(population_coverage.find_element(By.XPATH, '/html/body/div[3]/table[2]/tbody/tr[3]/td/img').screenshot_as_png)
    population_coverage.close()
    print("Population coverage done")
    result_in_text = ''
    for index, row in enumerate(population_coverage_results_table):
        new_row = row.split(' ')
        if index == 0:
            result_in_text = new_row[1] + ': ' + new_row[2] + '\n' + new_row[0] + '\t'
        elif index == 1:
            result_in_text = result_in_text + new_row[0][:-1] + '\t' + new_row[1][:-1] + '\t' + new_row[2][:-1] + '\n'
        elif index == 4:
            std_dev = new_row[0] + '_' + new_row[1]
            result_in_text = result_in_text + std_dev + '\t' + new_row[2] + '\t' + new_row[3] + '\t' + new_row[4]
        else:
            result_in_text = result_in_text + new_row[0] + '\t' + new_row[1] + '\t' + new_row[2] + '\t' + new_row[3] + '\n'
    pop_cov_results = open('results/pop_cov_results.txt', 'a')
    pop_cov_results.write(result_in_text)
    pop_cov_results.close()



    wait_immunogenicity_mhci = WebDriverWait(immunogenicity_mhci, 6000)
    wait_immunogenicity_mhci.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div[3]/table")))
    immunogenicity_mhci_results_table = immunogenicity_mhci.find_element(By.XPATH, '/html/body/div[3]/table/tbody').text.splitlines()
    immunogenicity_mhci.close()
    print("Immunogenicity done")
    immunogenicity_mhci_results_table = list(dict.fromkeys(immunogenicity_mhci_results_table))
    for row in immunogenicity_mhci_results_table:
        immunogenicity_mhci_result_row = row.split(' ')
        for i in range(immunogenicity_indexes):
            if immunogenicity_mhci_result_row[0] in analysis_results[i]:
                immunogenicity_score = round(float(immunogenicity_mhci_result_row[2]), 3)
                analysis_results[i].append(immunogenicity_score)



    def algpred_try_until_it_works(chunk_of_400, counter = 10):
        if counter == 0:
            print("Failed algpred after 10 tries")
            return

        algpred_seq_file = open('algpred_seq_file.txt', 'a')
        algpred_seq_file.write(chunk_of_400)
        algpred_seq_file.close()


        try:
            print("Runnning algpred2")
            algpred2 = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
            algpred2.get(algpred2_url)
            algpred2.find_element(By.XPATH, "/html/body/header/div[3]/section/form/table/tbody/tr/td/font/p/font[2]/input").send_keys(os.getcwd()+"/algpred_seq_file.txt")
            algpred2.find_element(By.XPATH, "/html/body/header/div[3]/section/form/table/tbody/tr/td/font/font/p[3]/font/font[2]/input[2]").click()
            wait_algpred2 = WebDriverWait(algpred2, 1200)
            wait_algpred2.until(ec.visibility_of_element_located((By.CLASS_NAME, "scrollable")))
        except:
            algpred2.close()
            os.remove(os.getcwd()+"/algpred_seq_file.txt")
            print("Retrying algpred ...")
            algpred_try_until_it_works(counter-1)
            return
        else:
            algpred2_results_table = algpred2.find_element(By.XPATH, '/html/body/header/div[3]/main/div/table[2]/tbody').text.splitlines()
            algpred2.close()
            print("Algpred worked!!!")
            os.remove(os.getcwd()+"/algpred_seq_file.txt")
            input_for_results = []
            for index, row in enumerate(algpred2_results_table):
                algpred2_result_row = row.split(' ')
                seq_nr = int(algpred2_result_row[0][3:])
                if index != len(algpred2_results_table) - 1:
                    algpred2_result_next_row = algpred2_results_table[index+1].split(' ')
                    seq_nr_next = int(algpred2_result_next_row[0][3:])
                    if seq_nr != seq_nr_next:
                        input_for_results.append(algpred2_result_row[-1])
                    else:
                        continue
                else:
                    input_for_results.append(algpred2_result_row[-1])
            return input_for_results

    algpred_all_results = []
    for seq_of_400_epitopes in algpred_chunks:
        algpred_results = algpred_try_until_it_works(seq_of_400_epitopes)
        if algpred_results is not None:
            for parameter in algpred_results:
                algpred_all_results.append(parameter)
        elif algpred_results is None:
            for empty_index in range(400):
                algpred_all_results.append("")

    for i in range(len(analysis_results)):
        analysis_results[i].append(algpred_all_results[i])



    def toxinpred_try_until_it_works(chunk_of_400, counter=10):
        if counter == 0:
            print("Failed after 10 tries")
            return
        try:
            print("Running toxinpred")
            toxinpred = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
            wait_toxinpred = WebDriverWait(toxinpred, 600)
            toxinpred.get(toxinpred_url)
            toxinpred.find_element(By.XPATH, "/html/body/table[2]/tbody/tr/td/form/fieldset/table[1]/tbody/tr[2]/td/textarea").send_keys(chunk_of_400)
            time.sleep(5)
            toxinpred.find_element(By.NAME, "checkAll").click()
        except:
            toxinpred.close()
            print("Retrying toxinpred...")
            toxinpred_try_until_it_works(chunk_of_400, counter-1)
            return
        else:
            print("Toxinpred worked!!!")
            toxinpred.find_element(By.XPATH, "/html/body/table[2]/tbody/tr/td/form/fieldset/table[2]/tbody/tr[3]/td/input[2]").click()

        def check_400(other_counter=10):
            if other_counter == 0:
                print("Failed selection after 10 tries")
                return

            try:
                toxinpred.find_element(By.XPATH, "/html/body/div[2]/table/tfoot/tr/td/select").click()
                toxinpred.find_element(By.XPATH, "/html/body/div[2]/table/tfoot/tr/td/select/option[8]").click()
            except:
                toxinpred.refresh()
                print("Retrying selection")
                check_400(other_counter-1)
                return
            else:
                toxinpred_results_table = toxinpred.find_element(By.XPATH, '/html/body/div[2]/table/tbody').text.splitlines()
                toxinpred.close()
                return toxinpred_results_table

        try:
            wait_toxinpred.until(ec.visibility_of_element_located((By.ID, "tableTwo")))
            time.sleep(5)
        except:
            toxinpred.close()
            print("Nothing detected. Retrying...")
            toxinpred_try_until_it_works(chunk_of_400, counter-1)
            return
        else:
            table_of_epitopes = check_400()
            return table_of_epitopes

    for index, seq_of_400_epitopes in enumerate(toxinpred_chunks):
        toxinpred_results_table = toxinpred_try_until_it_works(seq_of_400_epitopes)
        if toxinpred_results_table is not None:
            for results_index, row in enumerate(toxinpred_results_table):
                toxinpred_result_row = row.split(' ')
                toxinpred_row_included_indexes = [2,3,4,5,6,8,9,10]
                toxinpred_new_result_row = [toxinpred_result_row[x] for x in toxinpred_row_included_indexes]
                for parameter in toxinpred_new_result_row:
                    toxinpred_epitopes[400*index + results_index].append(parameter)
        elif toxinpred_results_table is None:
            for empty_index in range(400):
                for i in range(8):
                    toxinpred_epitopes[400*index + empty_index].append(None)

    none_list = []
    for i in range(8):
        none_list.append(None)
    for index in toxinpred_excluded_indexes:
        toxinpred_epitopes.insert(index, none_list)

    for alist in toxinpred_epitopes:
        del alist[0]

    for i in range(len(analysis_results)):
        for n in range(len(toxinpred_epitopes[i])):
            analysis_results[i].append(toxinpred_epitopes[i][n])



    os.remove(os.getcwd()+"/seq_file.txt")
    os.remove(os.getcwd()+"/pop_cov_file.txt")
    os.remove(os.getcwd()+"/immunogenicity_file.txt")
    os.remove(os.getcwd()+"/cluster_file.txt")
    os.remove(os.getcwd()+"/conservancy_seq_file.txt")
    os.remove(os.getcwd()+"/conservancy_protein_file.txt")

    return analysis_results



def make_csv_from_results(results_from_prediction, results_from_analysis):

    """The returned results from the analysis are added to the lists of prediction and saved as csv files for each different method."""

    columns_to_add = ['Mol_Weight', 'Isoelectric_Point', 'Aromaticity','Instability_Index','Helix_2_Struc', 'Turn_2_Struc', 'Sheet_2_Struc', 'Reduces_Cys', 'Disulfide_Bridge',
    'Hydropathicity', 'Charge_at_pH7', 'Antigenicity_Score', 'Antigen_Prediction', 'Allergen_Prediction','SVM_Score', 'Toxicity_Prediction', 'Hydrophobicity', 'Steric_hinderance',
    'Sidebulk', 'Amphipathicity', 'Hydrophilicity', 'Net_Hydrogen']

    for i in range(len(results_from_prediction[:-2])):
        if i == 0 or i == 1:
            for column_title in columns_to_add:
                results_from_prediction[i][0].append(column_title)
            results_from_prediction[i][0].insert(len(results_from_prediction[i][0])-9,'Immunogenicity_Score')
        else:
            for column_title in columns_to_add:
                results_from_prediction[i][0].append(column_title)

        if len(results_from_prediction[i]) > 1:
            length_to_cut = len(results_from_prediction[i][1:])
            for e in range(length_to_cut):
                for parameter in results_from_analysis[e][1:]:
                    results_from_prediction[i][e+1].append(parameter)
            results_from_analysis = results_from_analysis[length_to_cut:]

    for i in range(len(results_from_prediction)):
        if i == 0:
            f = open("results/mhci.csv", "w", newline="")
        if i == 1:
            f = open("results/mhci_processing.csv", "w", newline="")
        if i == 2:
            f = open("results/mhcii.csv", "w", newline="")
        if i == 3:
            f = open("results/bepipred2.csv", "w", newline="")
        if i == 4:
            f = open("results/bepipred.csv", "w", newline="")
        if i == 5:
            f = open("results/emini.csv", "w", newline="")
        if i == 6:
            f = open("results/chou_fasman.csv", "w", newline="")
        if i == 7:
            f = open("results/karplus_schulz.csv", "w", newline="")
        if i == 8:
            f = open("results/kolaskar_tonagonkar.csv", "w", newline="")
        if i == 9:
            f = open("results/parker.csv", "w", newline="")
        if i == 10:
            f = open("results/ellipro_linear.csv", "w", newline="")
        if i == 11:
            f = open("results/ellipro_discontinous.csv", "w", newline="")
        if i == 12:
            f = open("results/discotope.csv", "w", newline="")
        writer = csv.writer(f)
        writer.writerows(results_from_prediction[i])
        f.close()

# list_of_swissprot_ids = ['P59594', 'P0DTC2', 'K9N5Q8', 'P36334', 'Q0ZME7', 'P15423', 'Q6Q1S2', 'Q5MQD0', 'Q14EB0']
# prediction_results = read_prediction_results()
# analysis_input = make_inputs_for_analysis(prediction_results, list_of_swissprot_ids)
# analysis_results = analyse_all(analysis_input)
# make_csv_from_results(prediction_results, analysis_results)
