from ast import Continue
from multiprocessing.sharedctypes import Value
import re
import os
from xmlrpc.client import ProtocolError
import requests
import time
import csv
from csv import reader
import pandas as pd
import matplotlib.pyplot as plt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from urllib import request
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as ec
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.keys import Keys



def read_prediction_results():

    """This function will read all the results from the csv files created by the prediction methods in prediction.py, remove the duplicate peptides,
    then return them in a tuple of lists"""

    df = pd.read_csv('epitope_prediction_results/mhci_epitopes.csv')
    df.drop_duplicates(subset=['peptide'], keep=False, inplace=True)
    mhci_prediction_list = [df.columns.tolist()] + df.values.tolist()
    #protein_ids = df['protein_id'].tolist()
    #print(protein_ids)

    df = pd.read_csv('epitope_prediction_results/mhci_proc_epitopes.csv')
    df.drop_duplicates(subset=['peptide'], keep=False, inplace=True)
    mhci_proc_prediction_list = [df.columns.tolist()] + df.values.tolist()

    df = pd.read_csv('epitope_prediction_results/mhcii_epitopes.csv')
    df.drop_duplicates(subset=['peptide'], keep=False, inplace=True)
    mhcii_prediction_list = [df.columns.tolist()] + df.values.tolist()

    df = pd.read_csv('epitope_prediction_results/bepipred2.0_epitopes.csv')
    df.drop_duplicates(subset=['predicted_epitope'], keep=False, inplace=True)
    bepipred2_prediction_list = [df.columns.tolist()] + df.values.tolist()

    df = pd.read_csv('epitope_prediction_results/bepipred1.0_epitopes.csv')
    df.drop_duplicates(subset=['predicted_epitope'], keep=False, inplace=True)
    bepipred_prediction_list = [df.columns.tolist()] + df.values.tolist()

    df = pd.read_csv('epitope_prediction_results/emini_epitopes.csv')
    df.drop_duplicates(subset=['predicted_epitope'], keep=False, inplace=True)
    emini_prediction_list = [df.columns.tolist()] + df.values.tolist()

    df = pd.read_csv('epitope_prediction_results/choufasman_epitopes.csv')
    df.drop_duplicates(subset=['predicted_epitope'], keep=False, inplace=True)
    choufasman_prediction_list = [df.columns.tolist()] + df.values.tolist()

    df = pd.read_csv('epitope_prediction_results/karplusschulz_epitopes.csv')
    df.drop_duplicates(subset=['predicted_epitope'], keep=False, inplace=True)
    karplusschulz_prediction_list = [df.columns.tolist()] + df.values.tolist()

    df = pd.read_csv('epitope_prediction_results/kolaskartongaonkar_epitopes.csv')
    df.drop_duplicates(subset=['predicted_epitope'], keep=False, inplace=True)
    kolaskartongaonkar_prediction_list = [df.columns.tolist()] + df.values.tolist()

    df = pd.read_csv('epitope_prediction_results/parker_epitopes.csv')
    df.drop_duplicates(subset=['predicted_epitope'], keep=False, inplace=True)
    parker_prediction_list = [df.columns.tolist()] + df.values.tolist()

    df = pd.read_csv('epitope_prediction_results/ellipro_linear_epitopes.csv')
    df.drop_duplicates(subset=['peptide'], keep=False, inplace=True)
    ellipro_linear_prediction_list = [df.columns.tolist()] + df.values.tolist()

    df = pd.read_csv('epitope_prediction_results/ellipro_discontinous_epitopes.csv')
    df.drop_duplicates(subset=['peptide'], keep=False, inplace=True)
    ellipro_discontinous_prediction_list = [df.columns.tolist()] + df.values.tolist()
    #print(ellipro_discontinous_prediction_list)

    df = pd.read_csv('epitope_prediction_results/discotope_epitopes.csv')
    df.drop_duplicates(subset=['peptide'], keep=False, inplace=True)
    discotope_prediction_list = [df.columns.tolist()] + df.values.tolist()
    #print(discotope_prediction_list)

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


    list_of_all_nonlinear_epitopes = []
    for index, value in enumerate(results_from_prediction[-2]):
        if index > 0:
            list_of_aa_with_pos = value[0].split(',')
            list_of_aa_without_pos = []
            for i in list_of_aa_with_pos:
                new_aa = i[0]
                list_of_aa_without_pos.append(new_aa)
            nonlinear_epitope = ''.join(list_of_aa_without_pos)
            list_of_all_nonlinear_epitopes.append(nonlinear_epitope)

    for index, value in enumerate(results_from_prediction[-1]):
        if index > 0:
            list_of_aa_with_pos = value[0].split(',')
            list_of_aa_without_pos = []
            for i in list_of_aa_with_pos:
                new_aa = i[0]
                list_of_aa_without_pos.append(new_aa)
            nonlinear_epitope = ''.join(list_of_aa_without_pos)
            list_of_all_nonlinear_epitopes.append(nonlinear_epitope)
    # print(list_of_all_nonlinear_epitopes)



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
            algpred_input = algpred_input + '>seq' + str(n + 1) + '\n' + algpred_400e_chunk[n] + '\n'
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
    final_directory = os.path.join(current_directory, r'analysis_results')
    if not os.path.exists(final_directory):
        os.makedirs(final_directory)

    print("Inputs finished")

    return list_of_all_linear_epitopes, toxinpred_chunks, toxinpred_excluded_indexes, immunogenicity_indexes, algpred_chunks, list_of_all_nonlinear_epitopes

def protparam(list_of_linear_epitopes):

    columns_to_add = ['peptide', 'mol_weight', 'isoelectric_point', 'aromaticity','instability_index','helix_2_struc', 'turn_2_struc', 'sheet_2_struc', 'reduces_cys', 'disulfide_bridge',
    'hydropathicity', 'charge_at_pH7']
    with open('analysis_results/protparam_analysis.csv', 'a') as f:
            writer = csv.writer(f)
            writer.writerow(columns_to_add)

    output_list = []
    for index, epitope in enumerate(list_of_linear_epitopes):
        row = []
        x = ProteinAnalysis(epitope)
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
        output_list.insert(index, row)
        row.insert(0,epitope)
        with open('analysis_results/protparam_analysis.csv', 'a') as f:
            writer = csv.writer(f)
            writer.writerow(row)

    print("Protparam analysis done")
    return output_list

def immunogenicity():

    options = Options()
    options.headless = True

    try:
        immunogenicity_mhci_url = 'http://tools.iedb.org/immunogenicity/'
        immunogenicity_mhci = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
        immunogenicity_mhci.get(immunogenicity_mhci_url)
        immunogenicity_mhci.find_element(By.NAME, "sequence_file").send_keys(os.getcwd()+"/immunogenicity_file.txt")
        immunogenicity_mhci.find_element(By.NAME, "submit").click()

        wait_immunogenicity_mhci = WebDriverWait(immunogenicity_mhci, 1200)
        wait_immunogenicity_mhci.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div[3]/table")))
        immunogenicity_mhci_results_table = immunogenicity_mhci.find_element(By.XPATH, '/html/body/div[3]/table/tbody').text.splitlines()
        immunogenicity_mhci.close()
        immunogenicity_mhci_results_table = list(dict.fromkeys(immunogenicity_mhci_results_table))

        print("Immunogenicity analysis done")
        return immunogenicity_mhci_results_table

    except:
        print("Immunogenicity analysis failed")

    

def vaxijen():

    columns_to_add = ['antigenicity_score', 'antigen_prediction']
    with open('analysis_results/vaxijen_analysis.csv', 'a') as f:
            writer = csv.writer(f)
            writer.writerow(columns_to_add)

    options = Options()
    options.headless = True
    try:
        vaxijen_url = 'http://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html'
        vaxijen = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
        vaxijen.get(vaxijen_url)
        time.sleep(5)
        vaxijen.find_element(By.XPATH, "/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[1]/td[2]/p/input").send_keys(os.getcwd()+"/seq_file.txt")
        vaxijen.find_element(By.XPATH, "/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[2]/td[2]/p/select/option[2]").click()
        vaxijen.find_element(By.XPATH, "/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[2]/td[3]/input").send_keys('0.5')
        vaxijen.find_element(By.XPATH, "/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[3]/td[2]/input[1]").click()

        wait_vaxijen = WebDriverWait(vaxijen, 1200)
        wait_vaxijen.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div/table/tbody/tr[3]/td[3]/font/b")))
        vaxijen_results_body = vaxijen.find_element(By.XPATH, '/html/body/div/table/tbody/tr[4]/td[3]/table/tbody').text.splitlines()
        vaxijen.close()
        
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

        with open('analysis_results/vaxijen_analysis.csv', 'a') as f:
            writer = csv.writer(f)
            writer.writerows(vaxijen_results)
        print("Vaxijen antigenicity analysis done")
        return vaxijen_results
    except:
        print("Vaxijen antigenicity analysis failed")

def cluster(list_of_linear_epitopes):

    options = Options()
    options.headless = True
    cluster_analysis_url = 'http://tools.iedb.org/cluster/'
    try:
        if len(list_of_linear_epitopes) < 3000:
            cluster = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
            cluster.maximize_window()
            cluster.get(cluster_analysis_url)
            cluster.find_element(By.NAME, "sequence_file").send_keys(os.getcwd()+"/cluster_file.txt")
            cluster.find_element(By.NAME, "submit").click()
            wait_cluster = WebDriverWait(cluster, 1200)
            wait_cluster.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div[3]/form/table/tbody/tr[2]/th")))
            cluster.find_element(By.NAME, "submit").click()
            wait_cluster.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div[3]/form/table/tbody")))
            cluster.find_element(By.NAME, "submit").click()

            wait_cluster.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div[3]/div/div[3]/div[2]")))
            cluster_results_table = cluster.find_element(By.XPATH, '/html/body/div[3]/div/div[3]/div[2]/div[1]/table/tbody').text.splitlines()
            cluster.close()
            
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
            with open('analysis_results/cluster_analysis.csv', 'w') as f:
                writer = csv.writer(f)
                writer.writerows(cluster_list_of_rows)
            print("Cluster analysis done")
        else:
            print("Too many epitopes for cluster analysis (>3000)")
    except:
        print("Cluster analysis failed")
    
def conservancy():

    options = Options()
    options.headless = True
    try:
        conservancy_analysis_url = 'http://tools.iedb.org/conservancy/'
        conservancy = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
        conservancy.maximize_window()
        conservancy.get(conservancy_analysis_url)
        conservancy.find_element(By.NAME, "epitope_file").send_keys(os.getcwd()+"/conservancy_seq_file.txt")
        conservancy.find_element(By.NAME, "protein_file").send_keys(os.getcwd()+"/conservancy_protein_file.txt")
        conservancy.find_element(By.NAME, "submit").click()

        wait_conservancy = WebDriverWait(conservancy, 6000)
        wait_conservancy.until(ec.visibility_of_element_located((By.ID, "result_table")))
        conservancy_results_columns = conservancy.find_element(By.XPATH, '/html/body/div[3]/table/thead/tr').text.splitlines()
        conservancy_results_table = conservancy.find_element(By.XPATH, '/html/body/div[3]/table/tbody').text.splitlines()
        conservancy.close()
        
        conservancy_list_of_rows = []
        conservancy_list_of_rows.append(conservancy_results_columns)
        for row in conservancy_results_table:
            conservancy_row = row.split(' ')
            conservancy_list_of_rows.append(conservancy_row)
        with open('analysis_results/conservancy_analysis.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerows(conservancy_list_of_rows)
        print("Conservancy analysis done")
    except:
        print("Conservancy analysis failed")

def try_population_coverage(counter=10):

    """The following function tries 10 times to make population coverage work, since it can run into errors. It uploads the file, then 
    selects 'World' option, select both MHC classes and clicks submit. It then waits for the results, saves the graph in a .png file and 
    returns the results table"""

    population_coverage_url = 'http://tools.iedb.org/population/'

    if counter == 0:
        print("Population coverage analysis failed  after 10 tries")
        return
    
    try:
        population_coverage = webdriver.Firefox(executable_path = '../ScrapyEpitope/geckodriver')
        population_coverage.maximize_window()
        population_coverage.get(population_coverage_url)
        population_coverage.find_element(By.ID, "id_epitope_allele_file").send_keys(os.getcwd()+"/pop_cov_file.txt")
        population_coverage.find_element(By.XPATH, "/html/body/div[3]/form/table[2]/tbody/tr[2]/td[1]/select/option[1]").click()
        population_coverage.find_element(By.XPATH, "/html/body/div[3]/form/table[2]/tbody/tr[4]/td/input").click()
        time.sleep(3)
        population_coverage.find_element(By.NAME, "submit").click()
        wait_population_coverage = WebDriverWait(population_coverage, 1200)
        wait_population_coverage.until(ec.visibility_of_element_located((By.CLASS_NAME, "popcov")))
        population_coverage_results_table = population_coverage.find_element(By.XPATH, '/html/body/div[3]/table[1]/tbody').text.splitlines()
        if population_coverage_results_table is None:
            try_population_coverage(counter-1)
        time.sleep(3)
        with open('analysis_results/population_coverage_graph.png', 'wb') as file:
            file.write(population_coverage.find_element(By.XPATH, '/html/body/div[3]/table[2]/tbody/tr[3]/td/img').screenshot_as_png)
        time.sleep(10)
        population_coverage.close()
        
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
        
        pop_cov_results = open('analysis_results/pop_cov_results.txt', 'a') 
        pop_cov_results.write(result_in_text)
        pop_cov_results.close()
        print("Population coverage analysis done")

    except:
        try:
            population_coverage.close()
        except:
            pass
        print("Retrying population coverage ...")
        try_population_coverage(counter-1)
        return
        
def algpred(algpred_chunks):

    def algpred_try_until_it_works(chunk_of_400, counter = 10):

        """The following function tries 10 times to make algpred2 work, since it can run into errors. It uploads the sequences in chunks of 400
        and clicks submit. It then waits for the results and saves them in their respective order. After all sequences are analysed it returns a
        list of lists with all the results"""

        if counter == 0:
            print("Algpred failed after 10 tries")
            return
        
        algpred_seq_file = open('algpred_seq_file.txt', 'a') 
        algpred_seq_file.write(chunk_of_400)
        algpred_seq_file.close()

        options = Options()
        options.headless = True
        algpred2_url = 'https://webs.iiitd.edu.in/raghava/algpred2/batch.html'

        try:
            algpred2 = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
            algpred2.get(algpred2_url)
            algpred2.find_element(By.XPATH, "/html/body/header/div[3]/section/form/table/tbody/tr/td/font/p/font[2]/input").send_keys(os.getcwd()+"/algpred_seq_file.txt")
            algpred2.find_element(By.XPATH, "/html/body/header/div[3]/section/form/table/tbody/tr/td/font/font/p[3]/font/font[2]/input[2]").click()
            wait_algpred2 = WebDriverWait(algpred2, 1200)
            wait_algpred2.until(ec.visibility_of_element_located((By.CLASS_NAME, "scrollable")))
            algpred2_results_table = algpred2.find_element(By.XPATH, '/html/body/header/div[3]/main/div/table[2]/tbody').text.splitlines()
            algpred2.close()
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

        except:
            try:
                algpred2.close()
            except:
                pass
            os.remove(os.getcwd()+"/algpred_seq_file.txt")
            print("Retrying algpred ...")
            algpred_try_until_it_works(counter-1)
            return
    
    algpred_all_results = []
    for seq_of_400_epitopes in algpred_chunks:
        algpred_results = algpred_try_until_it_works(seq_of_400_epitopes)
        if algpred_results is not None:
            for parameter in algpred_results:
                algpred_all_results.append(parameter)
        elif algpred_results is None:
            for empty_index in range(400):
                algpred_all_results.append(None)

    print("Algpred2 allergenicity analysis done")
    return algpred_all_results

def toxinpred(toxinpred_chunks, toxinpred_excluded_indexes):

    def toxinpred_try_until_it_works(chunk_of_400, counter=10):

        """The following function tries 10 times to make toxinpred work, since it can run into errors. It uploads the sequences in chunks of 400
        and clicks submit. It then waits for the results and saves them in their respective order. After all sequences are analysed it returns a
        list of lists with all the results"""

        if counter == 0:
            print("Toxinpred failed after 10 tries")
            return

        def check_400(other_counter=10):
            if other_counter == 0:
                return

            try:
                toxinpred.find_element(By.XPATH, "/html/body/div[2]/table/tfoot/tr/td/select").click()
                toxinpred.find_element(By.XPATH, "/html/body/div[2]/table/tfoot/tr/td/select/option[8]").click()
            except:
                toxinpred.refresh()
                check_400(other_counter-1)
                return
            else:
                toxinpred_results_table = toxinpred.find_element(By.XPATH, '/html/body/div[2]/table/tbody').text.splitlines()
                toxinpred.close()
                return toxinpred_results_table

        options = Options()
        options.headless = True
        toxinpred_url = 'https://webs.iiitd.edu.in/raghava/toxinpred/multi_submit.php'

        try:
            toxinpred = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
            wait_toxinpred = WebDriverWait(toxinpred, 1200)
            toxinpred.get(toxinpred_url)
            toxinpred.find_element(By.XPATH, "/html/body/table[2]/tbody/tr/td/form/fieldset/table[1]/tbody/tr[2]/td/textarea").send_keys(chunk_of_400)
            time.sleep(5)
            toxinpred.find_element(By.NAME, "checkAll").click()
            toxinpred.find_element(By.XPATH, "/html/body/table[2]/tbody/tr/td/form/fieldset/table[2]/tbody/tr[3]/td/input[2]").click()
            wait_toxinpred.until(ec.visibility_of_element_located((By.ID, "tableTwo")))
            time.sleep(5)
            table_of_epitopes = check_400()
            return table_of_epitopes
        except:
            try:
                toxinpred.close()
            except:
                pass
            print("Retrying toxinpred...")
            toxinpred_try_until_it_works(chunk_of_400, counter-1)
            return   
        
    # Gets toxinpred results and stores them in analysis_results:
    toxinpred_results = []
    for index, seq_of_400_epitopes in enumerate(toxinpred_chunks):
        toxinpred_results_table = toxinpred_try_until_it_works(seq_of_400_epitopes)
        if toxinpred_results_table is not None:
            for results_index, row in enumerate(toxinpred_results_table):
                toxinpred_result_row = row.split(' ')
                toxinpred_row_included_indexes = [2,3,4,5,6,8,9,10]
                toxinpred_new_result_row = [toxinpred_result_row[x] for x in toxinpred_row_included_indexes]
                toxinpred_results.append(toxinpred_new_result_row)
                # for parameter in toxinpred_new_result_row:
                #     toxinpred_results[400*index + results_index].append(parameter)
        elif toxinpred_results_table is None:
            empty_list = [None] * 8
            for i in range(400):
                #toxinpred_results[400*index + empty_index].append(empty_list)
                toxinpred_results.append(empty_list)
    
    print("Toxinpred toxicity analysis done")

    none_list = [None] * 8
    for index in toxinpred_excluded_indexes:
        toxinpred_results.insert(index, none_list)
    
    return toxinpred_results

def expasy_and_solubility(list_of_epitopes, linear = 'Yes'):

    expasy_results = []
    protein_sol_results = []

    columns_to_add = ['peptide', '(-)_charged_residues (asp+glu)','(+)_charged_residues (arg+lys)', 'half_life_hours (mammalian reticulocytes, in vitro)', 
    'aliphatic_index']
    with open('analysis_results/expasy_analysis.csv', 'a') as f:
        writer = csv.writer(f)
        writer.writerow(columns_to_add)

    aa_composition_results = [['peptide', 'Ala (A)', 'Arg (R)', 'Asn (N)', 'Asp (D)', 'Cys (C)', 'Gln (Q)', 'Glu (E)', 'Gly (G)', 'His (H)', 'Ile (I)', 'Leu (L)'
    , 'Lys (K)', 'Met (M)', 'Phe (F)', 'Pro (P)', 'Ser (S)', 'Thr (T)', 'Trp (W)', 'Tyr (Y)', 'Val (V)', 'Pyl (O)', 'Sec (U)']]
    atomic_composition_results = [['peptide', 'Carbon (C)', 'Hydrogen (H)', 'Nitrogen (N)', 'Oxygen (O)', 'Sulfur (S)']]
    with open('analysis_results/iter_aa_composition_analysis.csv', 'a') as f:
        writer = csv.writer(f)
        writer.writerow(aa_composition_results[0])
    with open('analysis_results/iter_atomic_composition_analysis.csv', 'a') as f:
        writer = csv.writer(f)
        writer.writerow(atomic_composition_results[0])

    options = Options()
    options.headless = True
    expasy_url = 'https://web.expasy.org/protparam/'
    # protein_sol_url = 'https://protein-sol.manchester.ac.uk/'

    for index, seq in enumerate(list_of_epitopes):
        try:
            print("Analysing sequence with expasy: ", seq)
            expasy = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
            wait_expasy = WebDriverWait(expasy, 60)
            expasy.get(expasy_url)
            expasy.find_element(By.XPATH, "/html/body/div[2]/div[2]/form/textarea").send_keys(seq)
            expasy.find_element(By.XPATH, "/html/body/div[2]/div[2]/form/p[1]/input[2]").click()

            # if len(seq) >= 21:
            #     protein_sol = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
            #     protein_sol.get(protein_sol_url)
            #     protein_sol.find_element(By.NAME, "sequence-input").send_keys(seq)
            #     protein_sol.find_element(By.NAME, "singleprediction").click()
            #     time.sleep(3)
            #     protein_sol_result = protein_sol.find_element(By.XPATH, '/html/body/div[3]/div/div[1]/p[2]').text
            #     protein_sol_results.append(protein_sol_result)
            #     protein_sol.close()
            # else:
            #     protein_sol_results.append(None)
            
            wait_expasy.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div[2]/div[2]/pre[2]/form/input[1]")))
            expasy_text = expasy.find_element(By.XPATH, '/html/body/div[2]/div[2]/pre[2]').text.splitlines()
            expasy.close()

            aa_composition = expasy_text[6:28]
            neg_residues = expasy_text[32]
            pos_residues = expasy_text[33]
            atomic_composition = expasy_text[37:42]
            half_life = expasy_text[56:64]
            aliphatic_index = expasy_text[-3]
            if len(expasy_text) == 75:
                half_life = expasy_text[56:64]
            elif len(expasy_text) == 71:
                half_life = expasy_text[52:60]
            elif len(expasy_text) == 72:
                half_life = expasy_text[53:61]
            elif len(expasy_text) == 79:
                half_life = expasy_text[60:68]
            elif len(expasy_text) == 76:
                half_life = expasy_text[57:65]


            aa_row = [seq]
            for i in range(len(aa_composition)):
                cell = aa_composition[i].split()
                aa_row.append(cell[2])
            aa_composition_results.append(aa_row)

            with open('analysis_results/iter_aa_composition_analysis.csv', 'a') as f:
                writer = csv.writer(f)
                writer.writerow(aa_row)

            atomic_row = [seq]
            for i in range(len(atomic_composition)):
                cell = atomic_composition[i].split()
                atomic_row.append(cell[2])
            atomic_composition_results.append(atomic_row)

            with open('analysis_results/iter_atomic_composition_analysis.csv', 'a') as f:
                writer = csv.writer(f)
                writer.writerow(aa_row)

            expasy_row = []
            expasy_row.append(neg_residues[-1])
            expasy_row.append(pos_residues[-1])
            expasy_row.append(half_life[4].split()[4])
            expasy_row.append(float(aliphatic_index.split()[2]))
            expasy_results.append(expasy_row)

            expasy_row.insert(0, seq)
            with open('analysis_results/expasy_analysis.csv', 'a') as f:
                writer = csv.writer(f)
                writer.writerow(expasy_row)
            
        except:
            print(seq + " expasy failed")

            none_list = [None] * 4
            expasy_results.append(none_list)

            none_atomic_row = [seq]
            for i in range(6):
                none_atomic_row.append(None)
            atomic_composition_results.append(none_atomic_row)

            none_aa_row = [seq]
            for i in range(23):
                none_aa_row.append(None)
            aa_composition_results.append(none_aa_row)
            continue
        
        

    # print("Expasy analysis done")
    # print("Protein solubility analysis done")

    if linear == 'Yes':
        print("Expasy analysis on linear epitopes done")
        with open('analysis_results/aa_composition.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerows(aa_composition_results)

        with open('analysis_results/atomic_composition.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerows(atomic_composition_results)

    elif linear == 'No':
        print("Expasy analysis on nonlinear epitopes done")
        with open('analysis_results/nonlinear_aa_composition.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerows(aa_composition_results)

        with open('analysis_results/nonlinear_atomic_composition.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerows(atomic_composition_results)

    return expasy_results, protein_sol_results


def pepstats(list_of_sequences):
    
    columns_to_add = ['peptide', 'tiny_aa_percentage', 'small_aa_percentage', 'aliphatic_aa_percentage', 'aromatic_aa_percentage', 'non_polar_aa_percentage', 'polar_aa_percentage', 
    'charged_aa_percentage', 'basic_aa_percentage', 'acidic_aa_percentage']
    with open('analysis_results/pepstats_analysis.csv', 'a') as f:
        writer = csv.writer(f)
        writer.writerow(columns_to_add)

    pepstats_results = []
    for seq in list_of_sequences:

        try:
            print("Analysing sequence with pepstats: ", seq)
            os.system(
                'python embosspepstats.py --email erald.bb@gmail.com --sequence ' + seq + ' --outfile pepstats_results --quiet')
            # python embosspepstats.py --email erald.bb@gmail.com --sequence SVDCNMYICGDSTEC --outfile pepstats_results --quiet
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

                pepstats_results.append(results_row)
            
            results_row.insert(0, seq)
            with open('analysis_results/pepstats_analysis.csv', 'a') as f:
                writer = csv.writer(f)
                writer.writerow(results_row)

            os.remove(os.getcwd()+"/pepstats_results.out.txt")
            os.remove(os.getcwd()+"/pepstats_results.sequence.txt")

        except:
            print(seq + " pepstats failed")
            none_pepstats_row = []
            for i in range(10):
                none_pepstats_row.append(None)
            pepstats_results.append(none_pepstats_row)
            continue

    return pepstats_results

def analyse_all(tuple_inputs):

    """This function uses Selenium to access the following tools: Toxinpred, Algpred2, Vaxijen and IEDB tools such Immunogenicity for
    MHCI, Population Coverage, Cluster and Conservancy analysis. It also uses Protparam module from Biopython to analyse each sequence.
    Population coverage, cluster and conservancy results are stored in files. The other ones are returned as list of lists."""

    options = Options()
    options.headless = True

    list_of_linear_epitopes = tuple_inputs[0]
    toxinpred_chunks = tuple_inputs[1]
    toxinpred_excluded_indexes = tuple_inputs[2]
    immunogenicity_indexes = tuple_inputs[3]
    algpred_chunks = tuple_inputs[4]
    list_of_all_nonlinear_epitopes = tuple_inputs[5]

    # toxinpred_url = 'https://webs.iiitd.edu.in/raghava/toxinpred/multi_submit.php'
    # algpred2_url = 'https://webs.iiitd.edu.in/raghava/algpred2/batch.html'
    # vaxijen_url = 'http://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html'
    # immunogenicity_mhci_url = 'http://tools.iedb.org/immunogenicity/'
    # expasy_url = 'https://web.expasy.org/protparam/'
    # population_coverage_url = 'http://tools.iedb.org/population/'
    # cluster_analysis_url = 'http://tools.iedb.org/cluster/'
    # conservancy_analysis_url = 'http://tools.iedb.org/conservancy/'
    # protein_sol_url = 'https://protein-sol.manchester.ac.uk/'

    analysis_results = []
    for index, value in enumerate(list_of_linear_epitopes):
        row = []
        row.insert(0,list_of_linear_epitopes[index])
        analysis_results.append(row)

    protparam_results = protparam(list_of_linear_epitopes)
    try:
        for i in range(len(analysis_results)):
            for n in range(len(protparam_results[i])):
                analysis_results[i].append(protparam_results[i][n])
    except:
        pass


    vaxijen_results = vaxijen()
    try:
        for i in range(len(vaxijen_results)):
            analysis_results[i].append(vaxijen_results[i][0])
            analysis_results[i].append(vaxijen_results[i][1])
    except:
        pass

    # immunogenicity_mhci_results_table = immunogenicity()
    # try:
    #     for row in immunogenicity_mhci_results_table:
    #         immunogenicity_mhci_result_row = row.split(' ')
    #         for i in range(immunogenicity_indexes):
    #             if immunogenicity_mhci_result_row[0] in analysis_results[i]:
    #                 immunogenicity_score = round(float(immunogenicity_mhci_result_row[2]), 3)
    #                 analysis_results[i].append(immunogenicity_score)
    # except:
    #     print("Immunogenicity failed")

    # cluster(list_of_linear_epitopes)
    # conservancy()
    # try_population_coverage()

    algpred_results = algpred(algpred_chunks)
    try:
        for i in range(len(analysis_results)):
            analysis_results[i].append(algpred_results[i])
    except:
        pass

    toxinpred_results = toxinpred(toxinpred_chunks, toxinpred_excluded_indexes)
    try:
        for i in range(len(analysis_results)):
            for n in range(len(toxinpred_results[i])):
                analysis_results[i].append(toxinpred_results[i][n])
    except:
        pass

    discontinous_expasy = expasy_and_solubility(list_of_all_nonlinear_epitopes, linear = 'No')
    expasy_linear_results = expasy_and_solubility(list_of_linear_epitopes)
    try:
        for i in range(len(analysis_results)):
            for n in range(len(expasy_linear_results[0][i])):
                analysis_results[i].append(expasy_linear_results[0][i][n])
    except:
        print("Expasy failed")

    # for i in range(len(analysis_results)):
    #     analysis_results[i].append(expasy_and_solubility_results[1][i])

    pepstats_results = pepstats(list_of_linear_epitopes)
    try:
        for i in range(len(analysis_results)):
            for n in range(len(pepstats_results[i])):
                analysis_results[i].append(pepstats_results[i][n])
        print("Pepstats analysis done")
    except:
        print("Pepstats failed")
    print(analysis_results)





    # Go to immunogenicity url, upload immunogenicity file and click submit:
    # immunogenicity_mhci = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
    # immunogenicity_mhci.get(immunogenicity_mhci_url)
    # immunogenicity_mhci.find_element(By.NAME, "sequence_file").send_keys(os.getcwd()+"/immunogenicity_file.txt")
    # immunogenicity_mhci.find_element(By.NAME, "submit").click()

    # Go to vaxijen url, upload sequence file, set option viruses, set threshold at 0.5 and click submit:
    # vaxijen = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
    # vaxijen.get(vaxijen_url)
    # vaxijen.find_element(By.XPATH, "/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[1]/td[2]/p/input").send_keys(os.getcwd()+"/seq_file.txt")
    # vaxijen.find_element(By.XPATH, "/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[2]/td[2]/p/select/option[2]").click()
    # vaxijen.find_element(By.XPATH, "/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[2]/td[3]/input").send_keys('0.5')
    # vaxijen.find_element(By.XPATH, "/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[3]/td[2]/input[1]").click()

    # Go to cluster analysis url, upload cluster file, click submit and doesnt change any other option on the next pages until results show:
    # if len(list_of_linear_epitopes) < 3000:
    #     cluster = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
    #     cluster.maximize_window()
    #     cluster.get(cluster_analysis_url)
    #     cluster.find_element(By.NAME, "sequence_file").send_keys(os.getcwd()+"/cluster_file.txt")
    #     cluster.find_element(By.NAME, "submit").click()
    #     wait_cluster = WebDriverWait(cluster, 6000)
    #     wait_cluster.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div[3]/form/table/tbody/tr[2]/th")))
    #     cluster.find_element(By.NAME, "submit").click()
    #     wait_cluster.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div[3]/form/table/tbody")))
    #     cluster.find_element(By.NAME, "submit").click()
    
    # Go to conservancy analysis url, upload protein and epitope files and click submit:
    # conservancy = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
    # conservancy.maximize_window()
    # conservancy.get(conservancy_analysis_url)
    # conservancy.find_element(By.NAME, "epitope_file").send_keys(os.getcwd()+"/conservancy_seq_file.txt")
    # conservancy.find_element(By.NAME, "protein_file").send_keys(os.getcwd()+"/conservancy_protein_file.txt")
    # conservancy.find_element(By.NAME, "submit").click()

    # Population coverage:
    # def retry_pop_cov(counter=10):

    #     """The following function tries 10 times to make population coverage work, since it can run into errors. It uploads the file, then 
    #     selects 'World' option, select both mhc classes and clicks submit. It then waits for the results, saves the graph in a .png file and 
    #     returns the results table"""

    #     if counter == 0:
    #         print("Population coverage analysis failed  after 10 tries")
    #         return
        
    #     try:
    #         population_coverage = webdriver.Firefox(executable_path = '../ScrapyEpitope/geckodriver')
    #         population_coverage.maximize_window()
    #         population_coverage.get(population_coverage_url)
    #         population_coverage.find_element(By.ID, "id_epitope_allele_file").send_keys(os.getcwd()+"/pop_cov_file.txt")
    #         population_coverage.find_element(By.XPATH, "/html/body/div[3]/form/table[2]/tbody/tr[2]/td[1]/select/option[1]").click()
    #         population_coverage.find_element(By.XPATH, "/html/body/div[3]/form/table[2]/tbody/tr[4]/td/input").click()
    #         time.sleep(3)
    #         population_coverage.find_element(By.NAME, "submit").click()
    #         wait_population_coverage = WebDriverWait(population_coverage, 6000)
    #         wait_population_coverage.until(ec.visibility_of_element_located((By.CLASS_NAME, "popcov")))
    #         population_coverage_results_table = population_coverage.find_element(By.XPATH, '/html/body/div[3]/table[1]/tbody').text.splitlines()
    #         if population_coverage_results_table is None:
    #             retry_pop_cov(counter-1)
    #     except:
    #         population_coverage.close()
    #         print("Retrying population coverage ...")
    #         retry_pop_cov(counter-1)
    #         return
    #     else:
    #         with open('results/population_coverage_graph.png', 'wb') as file:
    #             file.write(population_coverage.find_element(By.XPATH, '/html/body/div[3]/table[2]/tbody/tr[3]/td/img').screenshot_as_png)
    #         population_coverage.close()
    #         print("Population coverage done")
    #         return population_coverage_results_table


    # Waits for vaxijen to load and stores the results in analysis_results
    # wait_vaxijen = WebDriverWait(vaxijen, 6000)
    # wait_vaxijen.until(ec.visibility_of_element_located((By.CLASS_NAME, "boilerplate")))
    # vaxijen_results_body = vaxijen.find_element(By.XPATH, '/html/body/div/table/tbody/tr[4]/td[3]/table/tbody').text.splitlines()
    # vaxijen.close()
    # print("Vaxijen done")
    # vaxijen_results = []
    # for value in vaxijen_results_body:
    #     if value.startswith('Overall') == True:
    #         result = []
    #         antigenicity_val = float(re.search(r'[-+]?\d*\.*\d+', value).group())
    #         result.append(antigenicity_val)
    #         if antigenicity_val >= 0.5:
    #             result.append('Probable Antigen')
    #         else:
    #             result.append('Probable Non-Antigen')
    #         vaxijen_results.append(result)
    # for i in range(len(vaxijen_results)):
    #    analysis_results[i].append(vaxijen_results[i][0])
    #    analysis_results[i].append(vaxijen_results[i][1])


    # Waits for conservancy to load and stores the results in a .csv file
    # wait_conservancy = WebDriverWait(conservancy, 6000)
    # wait_conservancy.until(ec.visibility_of_element_located((By.ID, "result_table")))
    # conservancy_results_columns = conservancy.find_element(By.XPATH, '/html/body/div[3]/table/thead/tr').text.splitlines()
    # conservancy_results_table = conservancy.find_element(By.XPATH, '/html/body/div[3]/table/tbody').text.splitlines()
    # conservancy.close()
    # print("Conservancy done")
    # conservancy_list_of_rows = []
    # conservancy_list_of_rows.append(conservancy_results_columns)
    # for row in conservancy_results_table:
    #     conservancy_row = row.split(' ')
    #     conservancy_list_of_rows.append(conservancy_row)
    # with open('results/conservancy_analysis.csv', 'w') as f:
    #     writer = csv.writer(f)
    #     writer.writerows(conservancy_list_of_rows)


    # Waits for cluster to load and stores the results in a .csv file
    # if len(list_of_linear_epitopes) < 3000:
    #     wait_cluster.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div[3]/div/div[3]/div[2]")))
    #     cluster_results_table = cluster.find_element(By.XPATH, '/html/body/div[3]/div/div[3]/div[2]/div[1]/table/tbody').text.splitlines()
    #     cluster.close()
    #     print("Cluster done")
    #     cluster_list_of_rows = []
    #     cluster_columns = ['Cluster.Sub-Cluster Number','Peptide Number','Alignment','Position','Description','Peptide']
    #     cluster_list_of_rows.append(cluster_columns)
    #     for i in range(len(cluster_results_table)):
    #         if i == 0:
    #             continue
    #         else:
    #             cluster_row = cluster_results_table[i].split(' ')
    #             if len(cluster_row) == 6:    
    #                 cluster_list_of_rows.append(cluster_row)
    #             else:
    #                 indexes_to_remove = len(cluster_row) - 6
    #                 seq = cluster_row.pop(4)
    #                 for i in range(indexes_to_remove):
    #                     removed_duplicate = cluster_row.pop(4)
    #                     seq = seq + removed_duplicate
    #                 cluster_row.insert(4, seq)
    #                 cluster_list_of_rows.append(cluster_row)
    #     with open('results/cluster_analysis.csv', 'w') as f:
    #         writer = csv.writer(f)
    #         writer.writerows(cluster_list_of_rows)


    # Gets results from population coverage function, and stores them in a .txt file
    # population_coverage_results_table = retry_pop_cov()
    # if population_coverage_results_table is None:
    #     print("Population Coverage failed")
    # else:
    #     result_in_text = ''
    #     for index, row in enumerate(population_coverage_results_table):
    #         new_row = row.split(' ')
    #         if index == 0:
    #             result_in_text = new_row[1] + ': ' + new_row[2] + '\n' + new_row[0] + '\t'
    #         elif index == 1:
    #             result_in_text = result_in_text + new_row[0][:-1] + '\t' + new_row[1][:-1] + '\t' + new_row[2][:-1] + '\n'
    #         elif index == 4:
    #             std_dev = new_row[0] + '_' + new_row[1]
    #             result_in_text = result_in_text + std_dev + '\t' + new_row[2] + '\t' + new_row[3] + '\t' + new_row[4]
    #         else:
    #             result_in_text = result_in_text + new_row[0] + '\t' + new_row[1] + '\t' + new_row[2] + '\t' + new_row[3] + '\n'
    #     pop_cov_results = open('results/pop_cov_results.txt', 'a') 
    #     pop_cov_results.write(result_in_text)
    #     pop_cov_results.close()


    # Waits for immunogenicity to load and stores the results in analysis_results
    # wait_immunogenicity_mhci = WebDriverWait(immunogenicity_mhci, 6000)
    # wait_immunogenicity_mhci.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div[3]/table")))
    # immunogenicity_mhci_results_table = immunogenicity_mhci.find_element(By.XPATH, '/html/body/div[3]/table/tbody').text.splitlines()
    # immunogenicity_mhci.close()
    # print("Immunogenicity done")
    # immunogenicity_mhci_results_table = list(dict.fromkeys(immunogenicity_mhci_results_table))
    # for row in immunogenicity_mhci_results_table:
    #     immunogenicity_mhci_result_row = row.split(' ')
    #     for i in range(immunogenicity_indexes):
    #         if immunogenicity_mhci_result_row[0] in analysis_results[i]:
    #             immunogenicity_score = round(float(immunogenicity_mhci_result_row[2]), 3)
    #             analysis_results[i].append(immunogenicity_score)



    # def algpred_try_until_it_works(chunk_of_400, counter = 10):

    #     """The following function tries 10 times to make algpred2 work, since it can run into errors. It uploads the sequences in chunks of 400
    #     and clicks submit. It then waits for the results and saves them in their respective order. After all sequences are analysed it returns a
    #     list of lists with all the results"""

    #     if counter == 0:
    #         print("Algpred failed after 10 tries")
    #         return
        
    #     algpred_seq_file = open('algpred_seq_file.txt', 'a') 
    #     algpred_seq_file.write(chunk_of_400)
    #     algpred_seq_file.close()

    #     try:
    #         algpred2 = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
    #         algpred2.get(algpred2_url)
    #         algpred2.find_element(By.XPATH, "/html/body/header/div[3]/section/form/table/tbody/tr/td/font/p/font[2]/input").send_keys(os.getcwd()+"/algpred_seq_file.txt")
    #         algpred2.find_element(By.XPATH, "/html/body/header/div[3]/section/form/table/tbody/tr/td/font/font/p[3]/font/font[2]/input[2]").click()
    #         wait_algpred2 = WebDriverWait(algpred2, 1200)
    #         wait_algpred2.until(ec.visibility_of_element_located((By.CLASS_NAME, "scrollable")))
    #         algpred2_results_table = algpred2.find_element(By.XPATH, '/html/body/header/div[3]/main/div/table[2]/tbody').text.splitlines()
    #     except:
    #         algpred2.close()
    #         os.remove(os.getcwd()+"/algpred_seq_file.txt")
    #         print("Retrying algpred ...")
    #         algpred_try_until_it_works(counter-1)
    #         return
    #     else:    
    #         algpred2.close()
    #         os.remove(os.getcwd()+"/algpred_seq_file.txt")
    #         input_for_results = []
    #         for index, row in enumerate(algpred2_results_table):
    #             algpred2_result_row = row.split(' ')
    #             seq_nr = int(algpred2_result_row[0][3:])
    #             if index != len(algpred2_results_table) - 1:
    #                 algpred2_result_next_row = algpred2_results_table[index+1].split(' ')
    #                 seq_nr_next = int(algpred2_result_next_row[0][3:])
    #                 if seq_nr != seq_nr_next:
    #                     input_for_results.append(algpred2_result_row[-1])
    #                 else:
    #                     continue
    #             else:
    #                 input_for_results.append(algpred2_result_row[-1])
    #         return input_for_results
    
    # # The returned results of algpred2 are transfered into analysis_results
    # algpred_all_results = []
    # for seq_of_400_epitopes in algpred_chunks:
    #     algpred_results = algpred_try_until_it_works(seq_of_400_epitopes)
    #     if algpred_results is not None:
    #         for parameter in algpred_results:
    #             algpred_all_results.append(parameter)
    #     elif algpred_results is None:
    #         for empty_index in range(400):
    #             algpred_all_results.append(None)
    
    # print("Algpred2 done")

    # for i in range(len(analysis_results)):
    #     analysis_results[i].append(algpred_all_results[i])


    # def toxinpred_try_until_it_works(chunk_of_400, counter=10):

    #     """The following function tries 10 times to make toxinpred work, since it can run into errors. It uploads the sequences in chunks of 400
    #     and clicks submit. It then waits for the results and saves them in their respective order. After all sequences are analysed it returns a
    #     list of lists with all the results"""

    #     if counter == 0:
    #         print("Failed after 10 tries")
    #         return

    #     def check_400(other_counter=10):
    #         if other_counter == 0:
    #             return

    #         try:
    #             toxinpred.find_element(By.XPATH, "/html/body/div[2]/table/tfoot/tr/td/select").click()
    #             toxinpred.find_element(By.XPATH, "/html/body/div[2]/table/tfoot/tr/td/select/option[8]").click()
    #         except:
    #             toxinpred.refresh()
    #             check_400(other_counter-1)
    #             return
    #         else:
    #             toxinpred_results_table = toxinpred.find_element(By.XPATH, '/html/body/div[2]/table/tbody').text.splitlines()
    #             toxinpred.close()
    #             return toxinpred_results_table

    #     try:
    #         toxinpred = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
    #         wait_toxinpred = WebDriverWait(toxinpred, 600)
    #         toxinpred.get(toxinpred_url)
    #         toxinpred.find_element(By.XPATH, "/html/body/table[2]/tbody/tr/td/form/fieldset/table[1]/tbody/tr[2]/td/textarea").send_keys(chunk_of_400)
    #         time.sleep(5)
    #         toxinpred.find_element(By.NAME, "checkAll").click()
    #         toxinpred.find_element(By.XPATH, "/html/body/table[2]/tbody/tr/td/form/fieldset/table[2]/tbody/tr[3]/td/input[2]").click()
    #         wait_toxinpred.until(ec.visibility_of_element_located((By.ID, "tableTwo")))
    #         time.sleep(5)
    #     except:
    #         toxinpred.close()
    #         print("Retrying toxinpred...")
    #         toxinpred_try_until_it_works(chunk_of_400, counter-1)
    #         return
    #     else:
    #         table_of_epitopes = check_400()
    #         return table_of_epitopes
        
    # # Gets toxinpred results and stores them in analysis_results:
    # for index, seq_of_400_epitopes in enumerate(toxinpred_chunks):
    #     toxinpred_results_table = toxinpred_try_until_it_works(seq_of_400_epitopes)
    #     if toxinpred_results_table is not None:
    #         for results_index, row in enumerate(toxinpred_results_table):
    #             toxinpred_result_row = row.split(' ')
    #             toxinpred_row_included_indexes = [2,3,4,5,6,8,9,10]
    #             toxinpred_new_result_row = [toxinpred_result_row[x] for x in toxinpred_row_included_indexes]
    #             for parameter in toxinpred_new_result_row:
    #                 toxinpred_epitopes[400*index + results_index].append(parameter)
    #     elif toxinpred_results_table is None:
    #         empty_list = [None] * 9
    #         for empty_index in range(400):
    #             toxinpred_epitopes[400*index + empty_index].append(empty_list)
    
    # print("Toxinpred done")

    # none_list = [None] * 9
    # for index in toxinpred_excluded_indexes:
    #     toxinpred_epitopes.insert(index, none_list)
    
    # for alist in toxinpred_epitopes:
    #     del alist[0]

    # for i in range(len(analysis_results)):
    #     for n in range(len(toxinpred_epitopes[i])):
    #         analysis_results[i].append(toxinpred_epitopes[i][n])


    # Expasy and protein solubility:
    # expasy_results = []
    # protein_sol_results = []
    # aa_composition_results = [['peptide', 'Ala (A)', 'Arg (R)', 'Asn (N)', 'Asp (D)', 'Cys (C)', 'Gln (Q)', 'Glu (E)', 'Gly (G)', 'His (H)', 'Ile (I)', 'Leu (L)'
    # , 'Lys (K)', 'Met (M)', 'Phe (F)', 'Pro (P)', 'Ser (S)', 'Thr (T)', 'Trp (W)', 'Tyr (Y)', 'Val (V)', 'Pyl (O)', 'Sec (U)']]
    # atomic_composition_results = [['peptide', 'Carbon (C)', 'Hydrogen (H)', 'Nitrogen (N)', 'Oxygen (O)', 'Sulfur (S)']]

    # for index, seq in enumerate(list_of_linear_epitopes):
    #     try:
    #         expasy = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
    #         wait_expasy = WebDriverWait(expasy, 600)
    #         expasy.get(expasy_url)
    #         expasy.find_element(By.XPATH, "/html/body/div[2]/div[2]/form/textarea").send_keys(seq)
    #         expasy.find_element(By.XPATH, "/html/body/div[2]/div[2]/form/p[1]/input[2]").click()

    #         if len(seq) >= 21:
    #             protein_sol = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
    #             protein_sol.get(protein_sol_url)
    #             protein_sol.find_element(By.NAME, "sequence-input").send_keys(seq)
    #             protein_sol.find_element(By.NAME, "singleprediction").click()
    #             time.sleep(3)
    #             protein_sol_result = protein_sol.find_element(By.XPATH, '/html/body/div[3]/div/div[1]/p[2]').text
    #             protein_sol_results.append(protein_sol_result)
    #             protein_sol.close()
    #         else:
    #             protein_sol_results.append(None)
            
    #         wait_expasy.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div[2]/div[2]/pre[2]/form/input[1]")))
    #         expasy_text = expasy.find_element(By.XPATH, '/html/body/div[2]/div[2]/pre[2]').text.splitlines()
    #         expasy.close()
            
    #     except Exception as e:
    #         break
        
    #     aa_composition = expasy_text[6:28]
    #     neg_residues = expasy_text[32]
    #     pos_residues = expasy_text[33]
    #     atomic_composition = expasy_text[37:42]
    #     half_life = expasy_text[56:64]
    #     aliphatic_index = expasy_text[-3]
    #     if len(expasy_text) == 75:
    #         half_life = expasy_text[56:64]
    #     elif len(expasy_text) == 71:
    #         half_life = expasy_text[52:60]
    #     elif len(expasy_text) == 72:
    #         half_life = expasy_text[53:61]
    #     elif len(expasy_text) == 79:
    #         half_life = expasy_text[60:68]
    #     elif len(expasy_text) == 76:
    #         half_life = expasy_text[57:65]

    #     aa_row = [seq]
    #     for i in range(len(aa_composition)):
    #         cell = aa_composition[i].split()
    #         aa_row.append(cell[2])
    #     aa_composition_results.append(aa_row)

    #     atomic_row = [seq]
    #     for i in range(len(atomic_composition)):
    #         cell = atomic_composition[i].split()
    #         atomic_row.append(cell[2])
    #     atomic_composition_results.append(atomic_row)

    #     expasy_row = []
    #     expasy_row.append(neg_residues[-1])
    #     expasy_row.append(pos_residues[-1])
    #     expasy_row.append(half_life[4].split()[4])
    #     expasy_row.append(float(aliphatic_index.split()[2]))
    #     expasy_results.append(expasy_row)

    # print("Expasy analysis done")
    # print("Protein solubility analysis done")

    # for i in range(len(analysis_results)):
    #     for n in range(len(expasy_results[i])):
    #         analysis_results[i].append(expasy_results[i][n])

    # if len(analysis_results) == len(protein_sol_results):
    #     for i in range(len(analysis_results)):
    #         analysis_results[i].append(protein_sol_results[i])

    # with open('results/aa_composition.csv', 'w') as f:
    #     writer = csv.writer(f)
    #     writer.writerows(aa_composition_results)

    # with open('results/atomic_composition.csv', 'w') as f:
    #     writer = csv.writer(f)
    #     writer.writerows(atomic_composition_results)


    # Pepstats:
    # def pepstats(list_of_sequences):
        
    #     pepstats_results = []
    #     for seq in list_of_sequences:
    #         os.system(
    #             'python embosspepstats.py --email erald.bb@gmail.com --sequence ' + seq + ' --outfile pepstats_results --quiet')

    #         results_row = []

    #         with open('pepstats_results.out.txt') as f:
    #             lines = f.readlines()
    #             aa_properties = lines[-11:-1]
    #             exp_inclusion_bodies = float(lines[7].split()[-1])

    #             tiny_aa = float(aa_properties[1].split('\t')[-1][:-1])
    #             small_aa = float(aa_properties[2].split('\t')[-1][:-1])
    #             aliphatic_aa = float(aa_properties[3].split('\t')[-1][:-1])
    #             aromatic_aa = float(aa_properties[4].split('\t')[-1][:-1])
    #             non_polar_aa = float(aa_properties[5].split('\t')[-1][:-1])
    #             polar_aa = float(aa_properties[6].split('\t')[-1][:-1])
    #             charged_aa = float(aa_properties[7].split('\t')[-1][:-1])
    #             basic_aa = float(aa_properties[8].split('\t')[-1][:-1])
    #             acidic_aa = float(aa_properties[9].split('\t')[-1][:-1])
                
    #             results_row.append(tiny_aa)
    #             results_row.append(small_aa)
    #             results_row.append(aliphatic_aa)
    #             results_row.append(aromatic_aa)
    #             results_row.append(non_polar_aa)
    #             results_row.append(polar_aa)
    #             results_row.append(charged_aa)
    #             results_row.append(basic_aa)
    #             results_row.append(acidic_aa)
    #             results_row.append(exp_inclusion_bodies)

    #             pepstats_results.append(results_row)

    #         os.remove(os.getcwd()+"/pepstats_results.out.txt")
    #         os.remove(os.getcwd()+"/pepstats_results.sequence.txt")

    #     return pepstats_results

    # pepstats_results = pepstats(list_of_linear_epitopes)
    # if len(pepstats_results) == len(analysis_results):
    #     for i in range(len(analysis_results)):
    #         for n in range(len(pepstats_results[i])):
    #             analysis_results[i].append(pepstats_results[i][n])
    #     print("Pepstats analysis done")
    # else:
    #     print("Not all epitopes analysed with Pepstats")



    os.remove(os.getcwd()+"/seq_file.txt")
    os.remove(os.getcwd()+"/pop_cov_file.txt")
    os.remove(os.getcwd()+"/immunogenicity_file.txt")
    os.remove(os.getcwd()+"/cluster_file.txt")
    os.remove(os.getcwd()+"/conservancy_seq_file.txt")
    os.remove(os.getcwd()+"/conservancy_protein_file.txt")

    return analysis_results



def make_csv_from_results(results_from_prediction, results_from_analysis):

    """The returned results from the analysis are added to the lists of prediction and saved as csv files for each different method."""

    columns_to_add = ['peptide', 'mol_weight', 'isoelectric_point', 'aromaticity','instability_index','helix_2_struc', 'turn_2_struc', 'sheet_2_struc', 'reduces_cys', 'disulfide_bridge',
    'hydropathicity', 'charge_at_pH7', 'antigenicity_score', 'antigen_prediction', 'allergen_prediction','toxinpred_svm_score', 'toxicity_prediction', 'hydrophobicity', 'steric_hinderance',
    'sidebulk', 'amphipathicity', 'hydrophilicity', 'net_hydrogen','(-)_charged_residues (asp+glu)','(+)_charged_residues (arg+lys)', 'half_life_hours (mammalian reticulocytes, in vitro)', 
    'aliphatic_index', 'tiny_aa_percentage', 'small_aa_percentage', 'aliphatic_aa_percentage', 'aromatic_aa_percentage', 'non_polar_aa_percentage', 'polar_aa_percentage', 
    'charged_aa_percentage', 'basic_aa_percentage', 'acidic_aa_percentage', 'improbability_of_expression_in_inclusion_bodies']

    # columns_to_add = ['peptide', 'mol_weight', 'isoelectric_point', 'aromaticity','instability_index','helix_2_struc', 'turn_2_struc', 'sheet_2_struc', 'reduces_cys', 'disulfide_bridge',
    # 'hydropathicity', 'charge_at_pH7', 'antigenicity_score', 'antigen_prediction', 'allergen_prediction','toxinpred_svm_score', 'toxicity_prediction', 'hydrophobicity', 'steric_hinderance',
    # 'sidebulk', 'amphipathicity', 'hydrophilicity', 'net_hydrogen','(-)_charged_residues (asp+glu)','(+)_charged_residues (arg+lys)', 'half_life_hours (mammalian reticulocytes, in vitro)', 
    # 'aliphatic_index', 'solubility', 'tiny_aa_percentage', 'small_aa_percentage', 'aliphatic_aa_percentage', 'aromatic_aa_percentage', 'non_polar_aa_percentage', 'polar_aa_percentage', 
    # 'charged_aa_percentage', 'basic_aa_percentage', 'acidic_aa_percentage', 'improbability_of_expression_in_inclusion_bodies']

    # Divide results_from_analysis into 11 lists
    results_lists = []
    for i in range(len(results_from_prediction)-2):
        list_to_append = []
        if i == 0 or i == 1:
            list_to_append.append(columns_to_add[:])
            list_to_append[0].insert(14,'Immunogenicity_Score')
        else:
            list_to_append.append(columns_to_add)

        if len(results_from_prediction[i]) > 1:
            lengths_to_cut = len(results_from_prediction[i])-1
            rows_for_list = results_from_analysis[:lengths_to_cut]
            results_from_analysis = results_from_analysis[lengths_to_cut:]
            list_to_append.extend(rows_for_list)
            results_lists.append(list_to_append)
        else:
            results_lists.append(list_to_append)

    for i in range(len(results_lists)):
        if i == 0:
            f = open("analysis_results/mhci.csv", "w", newline="")
        if i == 1:
            f = open("analysis_results/mhci_processing.csv", "w", newline="")
        if i == 2:
            f = open("analysis_results/mhcii.csv", "w", newline="")
        if i == 3:
            f = open("analysis_results/bepipred2.csv", "w", newline="")
        if i == 4:
            f = open("analysis_results/bepipred.csv", "w", newline="")
        if i == 5:
            f = open("analysis_results/emini.csv", "w", newline="")
        if i == 6:
            f = open("analysis_results/chou_fasman.csv", "w", newline="")
        if i == 7:
            f = open("analysis_results/karplus_schulz.csv", "w", newline="")
        if i == 8:
            f = open("analysis_results/kolaskar_tonagonkar.csv", "w", newline="")
        if i == 9:
            f = open("analysis_results/parker.csv", "w", newline="")
        if i == 10:
            f = open("analysis_results/ellipro_linear.csv", "w", newline="")
        # if i == 11:
        #     f = open("analysis_results/ellipro_discontinous.csv", "w", newline="")
        # if i == 12:
        #     f = open("analysis_results/discotope.csv", "w", newline="")
        writer = csv.writer(f)
        writer.writerows(results_lists[i])
        f.close()

# list_of_swissprot_ids = ['P59594', 'P0DTC2', 'K9N5Q8', 'P36334', 'Q0ZME7', 'P15423', 'Q6Q1S2', 'Q5MQD0', 'Q14EB0']
# prediction_results = read_prediction_results()
# analysis_input = make_inputs_for_analysis(prediction_results, list_of_swissprot_ids)
# analysis_results = analyse_all(analysis_input)
# make_csv_from_results(prediction_results, analysis_results)
