import os
import re
import requests
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as ec
from selenium.webdriver.firefox.options import Options

def prediction_choice():

    """"""

    prediction_choice = input("""The following options are available for epitope prediction:\n
        1. Analyse and predict epitopes from the whole protein sequences\n
        2. Perform multiple sequence analysis and predict epitopes only from the conserved sequences\n
    Please enter 1 or 2\n""")

    if prediction_choice != "1" or prediction_choice != "2":
        print("Please enter either 1 or 2")
        prediction_choice()

    return prediction_choice


def alignment(seq_list, matrix=None, gapopen=None, gapext=None, order=None, nbtree=None, treeout=None, maxiterate=None, ffts=None):

    """This function uses REST API services from https://www.ebi.ac.uk/ to output alignment files using the multiple
    sequence alignment methods "MAFFT" and "MUSCLE". The argument is a list of swissprot protein IDs. The files are
    saved in their respective folders in msa_results"""

    seq_string = ''
    for i in seq_list:
        seq_string = seq_string + 'sp:' + i + ','
    seq_string = seq_string[:-1]

    os.system(
        'python msa_algos/mafft.py --email erald.bb@gmail.com --stype protein --sequence ' + seq_string + ' --outfile msa_results/mafft/mafft --format fasta --matrix ' + matrix + ' --gapopen ' + str(gapopen) + ' --gapext ' + str(gapext) + ' --order ' + order + ' --nbtree ' + str(nbtree) + ' --treeout ' + treeout + ' --maxiterate ' + str(maxiterate) + ' --ffts ' + ffts)
    os.system(
        'python msa_algos/muscle.py --email erald.bb@gmail.com --sequence ' + seq_string + ' --outfile msa_results/muscle/muscle --format fasta')



def get_conserved_sequences(path_to_file, min_seq_conserved_pos = None, min_seq_flank_pos = None, max_contigous_nonconserved_pos = None, min_length_block = None, allowed_gap_pos = None):

    """This function uploads the alignment file and sends the parameters to the Gblocks NGPhylogeny website,
    then returns a dictionary with the conserved regions of each protein"""

    options = Options()
    options.headless = True

    gblocks_url = 'https://ngphylogeny.fr/tools/tool/276/form'
    driver = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
    driver.maximize_window()
    driver.get(gblocks_url)

    upload_file = driver.find_element(By.XPATH, "/html/body/div/div[1]/div[1]/div/div/form/fieldset/div/div/div[1]/div/input")
    upload_file.send_keys(path_to_file)

    data_type = driver.find_element(By.XPATH, "/html/body/div/div[1]/div[1]/div/div/form/fieldset/div/div/div[2]/div[1]/div/select")
    data_type.send_keys('Protein')

    b1 = driver.find_element(By.XPATH, "/html/body/div/div[1]/div[1]/div/div/form/fieldset/div/div/div[3]/div/input")
    b1.clear()
    b1.send_keys(min_seq_conserved_pos)

    b2 = driver.find_element(By.XPATH, "/html/body/div/div[1]/div[1]/div/div/form/fieldset/div/div/div[4]/div/input")
    b2.clear()
    b2.send_keys(min_seq_flank_pos)

    b3 = driver.find_element(By.XPATH, "/html/body/div/div[1]/div[1]/div/div/form/fieldset/div/div/div[5]/div/input")
    b3.clear()
    b3.send_keys(max_contigous_nonconserved_pos)

    b4 = driver.find_element(By.XPATH, "/html/body/div/div[1]/div[1]/div/div/form/fieldset/div/div/div[6]/div/input")
    b4.clear()
    b4.send_keys(min_length_block)

    b5 = driver.find_element(By.XPATH, "/html/body/div/div[1]/div[1]/div/div/form/fieldset/div/div/div[7]/div/select")
    b5.send_keys(allowed_gap_pos)

    submit = driver.find_element(By.XPATH, "/html/body/div/div[1]/div[1]/div/div/form/input")
    submit.click()
    wait = WebDriverWait(driver, 180)

    wait.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div/div[1]/div/div/div[6]/div/table/tbody/tr[2]/td[5]/a[5]")))
    sequence_info_html = driver.find_element(By.XPATH, "/html/body/div/div[1]/div/div/div[6]/div/table/tbody/tr[2]/td[5]/a[5]")
    sequence_info_html.click()

    wait.until(ec.visibility_of_element_located(
        (By.XPATH, "/html/body/div/div[1]/div[2]/pre/pre[3]/b[1]")))
    positions = driver.find_element(By.XPATH, '/html/body/div/div[1]/div[2]/pre/pre[3]').get_attribute("innerHTML").splitlines()[1]
    positions = positions[8:]
    list_positions = [int(s) for s in re.findall(r'\d+', positions)]

    driver.execute_script("window.history.go(-1)")
    wait.until(ec.visibility_of_element_located(
        (By.XPATH, "/html/body/div/div[1]/div/div/div[6]/div/table/tbody/tr[3]/td[4]/a[5]")))
    cleaned_sequences_fasta = driver.find_element(By.XPATH, "/html/body/div/div[1]/div/div/div[6]/div/table/tbody/tr[3]/td[4]/a[5]")
    cleaned_sequences_fasta.click()

    wait.until(ec.visibility_of_element_located(
        (By.XPATH, "/html/body/div/div[1]/div[2]/pre")))
    cleaned_seq = driver.find_element(By.XPATH, '/html/body/div/div[1]/div[2]/pre').text
    cleaned_seq_list = cleaned_seq.split('>')
    cleaned_seq_list.pop(0)

    conserved_sequences_dictionary = {}
    for i in range(len(cleaned_seq_list)):
        lines = cleaned_seq_list[i].splitlines()
        protein_id = lines[0].split(' ')[1]
        cleaned_fasta_all = ''.join(lines[1:])
        list_of_fasta = []
        for i in range(1, len(list_positions), 2):
            cut_point = list_positions[i] - list_positions[i - 1] + 1
            cons_seq = cleaned_fasta_all[:cut_point]
            list_of_fasta.append(cons_seq)
            cleaned_fasta_all = cleaned_fasta_all[cut_point:]
        conserved_sequences_dictionary[protein_id] = list_of_fasta
    driver.close()

    return conserved_sequences_dictionary


def get_protein_sequences(list_of_swissprot_ids):
    
    """"""
    
    protein_dict = {}

    for id in list_of_swissprot_ids:
        baseUrl="http://www.uniprot.org/uniprot/"
        currentUrl=baseUrl+id+".fasta"
        response = requests.post(currentUrl)
        response_lines = response.text.split("\n")[1:]
        cData=''.join(response_lines)
        protein_dict[id] = (cData)

    return protein_dict

list_of_swissprot_ids = ['P59594', 'P0DTC2', 'K9N5Q8', 'P36334', 'Q0ZME7', 'P15423', 'Q6Q1S2', 'Q5MQD0', 'Q14EB0']
get_protein_sequences(list_of_swissprot_ids)