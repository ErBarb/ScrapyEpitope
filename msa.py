import os
from selenium import webdriver
from selenium.webdriver.common.by import By
import time
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as ec
import re

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


#alignment(list_of_sequences)


path_to_muscle = 'C:/Users/barbu/PycharmProjects/pythonProject/Pipeline/msa_results/muscle/muscle_spike.aln-fasta.fasta'
path_to_mafft = 'C:/Users/barbu/PycharmProjects/pythonProject/Pipeline/msa_results/mafft/mafft_spike.aln-fasta.fasta'

def get_conserved_sequences(path_to_file):
    gblocks_url = 'https://ngphylogeny.fr/tools/tool/276/form'
    driver = webdriver.Chrome(executable_path="chromedriver")
    driver.maximize_window()
    driver.get(gblocks_url)

    driver.find_element(By.XPATH, "/html/body/div/div[1]/div[1]/div/div/form/fieldset/div/div/div[1]/div/input").send_keys(path_to_file)
    driver.find_element(By.XPATH, "/html/body/div/div[1]/div[1]/div/div/form/input").click()
    wait = WebDriverWait(driver, 180)

    wait.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div/div[1]/div/div/div[6]/div/table/tbody/tr[2]/td[5]/a[5]")))
    driver.find_element(By.XPATH, "/html/body/div/div[1]/div/div/div[6]/div/table/tbody/tr[2]/td[5]/a[5]").click()

    wait.until(ec.visibility_of_element_located(
        (By.XPATH, "/html/body/div/div[1]/div[2]/pre/pre[3]/b[1]")))
    positions = driver.find_element(By.XPATH, '/html/body/div/div[1]/div[2]/pre/pre[3]').get_attribute("innerHTML").splitlines()[1]
    positions = positions[8:]
    list_positions = [int(s) for s in re.findall(r'\d+', positions)]

    driver.execute_script("window.history.go(-1)")
    wait.until(ec.visibility_of_element_located(
        (By.XPATH, "/html/body/div/div[1]/div/div/div[6]/div/table/tbody/tr[3]/td[4]/a[5]")))
    driver.find_element(By.XPATH, "/html/body/div/div[1]/div/div/div[6]/div/table/tbody/tr[3]/td[4]/a[5]").click()

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

    for i in conserved_sequences_dictionary:
        print(i, conserved_sequences_dictionary[i])

    #time.sleep(60)
    driver.close()


get_conserved_sequences(path_to_muscle)
get_conserved_sequences(path_to_mafft)

