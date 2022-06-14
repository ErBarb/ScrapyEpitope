from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as ec
from selenium.webdriver.firefox.options import Options
import time
import os
import pandas as pd
from msa import alignment
from msa import get_conserved_sequences
import requests
import csv
import re

# import scrapy
# from scrapy import Spider
# from scrapy.http import FormRequest
# from scrapy.utils.project import get_project_settings
# from scrapy.utils.log import configure_logging
# from scrapy.crawler import CrawlerRunner
# from twisted.internet import reactor, defer

def get_pdb_from_swissprot(list_of_swissprot_ids):

    """"""

    options = Options()
    options.headless = True

    pdb_ids = []
    for id in list_of_swissprot_ids:
        try:
            uniprot_url = 'https://www.uniprot.org/uniprot/' + id
            uniprot = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
            uniprot.get(uniprot_url)
            wait_uniprot = WebDriverWait(uniprot, 30)
            wait_uniprot.until(ec.visibility_of_element_located((By.XPATH, '//div[@id="structure"]/protvista-uniprot-structure/div/div[2]/protvista-datatable/table/tbody')))
            uniprot_pdb_ids = uniprot.find_element(By.XPATH, '//div[@id="structure"]/protvista-uniprot-structure/div/div[2]/protvista-datatable/table/tbody').text.splitlines()
            for row in uniprot_pdb_ids:
                newrow = row.split()
                if (newrow[2] == 'EM' or newrow[2] == 'X-ray') and (newrow[6].startswith('1-') == True):
                    pdb_ids.append(newrow[1])
                    break
            else:
                if (newrow[2] == 'EM' or newrow[2] == 'X-ray') and re.search("^1.-", newrow[6]):
                    pdb_ids.append(newrow[1])
            uniprot.close()
        except:
            print("No PDB format for: " + id)
            uniprot.close()
    
    return pdb_ids

def mhci(conserved_sequences_dict, list_of_alleles, list_of_lengths):
    
    """This function uses the REST API from IEDB to access to MHC I Binding tool and predict the MHC Class I epitopes of the given conserved sequences. 
    It also needs a list of alleles and their respective lengths to run. It employs different methods to predict MHC Class I epitopes, including a 
    consensus approach which combines ANN, SMM, Comblib, NetMHCpan and Consensus. The results are returned as a list of lists"""

    alleles = ",".join(list_of_alleles)
    converted_list = [str(element) for element in list_of_lengths]
    lengths = ",".join(converted_list)

    mhci_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["allele"] + ["seq_num"] + ["start"] + ["end"] + ["length"] + \
        ["peptide"] + ["core"] + ["icore"] + ["score"] + ["percentile_rank"]
    mhci_results.append(columns)

    for key in conserved_sequences_dict:
        print("Predicting linear epitopes of protein " + key + " with MHCI")
        for conserved_sequence in conserved_sequences_dict[key]:
            if len(conserved_sequence) >= max(list_of_lengths):
                
                headers = {
                    'Content-Type': 'application/x-www-form-urlencoded',
                }
                
                data = 'method=recommended&sequence_text=' + conserved_sequence + '&allele=' + alleles + '&length=' + lengths
                
                response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/mhci/', headers=headers, data=data)
                response_data_split_by_line = response.content.decode('utf-8').splitlines()
                
                response_body = []
                for line in response_data_split_by_line:
                    split_line = line.split("\t")
                    response_body.append(split_line)

                try:
                    df = pd.DataFrame(response_body[1:], columns=response_body[0])

                except:

                    try:
                        headers = {
                            'Content-Type': 'application/x-www-form-urlencoded',
                        }
                
                        data = 'method=recommended&sequence_text=' + conserved_sequence + '&allele=HLA-A*01:01,HLA-A*02:01&length=8,9'
                
                        response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/mhci/', headers=headers, data=data)
                        response_data_split_by_line = response.content.decode('utf-8').splitlines()
                
                        response_body = []
                        for line in response_data_split_by_line:
                            split_line = line.split("\t")
                            response_body.append(split_line)

                        df = pd.DataFrame(response_body[1:], columns=response_body[0])
                        df["percentile_rank"] = pd.to_numeric(df["percentile_rank"], errors='coerce')
                        df = df.loc[df['percentile_rank'] <= 10]

                        if df.empty == True:
                            continue
                        else:
                            rows = [[i for i in row[1:]] for row in df.itertuples()]
                            for i in rows:
                                i = [key] + [conserved_sequence] + i
                                mhci_results.append(i)
                    except:
                        continue

                else:
                    df = pd.DataFrame(response_body[1:], columns=response_body[0])
                    df["percentile_rank"] = pd.to_numeric(df["percentile_rank"], errors='coerce')
                    df = df.loc[df['percentile_rank'] <= 1]

                    if df.empty == True:
                        continue
                    else:
                        rows = [[i for i in row[1:]] for row in df.itertuples()]
                        for i in rows:
                            i = [key] + [conserved_sequence] + i
                            mhci_results.append(i)
    with open('epitope_prediction_results/mhci_epitopes.csv', 'w', newline="") as f: 
        writer = csv.writer(f)
        writer.writerows(mhci_results)
    print("MHCI prediction done\n")

def mhci_proc(conserved_sequences_dict, list_of_alleles, list_of_lengths):

    """This function uses the REST API from IEDB to access to MHC I Processing tool and predict the MHC Class I epitopes of the given conserved sequences. 
    It also needs a list of alleles and their respective lengths to run. It combines predictors of proteasomal processing, TAP transport, and MHC binding 
    to produce an overall score for each peptide's intrinsic potential of being a T cell epitope. The results are returned as a list of lists"""


    alleles = ",".join(list_of_alleles)
    converted_list = [str(element) for element in list_of_lengths]
    lengths = ",".join(converted_list)

    mhci_proc_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["allele"] + ['#'] + ["start"] + ["end"] + ['peptide_length'] + \
        ["peptide"] + ["proteasome_score"] + ["tap_score"] + ["mhci_score"] + ["processing_score"] + ["total_score"] + ["mhci_ic50"]
    mhci_proc_results.append(columns)

    for key in conserved_sequences_dict:
        print("Predicting linear epitopes of protein " + key + " with MHCI Processing")
        for conserved_sequence in conserved_sequences_dict[key]:
            if len(conserved_sequence) >= max(list_of_lengths):

                headers = {
                    'Content-Type': 'application/x-www-form-urlencoded',
                }
                
                data = 'method=recommended&sequence_text=' + conserved_sequence + '&allele=' + alleles + '&length=' + lengths
                
                response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/processing/', headers=headers, data=data)
                response_data_split_by_line = response.content.decode('utf-8').splitlines()
                
                response_body = []
                for line in response_data_split_by_line:
                    split_line = line.split("\t")
                    response_body.append(split_line)
                
                
                df = pd.DataFrame(response_body[1:], columns=response_body[0])
                df["total_score"] = pd.to_numeric(df["total_score"])
                df = df.loc[df['total_score'] > 0]

                if df.empty == True:
                    continue
                else:
                    rows = [[i for i in row[1:]] for row in df.itertuples()]
                    for i in rows:
                        i = [key] + [conserved_sequence] + i
                        mhci_proc_results.append(i)
    with open('epitope_prediction_results/mhci_proc_epitopes.csv', 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(mhci_proc_results)                
    print("MHCI Processing prediction done\n")

def mhcii(conserved_sequences_dict, list_of_alleles, list_of_lengths):

    """This function uses the REST API from IEDB to access to MHC II Binding tool and predict the MHC Class II epitopes of the given conserved sequences. 
    It also needs a list of alleles and their respective lengths to run. It employs different methods to predict MHC Class II epitopes, including a consensus 
    approach which combines NN-align, SMM-align and Combinatorial library methods. The results are returned as a list of lists"""


    alleles = ",".join(list_of_alleles)
    converted_list = [str(element) for element in list_of_lengths]
    lengths = ",".join(converted_list)

    mhcii_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["allele"] + ["seq_num"] + ["start"] + ["end"] + ["length"] + \
        ["method"] + ["peptide"] + ["percentile_rank"] + ["adjusted_rank"] + \
        ["comblib_core"] + ['comblib_score'] + ['comblib_rank'] + ["comblib_adjusted_rank"] + \
        ["smm_align_core"] + ["smm_align_ic50"] + ["smm_align_rank"] + ["smm_align_adjusted_rank"] + \
        ['nn_align_core'] + ['nn_align_ic50'] + ["nn_align_rank"] + ["nn_align_adjusted_rank"] + \
        ["netmhciipan_core"] + ["netmhciipan_ic50"] + ["netmhciipan_rank"] + ['netmhciipan_adjusted_rank'] + \
        ['sturniolo_core'] + ["sturniolo_score"] + ["sturniolo_rank"] + ["sturniolo_adjusted_rank"]
    
    mhcii_results.append(columns)

    for key in conserved_sequences_dict:
        print("Predicting linear epitopes of protein " + key + " with MHCII")
        for conserved_sequence in conserved_sequences_dict[key]:
            if len(conserved_sequence) >= max(list_of_lengths):

                headers = {
                    'Content-Type': 'application/x-www-form-urlencoded',
                }
                
                data = 'method=recommended&sequence_text=' + conserved_sequence + '&allele=' + alleles + '&length=' + lengths
                
                response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/mhcii/', headers=headers, data=data)
                response_data_split_by_line = response.content.decode('utf-8').splitlines()

                response_body = []
                for line in response_data_split_by_line:
                    split_line = line.split("\t")
                    response_body.append(split_line)
                
                df = pd.DataFrame(response_body[1:], columns=response_body[0])
                df["adjusted_rank"] = pd.to_numeric(df["adjusted_rank"])
                df = df.loc[df['adjusted_rank'] <= 1]

                if df.empty == True:
                    continue
                else:
                    rows = [[i for i in row[1:]] for row in df.itertuples()]
                    for i in rows:
                        i = [key] + [conserved_sequence] + i
                        mhcii_results.append(i)
    with open('epitope_prediction_results/mhcii_epitopes.csv', 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(mhcii_results)   
    print("MHCII prediction done\n")

def bepipred2(conserved_sequences_dict):

    """This function uses the REST API from IEDB to access to Linear Antigen Prediction method Bepipred 2.0 and predict the linear epitopes of the given conserved 
    sequences using a Random Forest algorithm trained on epitopes and non-epitope amino acids determined from crystal structures. The results are returned as a 
    list of lists"""

    bepipred2_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["predicted_epitope"] + ["start_position"] + ["end_position"]
    bepipred2_results.append(columns)

    for key in conserved_sequences_dict:
        print("Predicting linear epitopes of protein " + key + " with Bepipred 2.0")
        for conserved_sequence in conserved_sequences_dict[key]:
            
            headers = {
                'Content-Type': 'application/x-www-form-urlencoded',
                }
            
            data = {
                'method': 'Bepipred-2.0',
                'sequence_text': conserved_sequence,
                'window_size': '7',
                }
            
            response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/bcell/', headers=headers, data=data)
            response_data_split_by_line = response.content.decode('utf-8').splitlines()
            
            response_body = []
            for line in response_data_split_by_line:
                split_line = line.split("\t")
                response_body.append(split_line)
            
            df = pd.DataFrame(response_body[1:], columns=response_body[0])

            df["Score"] = pd.to_numeric(df["Score"])
            df = df.loc[df['Score'] >= 0.5]
            positions = df.iloc[:, 0].to_list()
            residues = df.iloc[:, 1].to_list()

            predicted = []
            epitope = ''
            
            for index, value in enumerate(positions):
                if (index != len(positions) - 1 and int(value) + 1 == int(positions[index + 1])):
                    epitope = epitope + residues[index]
                elif (index != len(positions) - 1 and int(value) + 1 != int(positions[index + 1])):
                    epitope = epitope + residues[index]
                    if len(epitope) >= 7:
                        predicted.append(epitope)
                        predicted.append(int(value) - len(epitope) + 1)
                        predicted.append(int(value))
                        predicted = [key] + [conserved_sequence] + predicted
                        bepipred2_results.append(predicted)
                    epitope = ''
                    predicted = []
                elif index == len(positions) - 1:
                    epitope = epitope + residues[index]
                    if len(epitope) >= 7:
                        predicted.append(epitope)
                        predicted.append(int(value) - len(epitope) + 1)
                        predicted.append(int(value))
                        predicted = [key] + [conserved_sequence] + predicted
                        bepipred2_results.append(predicted)
                    epitope = ''
                    predicted = []
    with open('epitope_prediction_results/bepipred2.0_epitopes.csv', 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(bepipred2_results)   
    print("Bepipred 2.0 prediction done\n")

def bepipred(conserved_sequences_dict):
    
    """This function uses the REST API from IEDB to access to Linear Antigen Prediction method Bepipred and predict the linear epitopes of the given conserved 
    sequences using a combination of a hidden Markov model and a propensity scale method. The results are returned as a list of lists"""

    bepipred_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["predicted_epitope"] + ["start_position"] + ["end_position"]
    bepipred_results.append(columns)

    for key in conserved_sequences_dict:
        print("Predicting linear epitopes of protein " + key + " with Bepipred 1.0")
        for conserved_sequence in conserved_sequences_dict[key]:
            
            time.sleep(1)
            headers = {
                'Content-Type': 'application/x-www-form-urlencoded',
                }
            
            data = {
                'method': 'Bepipred',
                'sequence_text': conserved_sequence,
                'window_size': '7',
                }
            
            response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/bcell/', headers=headers, data=data)
            response_data_split_by_line = response.content.decode('utf-8').splitlines()
            
            response_body = []
            for line in response_data_split_by_line:
                split_line = line.split("\t")
                response_body.append(split_line)
            
            df = pd.DataFrame(response_body[1:], columns=response_body[0])

            df["Score"] = pd.to_numeric(df["Score"])
            df = df.loc[df['Score'] >= 0.35]
            positions = df.iloc[:, 0].to_list()
            residues = df.iloc[:, 1].to_list()

            predicted = []
            epitope = ''
            
            for index, value in enumerate(positions):
                if (index != len(positions) - 1 and int(value) + 1 == int(positions[index + 1])):
                    epitope = epitope + residues[index]
                elif (index != len(positions) - 1 and int(value) + 1 != int(positions[index + 1])):
                    epitope = epitope + residues[index]
                    if len(epitope) >= 7:
                        predicted.append(epitope)
                        predicted.append(int(value) - len(epitope) + 1)
                        predicted.append(int(value))
                        predicted = [key] + [conserved_sequence] + predicted
                        bepipred_results.append(predicted)
                    epitope = ''
                    predicted = []
                elif index == len(positions) - 1:
                    epitope = epitope + residues[index]
                    if len(epitope) >= 7:
                        predicted.append(epitope)
                        predicted.append(int(value) - len(epitope) + 1)
                        predicted.append(int(value))
                        predicted = [key] + [conserved_sequence] + predicted
                        bepipred_results.append(predicted)
                    epitope = ''
                    predicted = []
    with open('epitope_prediction_results/bepipred1.0_epitopes.csv', 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(bepipred_results) 
    print("Bepipred 1.0 prediction done\n")

def emini(conserved_sequences_dict):
    
    """This function uses the REST API from IEDB to access the Linear Antigen Prediction method Emini and predict the linear epitopes of the given conserved 
    sequences. The calculation is based on surface accessibility scale on a product. The results are returned as a list of lists"""

    emini_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["predicted_epitope"] + ["start_position"] + ["end_position"]
    emini_results.append(columns)

    for key in conserved_sequences_dict:
        print("Predicting linear epitopes of protein " + key + " with Emini")
        for conserved_sequence in conserved_sequences_dict[key]:
            
            time.sleep(1)
            headers = {
                'Content-Type': 'application/x-www-form-urlencoded',
                }
            
            data = {
                'method': 'Emini',
                'sequence_text': conserved_sequence,
                'window_size': '6',
                }
            
            response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/bcell/', headers=headers, data=data)
            response_data_split_by_line = response.content.decode('utf-8').splitlines()
            
            response_body = []
            for line in response_data_split_by_line:
                split_line = line.split("\t")
                response_body.append(split_line)

            df = pd.DataFrame(response_body[1:], columns=response_body[0])

            df["Score"] = pd.to_numeric(df["Score"])
            df = df.loc[df['Score'] >= 1]
            positions = df.iloc[:, 0].to_list()
            residues = df.iloc[:, 1].to_list()

            predicted = []
            epitope = ''
            
            for index, value in enumerate(positions):
                if (index != len(positions) - 1 and int(value) + 1 == int(positions[index + 1])):
                    epitope = epitope + residues[index]
                elif (index != len(positions) - 1 and int(value) + 1 != int(positions[index + 1])):
                    epitope = epitope + residues[index]
                    if len(epitope) >= 7:
                        predicted.append(epitope)
                        predicted.append(int(value) - len(epitope) + 1)
                        predicted.append(int(value))
                        predicted = [key] + [conserved_sequence] + predicted
                        emini_results.append(predicted)
                    epitope = ''
                    predicted = []
                elif index == len(positions) - 1:
                    epitope = epitope + residues[index]
                    if len(epitope) >= 7:
                        predicted.append(epitope)
                        predicted.append(int(value) - len(epitope) + 1)
                        predicted.append(int(value))
                        predicted = [key] + [conserved_sequence] + predicted
                        emini_results.append(predicted)
                    epitope = ''
                    predicted = []
    with open('epitope_prediction_results/emini_epitopes.csv', 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(emini_results) 
    print("Emini prediction done\n")

def choufasman(conserved_sequences_dict):

    """This function uses the REST API from IEDB to access the Linear Antigen Prediction method Chou-Fasman and predict the linear epitopes of the given 
    conserved sequences. It uses the Chou and Fasman scale which is commonly used to predict beta turns. The results are returned as a list of lists"""

    choufasman_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["predicted_epitope"] + ["start_position"] + ["end_position"]
    choufasman_results.append(columns)

    for key in conserved_sequences_dict:
        print("Predicting linear epitopes of protein " + key + " with Chou-Fasman")
        for conserved_sequence in conserved_sequences_dict[key]:
            
            time.sleep(1)
            headers = {
                'Content-Type': 'application/x-www-form-urlencoded',
                }
            
            data = {
                'method': 'Chou-Fasman',
                'sequence_text': conserved_sequence,
                'window_size': '7',
                }
            
            response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/bcell/', headers=headers, data=data)
            response_data_split_by_line = response.content.decode('utf-8').splitlines()
            
            response_body = []
            for line in response_data_split_by_line:
                split_line = line.split("\t")
                response_body.append(split_line)
            
            df = pd.DataFrame(response_body[1:], columns=response_body[0])
            df["Score"] = pd.to_numeric(df["Score"])
            average_of_scores = df["Score"].mean()
            df = df.loc[df['Score'] >= average_of_scores]

            positions = df.iloc[:, 0].to_list()
            residues = df.iloc[:, 1].to_list()

            predicted = []
            epitope = ''
            
            for index, value in enumerate(positions):
                if (index != len(positions) - 1 and int(value) + 1 == int(positions[index + 1])):
                    epitope = epitope + residues[index]
                elif (index != len(positions) - 1 and int(value) + 1 != int(positions[index + 1])):
                    epitope = epitope + residues[index]
                    if len(epitope) >= 7:
                        predicted.append(epitope)
                        predicted.append(int(value) - len(epitope) + 1)
                        predicted.append(int(value))
                        predicted = [key] + [conserved_sequence] + predicted
                        choufasman_results.append(predicted)
                    epitope = ''
                    predicted = []
                elif index == len(positions) - 1:
                    epitope = epitope + residues[index]
                    if len(epitope) >= 7:
                        predicted.append(epitope)
                        predicted.append(int(value) - len(epitope) + 1)
                        predicted.append(int(value))
                        predicted = [key] + [conserved_sequence] + predicted
                        choufasman_results.append(predicted)
                    epitope = ''
                    predicted = []
    with open('epitope_prediction_results/choufasman_epitopes.csv', 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(choufasman_results) 
    print("Chou-Fasman prediction done\n")

def karplusschulz(conserved_sequences_dict):

    """This function uses the REST API from IEDB to access the Linear Antigen Prediction method Karplus-Schulz and predict the linear epitopes of the given 
    conserved sequences. In this method, flexibility scale based on mobility of protein segments on the basis of the known temperature B factors of the 
    a-carbons of 31 proteins of known structure was constructed. The results are returned as a list of lists"""

    karplusschulz_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["predicted_epitope"] + ["start_position"] + ["end_position"]
    karplusschulz_results.append(columns)

    for key in conserved_sequences_dict:
        print("Predicting linear epitopes of protein " + key + " with Karplus-Schulz")
        for conserved_sequence in conserved_sequences_dict[key]:
            
            time.sleep(1)
            headers = {
                   'User-agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/47.0.2526.80 Safari/537.36',
                   'Content-Type': 'application/x-www-form-urlencoded'
                }
            
            data = {
                'method': 'Karplus-Schulz',
                'sequence_text': conserved_sequence,
                'window_size': '7',
                }
            
            response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/bcell/', headers=headers, data=data)
            response_data_split_by_line = response.content.decode('utf-8').splitlines()
            
            response_body = []
            for line in response_data_split_by_line:
                split_line = line.split("\t")
                response_body.append(split_line)

            df = pd.DataFrame(response_body[1:], columns=response_body[0])
            df["Score"] = pd.to_numeric(df["Score"])
            average_of_scores = df["Score"].mean()
            df = df.loc[df['Score'] >= average_of_scores]

            positions = df.iloc[:, 0].to_list()
            residues = df.iloc[:, 1].to_list()

            predicted = []
            epitope = ''
            
            for index, value in enumerate(positions):
                if (index != len(positions) - 1 and int(value) + 1 == int(positions[index + 1])):
                    epitope = epitope + residues[index]
                elif (index != len(positions) - 1 and int(value) + 1 != int(positions[index + 1])):
                    epitope = epitope + residues[index]
                    if len(epitope) >= 7:
                        predicted.append(epitope)
                        predicted.append(int(value) - len(epitope) + 1)
                        predicted.append(int(value))
                        predicted = [key] + [conserved_sequence] + predicted
                        karplusschulz_results.append(predicted)
                    epitope = ''
                    predicted = []
                elif index == len(positions) - 1:
                    epitope = epitope + residues[index]
                    if len(epitope) >= 7:
                        predicted.append(epitope)
                        predicted.append(int(value) - len(epitope) + 1)
                        predicted.append(int(value))
                        predicted = [key] + [conserved_sequence] + predicted
                        karplusschulz_results.append(predicted)
                    epitope = ''
                    predicted = []
    with open('epitope_prediction_results/karplusschulz_epitopes.csv', 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(karplusschulz_results) 
    print("Karplus-Schulz prediction done\n")

def kolaskartongaonkar(conserved_sequences_dict):

    """This function uses the REST API from IEDB to access the Linear Antigen Prediction method Kolaskar-Tongaonkar and predict the linear epitopes of the 
    given conserved sequences. It is a semi-empirical method which makes use of physicochemical properties of amino acid residues and their frequencies of 
    occurrence in experimentally known segmental epitopes was developed to predict antigenic determinants on proteins. The results are returned as a list of lists"""

    kolaskartongaonkar_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["predicted_epitope"] + ["start_position"] + ["end_position"]
    kolaskartongaonkar_results.append(columns)

    for key in conserved_sequences_dict:
        print("Predicting linear epitopes of protein " + key + " with Kolaskar-Tongaonkar")
        for conserved_sequence in conserved_sequences_dict[key]:
            
            time.sleep(1)
            headers = {
                'Content-Type': 'application/x-www-form-urlencoded',
                }
            
            data = {
                'method': 'Kolaskar-Tongaonkar',
                'sequence_text': conserved_sequence,
                'window_size': '7',
                }
            
            response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/bcell/', headers=headers, data=data)
            response_data_split_by_line = response.content.decode('utf-8').splitlines()

            response_body = []
            for line in response_data_split_by_line:
                split_line = line.split("\t")
                response_body.append(split_line)

            df = pd.DataFrame(response_body[1:], columns=response_body[0])
            df["Score"] = pd.to_numeric(df["Score"])
            average_of_scores = df["Score"].mean()
            df = df.loc[df['Score'] >= average_of_scores]

            positions = df.iloc[:, 0].to_list()
            residues = df.iloc[:, 1].to_list()

            predicted = []
            epitope = ''
            
            for index, value in enumerate(positions):
                if (index != len(positions) - 1 and int(value) + 1 == int(positions[index + 1])):
                    epitope = epitope + residues[index]
                elif (index != len(positions) - 1 and int(value) + 1 != int(positions[index + 1])):
                    epitope = epitope + residues[index]
                    if len(epitope) >= 7:
                        predicted.append(epitope)
                        predicted.append(int(value) - len(epitope) + 1)
                        predicted.append(int(value))
                        predicted = [key] + [conserved_sequence] + predicted
                        kolaskartongaonkar_results.append(predicted)
                    epitope = ''
                    predicted = []
                elif index == len(positions) - 1:
                    epitope = epitope + residues[index]
                    if len(epitope) >= 7:
                        predicted.append(epitope)
                        predicted.append(int(value) - len(epitope) + 1)
                        predicted.append(int(value))
                        predicted = [key] + [conserved_sequence] + predicted
                        kolaskartongaonkar_results.append(predicted)
                    epitope = ''
                    predicted = []
    with open('epitope_prediction_results/kolaskartongaonkar_epitopes.csv', 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(kolaskartongaonkar_results) 
    print("Kolaskar-Tongaonkar prediction done\n")

def parker(conserved_sequences_dict):

    """This function uses the REST API from IEDB to access the Linear Antigen Prediction method Parker and predict the linear epitopes of the given 
    conserved sequences. In this method, hydrophilic scale based on peptide retention times during high-performance liquid chromatography (HPLC) on a 
    reversed-phase column was constructed. A window of seven residues was used for analyzing epitope region. The results are returned as a list of lists"""

    parker_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["predicted_epitope"] + ["start_position"] + ["end_position"]
    parker_results.append(columns)

    for key in conserved_sequences_dict:
        print("Predicting linear epitopes of protein " + key + " with Parker")
        for conserved_sequence in conserved_sequences_dict[key]:
            
            time.sleep(1)
            headers = {
                'Content-Type': 'application/x-www-form-urlencoded',
                }
            
            data = {
                'method': 'Parker',
                'sequence_text': conserved_sequence,
                'window_size': '7',
                }
            
            response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/bcell/', headers=headers, data=data)
            response_data_split_by_line = response.content.decode('utf-8').splitlines()
            
            response_body = []
            for line in response_data_split_by_line:
                split_line = line.split("\t")
                response_body.append(split_line)

            df = pd.DataFrame(response_body[1:], columns=response_body[0])
            df["Score"] = pd.to_numeric(df["Score"])
            average_of_scores = df["Score"].mean()
            df = df.loc[df['Score'] >= average_of_scores]

            positions = df.iloc[:, 0].to_list()
            residues = df.iloc[:, 1].to_list()

            predicted = []
            epitope = ''
            
            for index, value in enumerate(positions):
                if (index != len(positions) - 1 and int(value) + 1 == int(positions[index + 1])):
                    epitope = epitope + residues[index]
                elif (index != len(positions) - 1 and int(value) + 1 != int(positions[index + 1])):
                    epitope = epitope + residues[index]
                    if len(epitope) >= 7:
                        predicted.append(epitope)
                        predicted.append(int(value) - len(epitope) + 1)
                        predicted.append(int(value))
                        predicted = [key] + [conserved_sequence] + predicted
                        parker_results.append(predicted)
                    epitope = ''
                    predicted = []
                elif index == len(positions) - 1:
                    epitope = epitope + residues[index]
                    if len(epitope) >= 7:
                        predicted.append(epitope)
                        predicted.append(int(value) - len(epitope) + 1)
                        predicted.append(int(value))
                        predicted = [key] + [conserved_sequence] + predicted
                        parker_results.append(predicted)
                    epitope = ''
                    predicted = []
    with open('epitope_prediction_results/parker_epitopes.csv', 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(parker_results) 
    print("Parker prediction done\n")

def ellipro(list_of_pdb_ids):
    
    """This function uses Selenium to access the Ellipro tool from IEDB. It requires a list of PDB IDs and returns two lists of lists, one with predicted linear 
    sequences, the other with predicted discontinous sequences. It implements a previously developed method that represents the protein structure as an 
    ellipsoid and calculates protrusion indexes for protein residues outside of the ellipsoid."""
    
    options = Options()
    options.headless = True

    ellipro_url = 'http://tools.iedb.org/ellipro/'
    linear_columns = ['pdb_id','chain','start','end','peptide','nr_of_residues','score']
    discontinous_columns = ['pdb_id','chain','peptide','start','end','nr_of_residues','score']
    linear_epitopes = [linear_columns]
    discontinous_epitopes = [discontinous_columns]

    for protein_id in list_of_pdb_ids:
        print("Predicting linear and discontinous epitopes of protein " + protein_id + " with Ellipro")
        driver = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
        driver.get(ellipro_url)
        driver.find_element(By.NAME, "pdb_id").send_keys(protein_id)
        driver.find_element(By.NAME, "submit").click()

        wait = WebDriverWait(driver, 180)
        wait.until(ec.visibility_of_element_located((By.ID, "result_table")))
        driver.find_element(By.NAME, "chain").click()
        chain = driver.find_element(By.XPATH, "/html/body/div[3]/form/table/tbody/tr[1]/td[3]").text
        driver.find_element(By.NAME, "submit").click()

        wait.until(ec.visibility_of_element_located((By.CLASS_NAME, "output_title")))
        table_of_linear_epitopes = driver.find_element(By.XPATH, '/html/body/div[3]/table[2]/tbody').text.splitlines()
        columns_of_linear_epitopes = driver.find_element(By.XPATH, '/html/body/div[3]/table[2]/thead').text.splitlines()

        table_of_discontinous_epitopes = driver.find_element(By.XPATH, '/html/body/div[3]/table[3]/tbody').text.splitlines()
        columns_of_discontinous_epitopes = driver.find_element(By.XPATH, '/html/body/div[3]/table[3]/thead').text.splitlines()

        linear_epitopes_split_rows = []
        for row in table_of_linear_epitopes:
            split_row = row.split(" ")
            linear_epitopes_split_rows.append(split_row)
        
        df = pd.DataFrame(linear_epitopes_split_rows, columns=columns_of_linear_epitopes[:-1])
        df["Score"] = pd.to_numeric(df["Score"])
        df = df.loc[df['Score'] >= 0.7]

        if df.empty == True:
            continue
        else:
            rows = [[i for i in row[1:]] for row in df.itertuples()]
            for i in rows:
                i = [protein_id] + i[1:]
                linear_epitopes.append(i)

        rows = []
        for i in range(len(table_of_discontinous_epitopes)):
            row = []
            for e in range(len(columns_of_discontinous_epitopes[:-1])):
                cell = driver.find_element(By.XPATH, '/html/body/div[3]/table[3]/tbody/tr[' + str(i+1) + ']/td[' + str(e+1) + ']').text
                row.append(cell)
            rows.append(row)
        driver.close()

        for arow in rows:
            if float(arow[3]) >= 0.7:
                unproc_aa_list = arow[1].split(', ')
                aa_list = []
                for position in unproc_aa_list:
                    aa = position[2:]
                    aa_list.append(aa)
                epitope = ','.join(aa_list)
                start_pos = unproc_aa_list[0][3:]
                end_pos = unproc_aa_list[-1][3:]
                nr_of_residues = len(aa_list)
                score = float(arow[3])
                row_to_append = [protein_id] + [chain] + [epitope] + [start_pos] + [end_pos] + [nr_of_residues] + [score]
                discontinous_epitopes.append(row_to_append)
    with open('epitope_prediction_results/ellipro_linear_epitopes.csv', 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(linear_epitopes) 
    with open('epitope_prediction_results/ellipro_discontinous_epitopes.csv', 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(discontinous_epitopes) 
    print("Ellipro prediction done\n")

def discotope(list_of_pdb_ids):
    
    """This function uses Selenium to access the Discotope 2.0 tool from IEDB. It requires a list of PDB IDs and returns 
    ..."""
    
    options = Options()
    options.headless = True

    discotope_url = 'http://tools.iedb.org/discotope/'
    discotope_columns = ['pdb_id','chain','peptide','start','end','nr_of_residues']
    discotope_epitopes = [discotope_columns]

    for protein_id in list_of_pdb_ids:
        print("Predicting discontinous epitopes of protein " + protein_id + " with Discotope")
        driver = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
        driver.get(discotope_url)
        driver.find_element(By.ID, "id_pdb").send_keys(protein_id)
        chain = 'A'
        driver.find_element(By.ID, "id_chain").send_keys(chain)
        driver.find_element(By.XPATH, "/html/body/div[3]/form/table/tbody/tr[3]/td[2]/select/option[2]").click()
        driver.find_element(By.NAME, "submit").click()

        wait = WebDriverWait(driver, 600)
        wait.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div[3]/form")))
        driver.find_element(By.XPATH, "/html/body/div[3]/form/a[1]/button").click()
        wait.until(ec.visibility_of_element_located((By.ID, "result_table")))
        discotope_table = driver.find_element(By.XPATH, '/html/body/div[3]/table/tbody').text.splitlines()
        discotope_table_columns = driver.find_element(By.XPATH, '/html/body/div[3]/table/thead').text.splitlines()
        driver.close()

        discotope_epitopes_split_rows = []
        for row in discotope_table:
            split_row = row.split(" ")
            discotope_epitopes_split_rows.append(split_row)
        
        df = pd.DataFrame(discotope_epitopes_split_rows, columns=discotope_table_columns)
        df["Discotope Score"] = pd.to_numeric(df["Discotope Score"])
        df = df.loc[df['Discotope Score'] >= -3.7]

        if df.empty == True:
            continue
        else:
            residue_id = df.iloc[:, 1].to_list()
            residue_name = df.iloc[:, 2].to_list()

        discontinous_epitope = []
        for i in range(len(residue_name)):
            if residue_name[i] == 'ARG':
                ARG = 'R'
                pos = ARG + str(residue_id[i])
                discontinous_epitope.append(pos)
            elif residue_name[i] == 'ASN':
                ASN = 'N'
                pos = ASN + str(residue_id[i])
                discontinous_epitope.append(pos)
            elif residue_name[i] == 'ASP':
                ASP = 'D'
                pos = ASP + str(residue_id[i])
                discontinous_epitope.append(pos)
            elif residue_name[i] == 'GLN':
                GLN = 'Q'
                pos = GLN + str(residue_id[i])
                discontinous_epitope.append(pos)
            elif residue_name[i] == 'GLU':
                GLU = 'E'
                pos = GLU + str(residue_id[i])
                discontinous_epitope.append(pos)
            elif residue_name[i] == 'LYS':
                LYS = 'K'
                pos = LYS + str(residue_id[i])
                discontinous_epitope.append(pos)
            elif residue_name[i] == 'PHE':
                PHE = 'F'
                pos = PHE + str(residue_id[i])
                discontinous_epitope.append(pos)
            elif residue_name[i] == 'TRP':
                TRP = 'W'
                pos = TRP + str(residue_id[i])
                discontinous_epitope.append(pos)
            elif residue_name[i] == 'TYR':
                TYR = 'Y'
                pos = TYR + str(residue_id[i])
                discontinous_epitope.append(pos)
            else:
                pos = residue_name[i][0] + str(residue_id[i])
                discontinous_epitope.append(pos)
        
        epitope = ','.join(discontinous_epitope)
        start_pos = residue_id[0]
        end_pos = residue_id[-1]
        row_to_append = [protein_id] + [chain] + [epitope] + [start_pos] + [end_pos] + [len(residue_name)]
        discotope_epitopes.append(row_to_append)
    with open('epitope_prediction_results/discotope_epitopes.csv', 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(discotope_epitopes)  
    print("Discotope prediction done\n")


def predict_all(dictionary_conserved_sequences, alleles_for_mhci, lengths_for_mhci, alleles_for_mhcii, lengths_for_mhcii, list_of_pdb_ids):
    
    """This function runs all the prediction methods above and returns a tuple with the lists of lists. It also \
        creates a folder where all the results are stored as csv files"""
    
    current_directory = os.getcwd()
    final_directory = os.path.join(current_directory, r'epitope_prediction_results')
    if not os.path.exists(final_directory):
        os.makedirs(final_directory)

    mhci(dictionary_conserved_sequences, alleles_for_mhci, lengths_for_mhci)
    mhci_proc(dictionary_conserved_sequences, alleles_for_mhci, lengths_for_mhci)
    mhcii(dictionary_conserved_sequences, alleles_for_mhcii, lengths_for_mhcii)
    
    bepipred2(dictionary_conserved_sequences)
    bepipred(dictionary_conserved_sequences)
    emini(dictionary_conserved_sequences)
    choufasman(dictionary_conserved_sequences)
    karplusschulz(dictionary_conserved_sequences)
    kolaskartongaonkar(dictionary_conserved_sequences)
    parker(dictionary_conserved_sequences)

    ellipro(list_of_pdb_ids)
    discotope(list_of_pdb_ids)



# class Netctlpan(Spider):
#     name = 'Netctlpan'
#     start_urls = ['http://tools.iedb.org/netchop/']
#     global netctlpan
#     global fasta

#     def parse(self, response):
#         yield FormRequest.from_response(
#             response,
#             url=self.start_urls[0],
#             formdata={
#                 'pred_tool': 'netchop',
#                 'pred_method': 'netctlpan',
#                 'sequence_text': fasta,
#                 'sequence_file': '(binary)',
#                 'method': '0',
#                 'netchop_threshold': '0.5',
#                 'netctl_cleavage': '0.15',
#                 'netctl_tap': '0.05',
#                 'supertype': 'A1',
#                 'netctl_threshold': '0.75',
#                 'species_list': 'human',
#                 'freq': 'freq',
#                 'allele_list': 'HLA-A01:01',
#                 'length_list': '9',
#                 'netctlpan_threshold': '-99.9',
#                 'netctlpan_cleavage': '0.225',
#                 'netctlpan_tap': '0.025',
#                 'epitope_threshold': '1.0'
#             },
#             callback=self.results_page)

#     def results_page(self, response):
#         url = 'http://tools.iedb.org/netchop/table/'
#         yield scrapy.Request(url=url, callback=self.get_results, dont_filter=True)

#     def get_results(self, results):
#         row = []
#         for cell in results.xpath('/html/body/div[3]/table/tbody/tr'):
#             row.append(cell.xpath('td[1]//text()').extract_first())
#             row.append(cell.xpath('td[2]//text()').extract_first())
#             row.append(cell.xpath('td[3]//text()').extract_first())
#             row.append(cell.xpath('td[4]//text()').extract_first())
#             row.append(cell.xpath('td[5]//text()').extract_first())
#             row.append(cell.xpath('td[6]//text()').extract_first())
#             netctlpan.append(row)
#             row = []
#         print(netctlpan)

#configure_logging()
#settings = get_project_settings()
#runner = CrawlerRunner(settings)

# @defer.inlineCallbacks
# def crawl():
    # yield runner.crawl(Ellipro)
    # yield runner.crawl(Discotope)
    # yield runner.crawl(Bepipred2)
    # yield runner.crawl(Bepipred)
    # yield runner.crawl(Emini)
    # yield runner.crawl(Kolaskar)
    # yield runner.crawl(Chou_Fasman)
    # yield runner.crawl(Karplus_Schulz)
    # yield runner.crawl(Parker)
    # yield runner.crawl(MhcII)
    # yield runner.crawl(Netctlpan)
    # reactor.stop()

# crawl()
# reactor.run()  # the script will block here until the last crawl call is finished


