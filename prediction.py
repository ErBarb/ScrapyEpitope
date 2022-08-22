from tkinter import E
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as ec
from selenium.webdriver.firefox.options import Options
import time
import os
import pandas as pd
import requests
import json
import csv
import re
from collections import defaultdict
import itertools
import matplotlib.pyplot as plt


def get_pdb_from_swissprot(list_of_swissprot_ids):

    """This function gets a PDB ID for each of the Swissprot IDs by looking for EM and X-Ray data on the Swissprot page of the ID and gets only
    PDB IDs that start with 1 or with 10 to 19. Ex: 1-2000 or 13-2000. An ID with 21-2000 is not accepted. The function outputs a list of the PDB
    IDs."""

    options = Options()
    options.headless = True

    pdb_ids = []
    for id in list_of_swissprot_ids:
        try:
            uniprot_url = 'https://www.uniprot.org/uniprot/' + id
            uniprot = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
            uniprot.get(uniprot_url)
            wait_uniprot = WebDriverWait(uniprot, 60)
            wait_uniprot.until(ec.visibility_of_element_located((By.XPATH, '//section[@id="structure"]/div/div[2]/protvista-uniprot-structure/div/div[2]/protvista-datatable/table/tbody')))
            uniprot_pdb_ids = uniprot.find_element(By.XPATH, '//section[@id="structure"]/div/div[2]/protvista-uniprot-structure/div/div[2]/protvista-datatable/table/tbody').text.splitlines()
            wait_uniprot.until(ec.visibility_of_element_located((By.XPATH, '//section[@id="sequence-container"]/ul/li[1]/div/div[2]')))
            for row in uniprot_pdb_ids:
                newrow = row.split()
                if (newrow[2] == 'EM' or newrow[2] == 'X-ray') and (newrow[6].startswith('1-') == True):
                    pdb_ids.append(newrow[1])
                    print("The following PDB ID was found for " + id + ": " + newrow[1])
                    break
            else:
                if (newrow[2] == 'EM' or newrow[2] == 'X-ray') and re.search("^1.-", newrow[6]):
                    print("The following PDB ID was found for " + id + ": " + newrow[1])
                    pdb_ids.append(newrow[1])
            uniprot.close()
        except:
            print("No PDB format for: " + id)
            uniprot.close()
    
    return pdb_ids

def getStartEndPositions(pattern, seq):
    all_instances = [m.start() for m in re.finditer(pattern, seq)]

    if len(all_instances) > 1:
        raise Exception("Too many found")
    elif len(all_instances) == 0:
        return None

    return (all_instances[0], all_instances[0] + len(pattern))

def swissprotIDSequenceLength(list_of_swissprot_ids):

    """This function outputs a dictionary with the sequences and length of the sequence (values) of Swissprot ID (key)"""

    protein_dict = {}
    for id in list_of_swissprot_ids:
        baseUrl="http://www.uniprot.org/uniprot/"
        currentUrl=baseUrl+id+".fasta"
        response = requests.post(currentUrl)
        response_lines = response.text.split("\n")[1:]
        cData=''.join(response_lines)
        protein_dict[id] = (cData, len(cData))
    return protein_dict

def mhci(conserved_sequences_dict, list_of_alleles, list_of_lengths):
    
    """This function uses the REST API from IEDB to access to MHC I Binding tool and predict the MHC Class I epitopes of the given conserved sequences. 
    It also needs a list of alleles and their respective lengths to run. It employs different methods to predict MHC Class I epitopes, including a 
    consensus approach which combines ANN, SMM, Comblib, NetMHCpan and Consensus. The results are returned as a list of lists"""

    #alleles = ",".join(list_of_alleles)
    converted_list = [str(element) for element in list_of_lengths]
    lengths = ",".join(converted_list)

    mhci_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["allele"] + ["seq_num"] + ["start"] + ["end"] + ["length"] + \
        ["peptide"] + ["core"] + ["icore"] + ["score"] + ["percentile_rank"]
    mhci_results.append(columns)
    #print(conserved_sequences_dict)
    for key, value in conserved_sequences_dict.items():
        print("Predicting linear epitopes of protein " + key + " with MHCI")
        for conserved_sequence in value:
            if len(conserved_sequence) <= 4000:
            
                for allele in list_of_alleles:
                    time.sleep(10)
                    print(key, allele)
                    input_allele = []
                    for i in range(len(list_of_lengths)):
                        input_allele.append(allele)
                    alleles = ",".join(input_allele)     

                    if len(conserved_sequence) >= max(list_of_lengths):
                        
                        try:

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

                            df = pd.DataFrame(response_body[1:], columns=response_body[0])
                            df["percentile_rank"] = pd.to_numeric(df["percentile_rank"], errors='coerce')
                            df = df.loc[df['percentile_rank'] <= 1]

                            if df.empty == True:
                                #print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence)) + " and allele " + allele)
                                continue
                            else:
                                rows = [[i for i in row[1:]] for row in df.itertuples()]
                                for i in rows:
                                    i = [key] + [conserved_sequence] + i
                                    mhci_results.append(i)

                        except:
                            
                            #print("Retrying prediction of linear epitopes of protein " + key + " and allele " + allele + " with MHCI")

                            try:

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

                                df = pd.DataFrame(response_body[1:], columns=response_body[0])
                                df["percentile_rank"] = pd.to_numeric(df["percentile_rank"], errors='coerce')
                                df = df.loc[df['percentile_rank'] <= 1]

                                if df.empty == True:
                                    #print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence))  + " and allele " + allele)
                                    continue
                                else:
                                    rows = [[i for i in row[1:]] for row in df.itertuples()]
                                    for i in rows:
                                        i = [key] + [conserved_sequence] + i
                                        mhci_results.append(i)
                            
                            except:
                                print("Epitope prediction for protein " + key + " with sequence length " + str(len(conserved_sequence)) + " and allele " + allele + " failed")
                                continue
            else:
                print("Protein too long.")
                continue

    with open('epitope_prediction_results/mhci_epitopes.csv', 'w', newline="") as f: 
        writer = csv.writer(f)
        writer.writerows(mhci_results)
    print("MHCI prediction done\n")

def mhci_proc(conserved_sequences_dict, list_of_alleles, list_of_lengths):

    """This function uses the REST API from IEDB to access to MHC I Processing tool and predict the MHC Class I epitopes of the given conserved sequences. 
    It also needs a list of alleles and their respective lengths to run. It combines predictors of proteasomal processing, TAP transport, and MHC binding 
    to produce an overall score for each peptide's intrinsic potential of being a T cell epitope. The results are returned as a list of lists"""


    #alleles = ",".join(list_of_alleles)
    converted_list = [str(element) for element in list_of_lengths]
    lengths = ",".join(converted_list)

    mhci_proc_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["allele"] + ['#'] + ["start"] + ["end"] + ['peptide_length'] + \
        ["peptide"] + ["proteasome_score"] + ["tap_score"] + ["mhci_score"] + ["processing_score"] + ["total_score"] + ["mhci_ic50"]
    mhci_proc_results.append(columns)

    for key, value in conserved_sequences_dict.items():
        print("Predicting linear epitopes of protein " + key + " with MHCI Processing")
        for conserved_sequence in value:
            if len(conserved_sequence) <= 4000:
                for allele in list_of_alleles:
                    time.sleep(10)
                    print(key, allele)
                    input_allele = []
                    for i in range(len(list_of_lengths)):
                        input_allele.append(allele)
                    alleles = ",".join(input_allele)

                    if len(conserved_sequence) >= max(list_of_lengths):
                        
                        try:
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
                                #print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence)) + " and allele " + allele)
                                continue
                            else:
                                rows = [[i for i in row[1:]] for row in df.itertuples()]
                                for i in rows:
                                    i = [key] + [conserved_sequence] + i
                                    mhci_proc_results.append(i)
                        except:

                            #print("Retrying prediction of linear epitopes of protein " + key + " and allele " + allele + " with MHCI")
                            
                            try:
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
                                    #print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence))  + " and allele " + allele)
                                    continue
                                else:
                                    rows = [[i for i in row[1:]] for row in df.itertuples()]
                                    for i in rows:
                                        i = [key] + [conserved_sequence] + i
                                        mhci_proc_results.append(i)

                            except:
                                print("Epitope prediction for protein " + key + " with sequence length " + str(len(conserved_sequence)) + " and allele " + allele + " failed")
                                continue
            else:
                print("Protein too long.")
                continue


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

    for key, value in conserved_sequences_dict.items():
        print("Predicting linear epitopes of protein " + key + " with MHCII")
        time.sleep(10)
        for conserved_sequence in value:
            if len(conserved_sequence) >= max(list_of_lengths):
                
                try:
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
                        print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence)))
                        continue
                    else:
                        rows = [[i for i in row[1:]] for row in df.itertuples()]
                        for i in rows:
                            i = [key] + [conserved_sequence] + i
                            mhcii_results.append(i)
                
                except:
                    
                    print("Retrying prediction of linear epitopes of protein " + key + " with MHCII")

                    try:
                        
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
                            print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence)))
                            continue
                        else:
                            rows = [[i for i in row[1:]] for row in df.itertuples()]
                            for i in rows:
                                i = [key] + [conserved_sequence] + i
                                mhcii_results.append(i)

                    except:
                        print("Epitope prediction for protein " + key + " and sequence length " + str(len(conserved_sequence)) + " failed")
                        continue


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

    for key, sequences in conserved_sequences_dict.items():
        print("Predicting linear epitopes of protein " + key + " with Bepipred 2.0")
        time.sleep(10)
        for conserved_sequence in sequences:
            try:
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

                if df.empty == True:
                    print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence)))
                    continue
                else:        
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
            
            except:

                print("Retrying prediction of linear epitopes of protein " + key + " with Bepipred 2.0")

                try:

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
                    
                    if df.empty == True:
                        print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence)))
                        continue
                    else:
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
                except:
                    print("Epitope prediction for protein " + key + " and sequence length " + str(len(conserved_sequence)) + " failed")
                    continue

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

    for key, sequences in conserved_sequences_dict.items():
        print("Predicting linear epitopes of protein " + key + " with Bepipred 1.0")
        for conserved_sequence in sequences:
            try:
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
                if df.empty == True:
                    print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence)))
                    continue
                else:
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
            except:

                print("Retrying prediction of linear epitopes of protein " + key + " with Bepipred 1.0")

                try:
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
                    if df.empty == True:
                        print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence)))
                        continue
                    else:
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
                except:
                    print("Epitope prediction for protein " + key + " and sequence length " + str(len(conserved_sequence)) + " failed")
                    continue
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

    for key, sequences in conserved_sequences_dict.items():
        print("Predicting linear epitopes of protein " + key + " with Emini")
        for conserved_sequence in sequences:
            try:
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
                if df.empty == True:
                    print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence)))
                    continue
                else:
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
            except:

                print("Retrying prediction of linear epitopes of protein " + key + " with Emini")

                try:
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
                    if df.empty == True:
                        print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence)))
                        continue
                    else:
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
                except:
                    print("Epitope prediction for protein " + key + " and sequence length " + str(len(conserved_sequence)) + " failed")
                    continue
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

    for key, sequences in conserved_sequences_dict.items():
        print("Predicting linear epitopes of protein " + key + " with Chou-Fasman")
        for conserved_sequence in sequences:
            try:
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
                
                if df.empty == True:
                    print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence)))
                    continue
                else:
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
            except:

                print("Retrying prediction of linear epitopes of protein " + key + " with Chou-Fasman")

                try:
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
                    
                    if df.empty == True:
                        print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence)))
                        continue
                    else:
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
                except:
                    print("Epitope prediction for protein " + key + " and sequence length " + str(len(conserved_sequence)) + " failed")
                    continue

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

    for key, sequences in conserved_sequences_dict.items():
        print("Predicting linear epitopes of protein " + key + " with Karplus-Schulz")
        for conserved_sequence in sequences:
            try:
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
                if df.empty == True:
                    print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence)))
                    continue
                else:
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
            except:

                print("Retrying prediction of linear epitopes of protein " + key + " with Karplus-Schulz")

                try:
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
                    if df.empty == True:
                        print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence)))
                        continue
                    else:
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
                except:
                    print("Epitope prediction for protein " + key + " and sequence length " + str(len(conserved_sequence)) + " failed")
                    continue
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

    for key, sequences in conserved_sequences_dict.items():
        print("Predicting linear epitopes of protein " + key + " with Kolaskar-Tongaonkar")
        for conserved_sequence in sequences:
            try:
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
                
                if df.empty == True:
                    print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence)))
                    continue
                else:
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
            except:

                print("Retrying prediction of linear epitopes of protein " + key + " with Kolaskar-Tongaonkar")

                try:
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
                    
                    if df.empty == True:
                        print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence)))
                        continue
                    else:
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
                except:
                    print("Epitope prediction for protein " + key + " and sequence length " + str(len(conserved_sequence)) + " failed")
                    continue

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

    for key, sequences in conserved_sequences_dict.items():
        print("Predicting linear epitopes of protein " + key + " with Parker")
        for conserved_sequence in sequences:
            try:
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

                if df.empty == True:
                    print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence)))
                    continue
                else:
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
            except:

                print("Retrying prediction of linear epitopes of protein " + key + " with Parker")

                try:
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

                    if df.empty == True:
                        print("No epitopes predicted for protein " + key + " with sequence length " + str(len(conserved_sequence)))
                        continue
                    else:
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
                except:
                    print("Epitope prediction for protein " + key + " and sequence length " + str(len(conserved_sequence)) + " failed")
                    continue
                
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
        try:
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
                print("No epitopes predicted for PDB ID " + protein_id)
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
        
        except:
            print("Epitope prediction for PDB ID " + protein_id + " failed")
            continue
    with open('epitope_prediction_results/ellipro_linear_epitopes.csv', 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(linear_epitopes) 
    with open('epitope_prediction_results/ellipro_discontinous_epitopes.csv', 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(discontinous_epitopes) 
    print("Ellipro prediction done\n")

def discotope(list_of_pdb_ids):
    
    """This function uses Selenium to access the Discotope 2.0 tool from IEDB. It takes a list of PDB IDs, inputs them into the server,
    processes the data and saves it into a csv file. Discotope is used to predict non-linear B-Cell epitopes from PDB IDs."""
    
    options = Options()
    options.headless = True

    discotope_url = 'http://tools.iedb.org/discotope/'
    discotope_columns = ['pdb_id','chain','peptide','start','end','nr_of_residues']
    discotope_epitopes = [discotope_columns]

    for protein_id in list_of_pdb_ids:
        try:
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
                print("No epitopes predicted for PDB ID " + protein_id)
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

        except:
            print("Epitope prediction for PDB ID " + protein_id + " failed")
            continue
    with open('epitope_prediction_results/discotope_epitopes.csv', 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(discotope_epitopes)  
    print("Discotope prediction done\n")


def predict_all(alleles_for_mhci, lengths_for_mhci, alleles_for_mhcii, lengths_for_mhcii, list_of_pdb_ids):
    
    """This function runs all the prediction methods above and returns a tuple with the lists of lists. It also creates 
    a folder where all the results are stored as csv files"""
    
    current_directory = os.getcwd()
    final_directory = os.path.join(current_directory, r'epitope_prediction_results')
    if not os.path.exists(final_directory):
        os.makedirs(final_directory)

    with open('temporary_files/protein_sequences_dict', 'r') as f:
        dictionary_conserved_sequences = json.load(f)
    
    # for key, value in dictionary_conserved_sequences.items():
    #     print(key)
    #     print(len(value[0]))

    mhci(dictionary_conserved_sequences, alleles_for_mhci, lengths_for_mhci)
    mhci_proc(dictionary_conserved_sequences, alleles_for_mhci, lengths_for_mhci)

    # mhcii(dictionary_conserved_sequences, alleles_for_mhcii, lengths_for_mhcii)
    # bepipred2(dictionary_conserved_sequences)
    # bepipred(dictionary_conserved_sequences)
    # emini(dictionary_conserved_sequences)
    # choufasman(dictionary_conserved_sequences)
    # karplusschulz(dictionary_conserved_sequences)
    # kolaskartongaonkar(dictionary_conserved_sequences)
    # parker(dictionary_conserved_sequences)
    

    #ellipro(list_of_pdb_ids)
    #discotope(list_of_pdb_ids)


def get_start_end_in_dict(rows_of_method_protein, method_dict, swissprot_id, protein_seq):

    """This function corrects the error of start and end positions when predicting epitopes of conserved sequences of proteins, since the tools
    output the start and end positions of the epitopes within the conserved sequence and not the whole protein sequence."""

    start_list = []
    end_list = []
    for index, row in rows_of_method_protein.iterrows():
        
        #start_end_cell = []
        try:
            # print(row['peptide'], row['start'], row['end'])
            row['start'] = getStartEndPositions(row['peptide'], protein_seq)[0]
            row['end'] = getStartEndPositions(row['peptide'], protein_seq)[1]
            # start_end_cell.append(row['start'])
            # start_end_cell.append(row['end'])
            start_list.append(row['start'])
            end_list.append(row['end'])
        except:
            try:
                # print(row['predicted_epitope'], row['start_position'], row['end_position'])
                row['start_position'] = getStartEndPositions(row['predicted_epitope'], protein_seq)[0]
                row['end_position'] = getStartEndPositions(row['predicted_epitope'], protein_seq)[1]
                # start_end_cell.append(row['start_position'])
                # start_end_cell.append(row['end_position'])
                start_list.append(row['start_position'])
                end_list.append(row['end_position'])
            except Exception as E:
                print(E)   

    method_dict[swissprot_id] = [start_list, end_list]


def epitope_distribution_plots(protein_dict, list_of_swissprot_ids):

    """This function creates epitope distribution plots based on the prediction data and the choice given by the user."""

    choice = input("""The following options are available for the epitope distribution anaysis:\n
    1. Plot the linear epitopes from a single method against the protein length of a protein  of your choice\n
    2. Plot the linear T-Cell and B-Cell epitopes against the protein length for each protein (ex. 9 proteins x 2 cell types --> 18 plots)\n
    3. Plot all linear epitopes against the protein length for each protein (ex. 9 proteins --> 9 plots)\n
Please enter 1, 2 or 3...\n""")

    def choose_method(counter = 10):
        
        if counter == 0:
            print("Failed after 10 tries")
            return

        method_list = ["mhci", "mhci_proc", "mhcii", "bepipred2", "bepipred1", "emini", "choufasman", "karplusschulz", "kolaskartongaonkar", "parker"]
        print(*method_list, sep = ", ")
        method_choice = input("""Please write of the following methods to be analysed for epitope distribution:\n""")
        if method_choice in method_list:
            return method_choice
        else:
            print("Please enter one of the methods mentioned...\n")
            return choose_method(counter - 1)
            

    def choose_protein(list_of_swissprot_ids, counter = 10):
        
        if counter == 0:
            print("Failed after 10 tries")
            return

        print(*list_of_swissprot_ids, sep = ", ")
        protein_choice = input("Please enter the swissprot id to be analysed for epitope distribution:\n")
        if protein_choice in list_of_swissprot_ids:
            return protein_choice
        else:
            print("Please enter one of the proteins mentioned...\n")
            return choose_protein(list_of_swissprot_ids, counter - 1)

    epitopes_dict = {}

    mhci_df = pd.read_csv('epitope_prediction_results/mhci_epitopes.csv')
    mhci_proc_df = pd.read_csv('epitope_prediction_results/mhci_proc_epitopes.csv')
    mhcii_df = pd.read_csv('epitope_prediction_results/mhcii_epitopes.csv')
    bepipred2_df = pd.read_csv('epitope_prediction_results/bepipred2.0_epitopes.csv')
    bepipred1_df = pd.read_csv('epitope_prediction_results/bepipred1.0_epitopes.csv')
    emini_df = pd.read_csv('epitope_prediction_results/emini_epitopes.csv')
    choufasman_df = pd.read_csv('epitope_prediction_results/choufasman_epitopes.csv')
    karplusschulz_df = pd.read_csv('epitope_prediction_results/karplusschulz_epitopes.csv')
    kolaskartongaonkar_df = pd.read_csv('epitope_prediction_results/kolaskartongaonkar_epitopes.csv')
    parker_df = pd.read_csv('epitope_prediction_results/parker_epitopes.csv')

    # ellipro_linear_df = pd.read_csv('epitope_prediction_results/ellipro_linear_epitopes.csv')
    # ellipro_discontinous_df = pd.read_csv('epitope_prediction_results/ellipro_discontinous_epitopes.csv')
    # discotope_df = pd.read_csv('epitope_prediction_results/discotope_epitopes.csv')
    # df.drop_duplicates(subset=['peptide'], keep=False, inplace=True)

    mhci_dict = {}
    mhci_proc_dict = {}
    mhcii_dict = {}
    bepipred2_dict = {}
    bepipred1_dict = {}
    emini_dict = {}
    choufasman_dict = {}
    karplusschulz_dict = {}
    kolaskartongaonkar_dict = {}
    parker_dict = {}
    # ellipro_linear_dict = {}
    # ellipro_discontinous_dict = {}
    # discotope_dict = {}



    for key, value in protein_dict.items():

        rows_of_mhci_protein = mhci_df.loc[mhci_df['protein_id'] == key]
        rows_of_mhci_proc_protein = mhci_proc_df.loc[mhci_proc_df['protein_id'] == key]
        rows_of_mhcii_protein = mhcii_df.loc[mhcii_df['protein_id'] == key]
        rows_of_bepipred2_protein = bepipred2_df.loc[bepipred2_df['protein_id'] == key]
        rows_of_bepipred1_protein = bepipred1_df.loc[bepipred1_df['protein_id'] == key]
        rows_of_emini_protein = emini_df.loc[emini_df['protein_id'] == key]
        rows_of_choufasman_protein = choufasman_df.loc[choufasman_df['protein_id'] == key]
        rows_of_karplusschulz_protein = karplusschulz_df.loc[karplusschulz_df['protein_id'] == key]
        rows_of_kolaskartongaonkar_protein = kolaskartongaonkar_df.loc[kolaskartongaonkar_df['protein_id'] == key]
        rows_of_parker_protein = parker_df.loc[parker_df['protein_id'] == key]

        # rows_of_ellipro_linear_protein = ellipro_linear_df.loc[ellipro_linear_df['protein_id'] == key]
        # rows_of_ellipro_discontinous_protein = ellipro_discontinous_df.loc[ellipro_discontinous_df['protein_id'] == key]
        # rows_of_discotope_protein = discotope_df.loc[discotope_df['protein_id'] == key]
        
        get_start_end_in_dict(rows_of_mhci_protein, mhci_dict, key, value[0])
        get_start_end_in_dict(rows_of_mhci_proc_protein, mhci_proc_dict, key, value[0])
        get_start_end_in_dict(rows_of_mhcii_protein, mhcii_dict, key, value[0])
        get_start_end_in_dict(rows_of_bepipred2_protein, bepipred2_dict, key, value[0])
        get_start_end_in_dict(rows_of_bepipred1_protein, bepipred1_dict, key, value[0])
        get_start_end_in_dict(rows_of_emini_protein, emini_dict, key, value[0])
        get_start_end_in_dict(rows_of_choufasman_protein, choufasman_dict, key, value[0])
        get_start_end_in_dict(rows_of_karplusschulz_protein, karplusschulz_dict, key, value[0])
        get_start_end_in_dict(rows_of_kolaskartongaonkar_protein, kolaskartongaonkar_dict, key, value[0])
        get_start_end_in_dict(rows_of_parker_protein, parker_dict, key, value[0])
    
    epitopes_dict["mhci"] = mhci_dict
    epitopes_dict["mhci_proc"] = mhci_proc_dict
    epitopes_dict["mhcii"] = mhcii_dict
    epitopes_dict["bepipred2"] = bepipred2_dict
    epitopes_dict["bepipred1"] = bepipred1_dict
    epitopes_dict["emini"] = emini_dict
    epitopes_dict["choufasman"] = choufasman_dict
    epitopes_dict["karplusschulz"] = karplusschulz_dict
    epitopes_dict["kolaskartongaonkar"] = kolaskartongaonkar_dict
    epitopes_dict["parker"] = parker_dict
    
    # for k, v in epitopes_dict.items():
    #     print(k, v)

    if choice == '1':

        method = choose_method()
        protein = choose_protein(list_of_swissprot_ids)
        #print(method, protein)

        for swissprot_id, seq_length in protein_dict.items():
            if swissprot_id == protein:
                length = seq_length[1]
                
        
        for key, value in epitopes_dict.items():
            if key == method:
                for k, v in value.items():
                    if k == protein:
                        start_pos = v[0]
                        end_pos = v[1]
                        count = [0]*length
                        for i in range(0,len(start_pos)):
                            for j in range(start_pos[i],end_pos[i]):
                                count[j] = count[j]+1
                        
                        x_axis = list(range(1,length+1))

                        plt.plot(x_axis,count)
                        plt.title(method + ' ' + protein)
                        plt.xlabel('Spike protein (' + protein + ') residues')
                        plt.ylabel('#Epitopes')
                        plt.savefig('epitope_distribution_plots/' + method + '_' + protein + '.png')
                        plt.clf()

        print("Epitope distribution analysis done.")
    
    elif choice == '2':
        
        t_cell_dict = defaultdict(list)
        for d in (mhci_dict, mhci_proc_dict, mhcii_dict):
            for key, value in d.items():
                t_cell_dict[key].append(value)
        for k, v in t_cell_dict.items():
            new_val = list(itertools.chain.from_iterable(v))
            t_cell_dict[k] = new_val
    
        b_cell_dict = defaultdict(list)
        for d in (bepipred2_dict, bepipred1_dict, emini_dict, choufasman_dict, karplusschulz_dict, kolaskartongaonkar_dict, parker_dict):
            for key, value in d.items():
                b_cell_dict[key].append(value)
        for k, v in b_cell_dict.items():
            new_val = list(itertools.chain.from_iterable(v))
            b_cell_dict[k] = new_val

        for key, value in t_cell_dict.items():
            start_pos = []
            end_pos = []
            for list_index, alist in enumerate(value):
                if list_index % 2 == 0:
                    for a_position in alist:
                        start_pos.append(a_position)
                else:
                    for a_position in alist:
                        end_pos.append(a_position)

            for swissprot_id, seq_length in protein_dict.items():
                if swissprot_id == key:
                    length = seq_length[1]
                
            count = [0]*length

            if len(start_pos) == len(end_pos):
                for i in range(0,len(start_pos)):
                    for j in range(start_pos[i],end_pos[i]):
                        count[j] = count[j]+1
                        
                x_axis = list(range(1,length+1))

                plt.plot(x_axis,count)
                plt.title("T-Cell " + key)
                plt.xlabel('Spike protein (' + key + ') residues')
                plt.ylabel('#Epitopes')
                plt.savefig('epitope_distribution_plots/' + "T-Cell" + '_' + key + '.png')
                plt.clf()
            else:
                print("Something went wrong")

        for key, value in b_cell_dict.items():
            start_pos = []
            end_pos = []
            for list_index, alist in enumerate(value):
                if list_index % 2 == 0:
                    for a_position in alist:
                        start_pos.append(a_position)
                else:
                    for a_position in alist:
                        end_pos.append(a_position)

            for swissprot_id, seq_length in protein_dict.items():
                if swissprot_id == key:
                    length = seq_length[1]
                
            count = [0]*length

            if len(start_pos) == len(end_pos):
                for i in range(0,len(start_pos)):
                    for j in range(start_pos[i],end_pos[i]):
                        count[j] = count[j]+1
                        
                x_axis = list(range(1,length+1))

                plt.plot(x_axis,count)
                plt.title("B-Cell " + key)
                plt.xlabel('Spike protein (' + key + ') residues')
                plt.ylabel('#Epitopes')
                plt.savefig('epitope_distribution_plots/' + "B-Cell" + '_' + key + '.png')
                plt.clf()
            else:
                print("Something went wrong")
            
        print("Epitope distribution analysis done.")

    elif choice == '3':

        all_dict = defaultdict(list)
        for d in (mhci_dict, mhci_proc_dict, mhcii_dict, bepipred2_dict, bepipred1_dict, emini_dict, choufasman_dict, karplusschulz_dict, kolaskartongaonkar_dict, parker_dict):
            for key, value in d.items():
                all_dict[key].append(value)
        for k, v in all_dict.items():
            new_val = list(itertools.chain.from_iterable(v))
            all_dict[k] = new_val

        for key, value in all_dict.items():
            start_pos = []
            end_pos = []
            for list_index, alist in enumerate(value):
                if list_index % 2 == 0:
                    for a_position in alist:
                        start_pos.append(a_position)
                else:
                    for a_position in alist:
                        end_pos.append(a_position)

            for swissprot_id, seq_length in protein_dict.items():
                if swissprot_id == key:
                    length = seq_length[1]
                
            count = [0]*length
        
            if len(start_pos) == len(end_pos):
                for i in range(0,len(start_pos)):
                    for j in range(start_pos[i],end_pos[i]):
                        count[j] = count[j]+1
                        
                x_axis = list(range(1,length+1))

                plt.plot(x_axis,count)
                plt.title("All of " + key)
                plt.xlabel('Spike protein (' + key + ') residues')
                plt.ylabel('#Epitopes')
                plt.savefig('epitope_distribution_plots/' + "all" + '_' + key + '.png')
                plt.clf()
            else:
                print("Something went wrong")

        print("Epitope distribution analysis done.")
    
    else:
        print("Please choose one of the options mentioned")
        epitope_distribution_plots(list_of_swissprot_ids)



def dssp_analysis(list_of_pdb_ids):
    
    """This function carries out DSSP analysis on a list of PDB IDs and returns plots that show the distribution of ACC values
    within the sequence of each of the PDB IDs."""

    options = Options()
    options.headless = True

    for id in list_of_pdb_ids:
        try:
            dssp_url = 'https://www3.cmbi.umcn.nl/xssp/'
            dssp = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
            dssp.get(dssp_url)
            #wait_dssp = WebDriverWait(dssp, 60)

            dssp.find_element(By.ID, "pdb_id").send_keys(id)
            time.sleep(3)
            dssp.find_element(By.XPATH, "/html/body/div/form/div[3]/button[2]").click()
            time.sleep(5)
            dssp_results = dssp.find_element(By.XPATH, '/html/body/div/div/div/textarea').text.splitlines()
            dssp.close()
            dssp_only_table = dssp_results[28:]

            acc_results = []
            for row in dssp_only_table:
                newrow = row.split()[:-9]
                try:
                    acc_results.append(int(newrow[13]))
                except:
                    try:
                        acc_results.append(int(newrow[12]))
                    except:
                        try:
                            acc_results.append(int(newrow[11]))
                        except:
                            try:
                                acc_results.append(int(newrow[10]))
                            except:
                                try:
                                    acc_results.append(int(newrow[9]))
                                except:
                                    try:
                                        acc_results.append(int(newrow[8]))
                                    except:
                                        try:
                                            acc_results.append(int(newrow[7]))
                                        except:
                                            try:
                                                acc_results.append(int(newrow[6]))
                                            except:
                                                try:
                                                    if newrow[1] == '!' or newrow[1] == '!*':
                                                        acc_results.append(0)
                                                except:
                                                    print("Failed to fetch DSSP output for: " + id)
                                                    break

            current_directory = os.getcwd()
            final_directory = os.path.join(current_directory, r'acc_results')
            if not os.path.exists(final_directory):
                os.makedirs(final_directory)

            lst = list(range(1,len(dssp_only_table)+1))
            plt.plot(lst, acc_results)
            plt.ylabel('ACC')
            plt.xlabel('Residue #')
            plt.savefig('acc_results/' + id + '.png')
            plt.clf()

        except:
            print("Failed to fetch DSSP output for: " + id)
            dssp.close()
    print("DSSP analysis done")

def analysis_choice(list_of_swissprot_ids, list_of_pdb_ids, counter = 10):

    """This function gives a choice to the user, whether DSSP analysis and epitope distribution analysis should be carried out."""

    if counter == 0:
        print("Failed after 10 tries")
        return

    choice = input("Do you want DSSP analysis and epitope distribution analysis?")

    if choice == 'y' or choice == 'Y' or choice == 'yes' or choice == 'Yes':
        protein_dict = swissprotIDSequenceLength(list_of_swissprot_ids)
        epitope_distribution_plots(protein_dict, list_of_swissprot_ids)
        dssp_analysis(list_of_pdb_ids)
    elif analysis_choice == 'n' or analysis_choice == 'N' or analysis_choice == 'no' or analysis_choice == 'No':
        pass
    else:
        print("Please enter either yes or no")
        analysis_choice(list_of_swissprot_ids, list_of_pdb_ids, counter - 1)


# list_of_swissprot_ids = ['P59594', 'P0DTC2', 'K9N5Q8', 'P36334', 'Q0ZME7', 'P15423', 'Q6Q1S2', 'Q5MQD0', 'Q14EB0']
# protein_dict = swissprotIDSequenceLength(list_of_swissprot_ids)
# epitope_distribution_plots(protein_dict, list_of_swissprot_ids)


