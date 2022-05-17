from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as ec
import scrapy
from scrapy import Spider
from scrapy.http import FormRequest
from scrapy.utils.project import get_project_settings
from scrapy.utils.log import configure_logging
from scrapy.crawler import CrawlerRunner
from twisted.internet import reactor, defer
import time
import pandas as pd
from msa import alignment
from msa import get_conserved_sequences
import requests


example_seq_dict = {'P0DTC2': ['SRSARSAIEDLLFDKVTIADPGYMQGYDDC','WSYTGSSFYAPEPITSLNTKY'],
'P36334': ['WMYTGSGYYYPEPITENNVVV','FKEELDQWFKNQTSVAPDL']}
example_seq = 'SYIVVGRGEQQINHHWHK'
list_of_swissprot_ids = ['P59594', 'P0DTC2', 'K9N5Q8', 'P36334', 'Q0ZME7', 'P15423', 'Q6Q1S2', 'Q5MQD0', 'Q14EB0']
path_to_mafft_alignment = '/home/erald/Desktop/ScrapyEpitope/msa_results/mafft/mafft.aln-fasta.fasta'
path_to_muscle_alignment = '/home/erald/Desktop/ScrapyEpitope/msa_results/muscle/muscle.aln-fasta.fasta'
mhci_alleles = ['HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*03:01']
mhci_lengths = [8, 9, 10]
mhcii_alleles = ['HLA-DRB1*01:01', 'HLA-DRB1*03:01']
mhcii_lengths = [13, 14, 16]
list_of_pdb_ids = ['6vxx', '6crz', '5x5f']


def mhci(conserved_sequences_dict, list_of_alleles, list_of_lengths):
    
    """This function uses the REST API from IEDB to access to MHC I Binding tool and predict the MHC Class I epitopes of the given conserved sequences. 
    It also needs a list of alleles and their respective lengths to run. It employs different methods to predict MHC Class I epitopes, including a 
    consensus approach which combines ANN, SMM, Comblib, NetMHCpan and Consensus. The results are returned as a list of lists"""

    alleles = ",".join(list_of_alleles)
    converted_list = [str(element) for element in list_of_lengths]
    lengths = ",".join(converted_list)

    mhci_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["allele"] + ["seq_num"] + ["start"] + ["end"] + ["length"] + ["peptide"] + ["core"] + ["icore"] + ["score"] + ["percentile_rank"]
    mhci_results.append(columns)

    for key in conserved_sequences_dict:
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
                
                        mhci_results.append(response_body[0])

                        df = pd.DataFrame(response_body[1:], columns=response_body[0])
                        df["percentile_rank"] = pd.to_numeric(df["percentile_rank"], errors='coerce')
                        df = df.loc[df['percentile_rank'] <= 10]

                        if df.empty == True:
                            pass
                        else:
                            rows = [[i for i in row[1:]] for row in df.itertuples()]
                            for i in rows:
                                i = [key] + [conserved_sequence] + i
                                mhci_results.append(i)
                    except:
                        pass

                else:
                    df = pd.DataFrame(response_body[1:], columns=response_body[0])
                    df["percentile_rank"] = pd.to_numeric(df["percentile_rank"], errors='coerce')
                    df = df.loc[df['percentile_rank'] <= 10]

                    if df.empty == True:
                        pass
                    else:
                        rows = [[i for i in row[1:]] for row in df.itertuples()]
                        for i in rows:
                            i = [key] + [conserved_sequence] + i
                            mhci_results.append(i)
    return mhci_results

def mhci_proc(conserved_sequences_dict, list_of_alleles, list_of_lengths):

    """This function uses the REST API from IEDB to access to MHC I Processing tool and predict the MHC Class I epitopes of the given conserved sequences. 
    It also needs a list of alleles and their respective lengths to run. It combines predictors of proteasomal processing, TAP transport, and MHC binding 
    to produce an overall score for each peptide's intrinsic potential of being a T cell epitope. The results are returned as a list of lists"""


    alleles = ",".join(list_of_alleles)
    converted_list = [str(element) for element in list_of_lengths]
    lengths = ",".join(converted_list)

    mhci_proc_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["allele"] + ["start"] + ["end"] + ['peptide_length'] + ["peptide"] + ["proteasome_score"] + ["tap_score"] + ["mhci_score"] + ["processing_score"] + ["total_score"] + ["mhci_ic50"]
    mhci_proc_results.append(columns)

    for key in conserved_sequences_dict:
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
                    pass
                else:
                    rows = [[i for i in row[1:]] for row in df.itertuples()]
                    for i in rows:
                        i = [key] + [conserved_sequence] + i
                        mhci_proc_results.append(i)

    return mhci_proc_results

def mhcii(conserved_sequences_dict, list_of_alleles, list_of_lengths):

    """This function uses the REST API from IEDB to access to MHC II Binding tool and predict the MHC Class II epitopes of the given conserved sequences. 
    It also needs a list of alleles and their respective lengths to run. It employs different methods to predict MHC Class II epitopes, including a consensus 
    approach which combines NN-align, SMM-align and Combinatorial library methods. The results are returned as a list of lists"""


    alleles = ",".join(list_of_alleles)
    converted_list = [str(element) for element in list_of_lengths]
    lengths = ",".join(converted_list)

    mhcii_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["allele"] + ["seq_num"] + ["start"] + ["end"] + ["length"] + ["core_peptide"] + ["peptide"] + ["ic50"] + ["rank"] + ["adjusted_rank"]
    mhcii_results.append(columns)

    for key in conserved_sequences_dict:
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
                df = df.loc[df['adjusted_rank'] <= 10]

                if df.empty == True:
                    pass
                else:
                    rows = [[i for i in row[1:]] for row in df.itertuples()]
                    for i in rows:
                        i = [key] + [conserved_sequence] + i
                        mhcii_results.append(i)

    return mhcii_results

def bepipred2(conserved_sequences_dict):

    """This function uses the REST API from IEDB to access to Linear Antigen Prediction method Bepipred 2.0 and predict the linear epitopes of the given conserved 
    sequences using a Random Forest algorithm trained on epitopes and non-epitope amino acids determined from crystal structures. The results are returned as a 
    list of lists"""

    bepipred2_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["predicted_epitope"] + ["start_position"] + ["end_position"]
    bepipred2_results.append(columns)

    for key in conserved_sequences_dict:
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
    
    return bepipred2_results

def bepipred(conserved_sequences_dict):
    
    """This function uses the REST API from IEDB to access to Linear Antigen Prediction method Bepipred and predict the linear epitopes of the given conserved 
    sequences using a combination of a hidden Markov model and a propensity scale method. The results are returned as a list of lists"""

    bepipred_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["predicted_epitope"] + ["start_position"] + ["end_position"]
    bepipred_results.append(columns)

    for key in conserved_sequences_dict:
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

    return bepipred_results

def emini(conserved_sequences_dict):
    
    """This function uses the REST API from IEDB to access the Linear Antigen Prediction method Emini and predict the linear epitopes of the given conserved 
    sequences. The calculation is based on surface accessibility scale on a product. The results are returned as a list of lists"""

    emini_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["predicted_epitope"] + ["start_position"] + ["end_position"]
    emini_results.append(columns)

    for key in conserved_sequences_dict:
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

    return emini_results

def choufasman(conserved_sequences_dict):

    """This function uses the REST API from IEDB to access the Linear Antigen Prediction method Chou-Fasman and predict the linear epitopes of the given 
    conserved sequences. It uses the Chou and Fasman scale which is commonly used to predict beta turns. The results are returned as a list of lists"""

    choufasman_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["predicted_epitope"] + ["start_position"] + ["end_position"]
    choufasman_results.append(columns)

    for key in conserved_sequences_dict:
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

    return choufasman_results

def karplusschulz(conserved_sequences_dict):

    """This function uses the REST API from IEDB to access the Linear Antigen Prediction method Karplus-Schulz and predict the linear epitopes of the given 
    conserved sequences. In this method, flexibility scale based on mobility of protein segments on the basis of the known temperature B factors of the 
    a-carbons of 31 proteins of known structure was constructed. The results are returned as a list of lists"""

    karplusschulz_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["predicted_epitope"] + ["start_position"] + ["end_position"]
    karplusschulz_results.append(columns)

    for key in conserved_sequences_dict:
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

    return karplusschulz_results

def kolaskartongaonkar(conserved_sequences_dict):

    """This function uses the REST API from IEDB to access the Linear Antigen Prediction method Kolaskar-Tongaonkar and predict the linear epitopes of the 
    given conserved sequences. It is a semi-empirical method which makes use of physicochemical properties of amino acid residues and their frequencies of 
    occurrence in experimentally known segmental epitopes was developed to predict antigenic determinants on proteins. The results are returned as a list of lists"""

    kolaskartongaonkar_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["predicted_epitope"] + ["start_position"] + ["end_position"]
    kolaskartongaonkar_results.append(columns)

    for key in conserved_sequences_dict:
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

    return kolaskartongaonkar_results

def parker(conserved_sequences_dict):

    """This function uses the REST API from IEDB to access the Linear Antigen Prediction method Parker and predict the linear epitopes of the given 
    conserved sequences. In this method, hydrophilic scale based on peptide retention times during high-performance liquid chromatography (HPLC) on a 
    reversed-phase column was constructed. A window of seven residues was used for analyzing epitope region. The results are returned as a list of lists"""

    parker_results = []
    columns = ["protein_id"] + ["conserved_sequence"] + ["predicted_epitope"] + ["start_position"] + ["end_position"]
    parker_results.append(columns)

    for key in conserved_sequences_dict:
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

    return parker_results

def ellipro(list_of_pdb_ids):
    
    """This function uses Selenium to access the Ellipro tool from IEDB. It requires a list of PDB IDs and returns two lists of lists, one with predicted linear 
    sequences, the other with predicted discontinous sequences. It implements a previously developed method that represents the protein structure as an 
    ellipsoid and calculates protrusion indexes for protein residues outside of the ellipsoid."""
    
    ellipro_url = 'http://tools.iedb.org/ellipro/'
    linear_columns = ['pdb_id','chain','start','end','peptide','nr_of_residues','score']
    discontinous_columns = ['pdb_id','chain','start','end','peptide','nr_of_residues','score']
    linear_epitopes = [linear_columns]
    discontinous_epitopes = [discontinous_columns]

    for protein_id in list_of_pdb_ids:
        driver = webdriver.Firefox(executable_path = '../ScrapyEpitope/geckodriver')
        driver.maximize_window()

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
            pass
        else:
            rows = [[i for i in row[1:]] for row in df.itertuples()]
            for i in rows:
                i = [protein_id] + i[1:]
                linear_epitopes.append(i)

        for i in range(len(table_of_discontinous_epitopes)):
            row= []
            for e in range(len(columns_of_discontinous_epitopes[:-1])):
                cell = driver.find_element(By.XPATH, '/html/body/div[3]/table[3]/tbody/tr[' + str(i+1) + ']/td[' + str(e+1) + ']').text
                row.append(cell)
            if float(row[3]) >= 0.7:
                unproc_discontinous_epitope = row[1].split(', ')
                processed_discontinous_epitope = ''
                start_pos = unproc_discontinous_epitope[0][3:]
                end_pos = unproc_discontinous_epitope[-1][3:]
                for position in unproc_discontinous_epitope:
                    aa = position[2]
                    processed_discontinous_epitope = processed_discontinous_epitope + aa
                row[1] = processed_discontinous_epitope
                row = [protein_id] + [chain]  + [start_pos] + [end_pos] + row[1:]
                discontinous_epitopes.append(row)
        driver.close()

    return linear_epitopes, discontinous_epitopes


def predict_all(list_of_swissprot_ids, path_to_alignment, alleles_for_mhci, lengths_for_mhci, alleles_for_mhcii, lengths_for_mhcii, list_of_pdb_ids):
    
    """This function runs all the prediction methods above and return a tuple with the lists of lists"""

    alignment(list_of_swissprot_ids, matrix='bl62', gapopen=1.53, gapext=0.123, order='aligned', nbtree=2, treeout='true', maxiterate=2, ffts='none')
    conserved_sequences_mafft = get_conserved_sequences(path_to_alignment, min_seq_conserved_pos='default', min_seq_flank_pos='default', max_contigous_nonconserved_pos = 8, min_length_block= 10, allowed_gap_pos='None')
    
    mhci_epitopes = mhci(example_seq_dict, alleles_for_mhci, lengths_for_mhci)
    mhci_proc_epitopes = mhci_proc(example_seq_dict, alleles_for_mhci, lengths_for_mhci)
    mhcii_epitopes = mhcii(example_seq_dict, alleles_for_mhcii, lengths_for_mhcii)
    
    bepipred2_epitopes = bepipred2(example_seq_dict)
    bepipred_epitopes = bepipred(example_seq_dict)
    emini_epitopes = emini(example_seq_dict)
    choufasman_epitopes = choufasman(example_seq_dict)
    karplusschulz_epitopes = karplusschulz(example_seq_dict)
    kolaskartongaonkar_epitopes = kolaskartongaonkar(example_seq_dict)
    parker_epitopes = parker(example_seq_dict)

    ellipro_epitopes = ellipro(list_of_pdb_ids)
    
    return mhci_epitopes, mhci_proc_epitopes, mhcii_epitopes, bepipred2_epitopes, bepipred_epitopes, emini_epitopes, choufasman_epitopes, karplusschulz_epitopes, kolaskartongaonkar_epitopes, parker_epitopes, ellipro_epitopes


#predicted_epitopes = predict_all(list_of_swissprot_ids, path_to_mafft_alignment, mhci_alleles, mhci_lengths, mhcii_alleles, mhcii_lengths, list_of_pdb_ids)
#for alist in predicted_epitopes:
#    print(*alist, sep='\n')


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


