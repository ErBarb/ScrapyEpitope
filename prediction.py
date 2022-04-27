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


example_seq_dict = {'P0DTC2': ['FNCLGMSNRDFLEAGHAHKVPRRLLKAHKVPRDADGSGSGSG','GATQAGRFSITP']}
example_seq = 'SYIVVGRGEQQINHHWHK'

def mhci_curl(conserved_sequences_dict):
    
    mhci_results = []
    for key in conserved_sequences_dict:
        for conserved_sequence in conserved_sequences_dict[key]:
            print(conserved_sequence)
            headers = {
                'Content-Type': 'application/x-www-form-urlencoded',
            }
            
            data = 'method=recommended&sequence_text=' + conserved_sequence + '&allele=HLA-A*01:01,HLA-A*02:01&length=8,9'
            #data = 'method=ann&sequence_text=%3Epeptide1%0AGHAHKVPRRLLKAAR%0A%3Epeptide2%0ALKAADASADADGSGSGSGSG&allele=HLA-A*01:01,HLA-A*02:01&length=8,9'
            
            response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/mhci/', headers=headers, data=data)
            response_data_split_by_line = response.content.decode('utf-8').splitlines()
            
            response_body = []
            for line in response_data_split_by_line:
                split_line = line.split("\t")
                response_body.append(split_line)
            
            mhci_results.append(response_body[0])

            try:
                df = pd.DataFrame(response_body[1:], columns=response_body[0])

            except:
                print('Retrying prediction for sequence: ' + conserved_sequence)

                try:
                    headers = {
                        'Content-Type': 'application/x-www-form-urlencoded',
                    }
            
                    data = 'method=recommended&sequence_text=' + conserved_sequence + '&allele=HLA-A*01:01,HLA-A*02:01&length=8,9'
                    #data = 'method=ann&sequence_text=%3Epeptide1%0AGHAHKVPRRLLKAAR%0A%3Epeptide2%0ALKAADASADADGSGSGSGSG&allele=HLA-A*01:01,HLA-A*02:01&length=8,9'
            
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
                        print('DataFrame is empty! No epitopes passed the threshold.')
                    else:
                        print(df)
                        rows = [[i for i in row[1:]] for row in df.itertuples()]
                        for i in rows:
                            i = [key] + [conserved_sequence] + i
                            mhci_results.append(i)
                except:
                    print("Epitope prediction for sequence " + conserved_sequence + " failed")

            else:
                df = pd.DataFrame(response_body[1:], columns=response_body[0])
                df["percentile_rank"] = pd.to_numeric(df["percentile_rank"], errors='coerce')
                df = df.loc[df['percentile_rank'] <= 10]

                if df.empty == True:
                    print('DataFrame is empty! No epitopes passed the threshold.')
                else:
                    print(df)
                    rows = [[i for i in row[1:]] for row in df.itertuples()]
                    for i in rows:
                        i = [key] + [conserved_sequence] + i
                        mhci_results.append(i)
    return mhci_results

def mhcii_curl(conserved_sequences_dict):
    #Length of input sequence needs to be greater than 16 aa
    mhcii_results = []

    for key in conserved_sequences_dict:
        for conserved_sequence in conserved_sequences_dict[key]:
            if len(conserved_sequence) > 16:
                headers = {
                    'Content-Type': 'application/x-www-form-urlencoded',
                }
                
                data = 'method=nn_align&sequence_text=' + conserved_sequence + '&allele=HLA-DRB1*01:01&length=13,14,16'
                #data = 'method=ann&sequence_text=%3Epeptide1%0AGHAHKVPRRLLKAAR%0A%3Epeptide2%0ALKAADASADADGSGSGSGSG&allele=HLA-A*01:01,HLA-A*02:01&length=8,9'
                
                response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/mhcii/', headers=headers, data=data)
                response_data_split_by_line = response.content.decode('utf-8').splitlines()
                
                response_body = []
                for line in response_data_split_by_line:
                    split_line = line.split("\t")
                    response_body.append(split_line)
                
                mhcii_results.append(response_body[0])

                df = pd.DataFrame(response_body[1:], columns=response_body[0])
                df["adjusted_rank"] = pd.to_numeric(df["adjusted_rank"])
                df = df.loc[df['adjusted_rank'] <= 100]

                rows = [[i for i in row[1:]] for row in df.itertuples()]
                for i in rows:
                    i = [key] + [conserved_sequence] + i
                    mhcii_results.append(i)

    return mhcii_results

def bepipred2_curl(conserved_sequences_dict):
    bepipred2_results = []

    for key in conserved_sequences_dict:
        for conserved_sequence in conserved_sequences_dict[key]:
            #if len(conserved_sequence) > 16:
            headers = {
                'Content-Type': 'application/x-www-form-urlencoded',
            }
            
            data = {
                'method': 'Bepipred-2.0',
                'sequence_text': conserved_sequence,
                'window_size': '7',
}
            #data = 'method=ann&sequence_text=%3Epeptide1%0AGHAHKVPRRLLKAAR%0A%3Epeptide2%0ALKAADASADADGSGSGSGSG&allele=HLA-A*01:01,HLA-A*02:01&length=8,9'
            
            response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/bcell/', headers=headers, data=data)
            response_data_split_by_line = response.content.decode('utf-8').splitlines()
            
            response_body = []
            for line in response_data_split_by_line:
                split_line = line.split("\t")
                response_body.append(split_line)
            
            bepipred2_results.append(response_body[0])

            df = pd.DataFrame(response_body[1:], columns=response_body[0])
            df["Score"] = pd.to_numeric(df["Score"])
            df = df.loc[df['Score'] >= 0.5]

            rows = [[i for i in row[1:]] for row in df.itertuples()]
            for i in rows:
                i = [key] + [conserved_sequence] + i
                bepipred2_results.append(i)

    return bepipred2_results

def bepipred_curl(conserved_sequences_dict):
    bepipred_results = []

    for key in conserved_sequences_dict:
        for conserved_sequence in conserved_sequences_dict[key]:
            #if len(conserved_sequence) > 16:
            headers = {
                'Content-Type': 'application/x-www-form-urlencoded',
            }
            
            data = {
                'method': 'Bepipred',
                'sequence_text': conserved_sequence,
                'window_size': '7',
}
            #data = 'method=ann&sequence_text=%3Epeptide1%0AGHAHKVPRRLLKAAR%0A%3Epeptide2%0ALKAADASADADGSGSGSGSG&allele=HLA-A*01:01,HLA-A*02:01&length=8,9'
            
            response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/bcell/', headers=headers, data=data)
            response_data_split_by_line = response.content.decode('utf-8').splitlines()
            
            response_body = []
            for line in response_data_split_by_line:
                split_line = line.split("\t")
                response_body.append(split_line)
            
            bepipred_results.append(response_body[0])

            df = pd.DataFrame(response_body[1:], columns=response_body[0])
            df["Score"] = pd.to_numeric(df["Score"])
            df = df.loc[df['Score'] >= 0.35]

            rows = [[i for i in row[1:]] for row in df.itertuples()]
            for i in rows:
                i = [key] + [conserved_sequence] + i
                bepipred_results.append(i)

    return bepipred_results

def emini_curl(conserved_sequences_dict):
    emini_results = []

    for key in conserved_sequences_dict:
        for conserved_sequence in conserved_sequences_dict[key]:
            #if len(conserved_sequence) > 16:
            headers = {
                'Content-Type': 'application/x-www-form-urlencoded',
            }
            
            data = {
                'method': 'Emini',
                'sequence_text': conserved_sequence,
                'window_size': '6',
}
            #data = 'method=ann&sequence_text=%3Epeptide1%0AGHAHKVPRRLLKAAR%0A%3Epeptide2%0ALKAADASADADGSGSGSGSG&allele=HLA-A*01:01,HLA-A*02:01&length=8,9'
            
            response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/bcell/', headers=headers, data=data)
            response_data_split_by_line = response.content.decode('utf-8').splitlines()
            
            response_body = []
            for line in response_data_split_by_line:
                split_line = line.split("\t")
                response_body.append(split_line)
            
            emini_results.append(response_body[0])

            df = pd.DataFrame(response_body[1:], columns=response_body[0])
            df["Score"] = pd.to_numeric(df["Score"])
            df = df.loc[df['Score'] >= 1]

            rows = [[i for i in row[1:]] for row in df.itertuples()]
            for i in rows:
                i = [key] + [conserved_sequence] + i
                emini_results.append(i)

    return emini_results

def choufasman_curl(conserved_sequences_dict):
    #Figure out how to find threshold, maybe just get average of lowest and highest score
    choufasman_results = []

    for key in conserved_sequences_dict:
        for conserved_sequence in conserved_sequences_dict[key]:
            #if len(conserved_sequence) > 16:
            headers = {
                'Content-Type': 'application/x-www-form-urlencoded',
            }
            
            data = {
                'method': 'Chou-Fasman',
                'sequence_text': conserved_sequence,
                'window_size': '7',
}
            #data = 'method=ann&sequence_text=%3Epeptide1%0AGHAHKVPRRLLKAAR%0A%3Epeptide2%0ALKAADASADADGSGSGSGSG&allele=HLA-A*01:01,HLA-A*02:01&length=8,9'
            
            response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/bcell/', headers=headers, data=data)
            response_data_split_by_line = response.content.decode('utf-8').splitlines()
            
            response_body = []
            for line in response_data_split_by_line:
                split_line = line.split("\t")
                response_body.append(split_line)
            
            choufasman_results.append(response_body[0])

            df = pd.DataFrame(response_body[1:], columns=response_body[0])
            df["Score"] = pd.to_numeric(df["Score"])
            df = df.loc[df['Score'] >= 1]

            rows = [[i for i in row[1:]] for row in df.itertuples()]
            for i in rows:
                i = [key] + [conserved_sequence] + i
                choufasman_results.append(i)

    return choufasman_results

def karplusschulz_curl(conserved_sequences_dict):
    #Figure out how to find threshold, maybe just get average of lowest and highest score
    karplusschulz_results = []

    for key in conserved_sequences_dict:
        for conserved_sequence in conserved_sequences_dict[key]:
            #if len(conserved_sequence) > 16:
            headers = {
                'Content-Type': 'application/x-www-form-urlencoded',
            }
            
            data = {
                'method': 'Karplus-Schulz',
                'sequence_text': conserved_sequence,
                'window_size': '7',
}
            #data = 'method=ann&sequence_text=%3Epeptide1%0AGHAHKVPRRLLKAAR%0A%3Epeptide2%0ALKAADASADADGSGSGSGSG&allele=HLA-A*01:01,HLA-A*02:01&length=8,9'
            
            response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/bcell/', headers=headers, data=data)
            response_data_split_by_line = response.content.decode('utf-8').splitlines()
            
            response_body = []
            for line in response_data_split_by_line:
                split_line = line.split("\t")
                response_body.append(split_line)
            
            karplusschulz_results.append(response_body[0])

            df = pd.DataFrame(response_body[1:], columns=response_body[0])
            df["Score"] = pd.to_numeric(df["Score"])
            df = df.loc[df['Score'] >= 1]

            rows = [[i for i in row[1:]] for row in df.itertuples()]
            for i in rows:
                i = [key] + [conserved_sequence] + i
                karplusschulz_results.append(i)

    return karplusschulz_results

def kolaskartongaonkar_curl(conserved_sequences_dict):
    #Figure out how to find threshold, maybe just get average of lowest and highest score
    kolaskartongaonkar_results = []

    for key in conserved_sequences_dict:
        for conserved_sequence in conserved_sequences_dict[key]:
            #if len(conserved_sequence) > 16:
            headers = {
                'Content-Type': 'application/x-www-form-urlencoded',
            }
            
            data = {
                'method': 'Kolaskar-Tongaonkar',
                'sequence_text': conserved_sequence,
                'window_size': '7',
}
            #data = 'method=ann&sequence_text=%3Epeptide1%0AGHAHKVPRRLLKAAR%0A%3Epeptide2%0ALKAADASADADGSGSGSGSG&allele=HLA-A*01:01,HLA-A*02:01&length=8,9'
            
            response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/bcell/', headers=headers, data=data)
            response_data_split_by_line = response.content.decode('utf-8').splitlines()
            
            response_body = []
            for line in response_data_split_by_line:
                split_line = line.split("\t")
                response_body.append(split_line)
            
            kolaskartongaonkar_results.append(response_body[0])

            df = pd.DataFrame(response_body[1:], columns=response_body[0])
            df["Score"] = pd.to_numeric(df["Score"])
            df = df.loc[df['Score'] >= 1]

            rows = [[i for i in row[1:]] for row in df.itertuples()]
            for i in rows:
                i = [key] + [conserved_sequence] + i
                kolaskartongaonkar_results.append(i)

    return kolaskartongaonkar_results

def parker_curl(conserved_sequences_dict):
    #Figure out how to find threshold, maybe just get average of lowest and highest score
    parker_results = []

    for key in conserved_sequences_dict:
        for conserved_sequence in conserved_sequences_dict[key]:
            #if len(conserved_sequence) > 16:
            headers = {
                'Content-Type': 'application/x-www-form-urlencoded',
            }
            
            data = {
                'method': 'Parker',
                'sequence_text': conserved_sequence,
                'window_size': '7',
}
            #data = 'method=ann&sequence_text=%3Epeptide1%0AGHAHKVPRRLLKAAR%0A%3Epeptide2%0ALKAADASADADGSGSGSGSG&allele=HLA-A*01:01,HLA-A*02:01&length=8,9'
            
            response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/bcell/', headers=headers, data=data)
            response_data_split_by_line = response.content.decode('utf-8').splitlines()
            
            response_body = []
            for line in response_data_split_by_line:
                split_line = line.split("\t")
                response_body.append(split_line)
            
            parker_results.append(response_body[0])

            df = pd.DataFrame(response_body[1:], columns=response_body[0])
            df["Score"] = pd.to_numeric(df["Score"])
            df = df.loc[df['Score'] >= 1.5]

            rows = [[i for i in row[1:]] for row in df.itertuples()]
            for i in rows:
                i = [key] + [conserved_sequence] + i
                parker_results.append(i)

    return parker_results

def predict_all():
    list_of_swissprot_ids = ['P59594', 'P0DTC2', 'K9N5Q8', 'P36334', 'Q0ZME7', 'P15423', 'Q6Q1S2', 'Q5MQD0', 'Q14EB0']
    #alignment(list_of_swissprot_ids)
    path_to_mafft_alignment = '/home/erald/Desktop/ScrapyEpitope/msa_results/mafft/mafft.aln-fasta.fasta'
    conserved_sequences = get_conserved_sequences(path_to_mafft_alignment, min_seq_conserved_pos='default', min_seq_flank_pos='default', max_contigous_nonconserved_pos = 8, min_length_block= 10, allowed_gap_pos='None')
    mhci_epitopes = mhci_curl(conserved_sequences)
    #mhcii_epitopes = mhcii_curl(conserved_sequences)
    #bepipred2_epitopes = bepipred2_curl(conserved_sequences)
    #bepipred_epitopes = bepipred_curl(conserved_sequences)
    #emini_epitopes = emini_curl(conserved_sequences)
    #choufasman_epitopes = choufasman_curl(conserved_sequences)
    #karplusschulz_epitopes = karplusschulz_curl(conserved_sequences)
    #kolaskartongaonkar_epitopes = kolaskartongaonkar_curl(conserved_sequences)
    #parker_epitopes = parker_curl(conserved_sequences)
    return mhci_epitopes 
    #mhcii_epitopes, bepipred2_epitopes, bepipred_epitopes, emini_epitopes, choufasman_epitopes, karplusschulz_epitopes, kolaskartongaonkar_epitopes, parker_epitopes


predicted_epitopes = predict_all()
#print(*predicted_epitopes, sep='\n')


def mhci(conserved_sequences_dict, mhci_alleles=None, mhci_lengths=None):

    """This function uses Selenium to access the MHCI prediction tool from IEDB and returns a list of lists of the
    predicted epitopes and their attributes. Its arguments are the dictionary with the conserved sequences and their
    protein IDs, and lists of alleles and lengths of epitopes for each allele(coming soon)"""

    mhci_url = 'http://tools.iedb.org/mhci/'
    columns = ['protein_id', 'conserved_sequence', 'method_used', 'allele', 'seq_num', 'start', 'end', 'length', 'peptide', 'score', 'percentile_rank']
    mhci_results = [columns]

    for key in conserved_sequences_dict:
        for conserved_sequence in conserved_sequences_dict[key]:
            driver = webdriver.Firefox(executable_path = '../ScrapyEpitope/geckodriver')
            driver.maximize_window()
            driver.get(mhci_url)

            driver.find_element(By.NAME, "sequence_text").send_keys(conserved_sequence)
            driver.find_element(By.ID, "id_refset").click()
            time.sleep(3)
            # driver.find_element(By.XPATH, "//select[@name='allele_list']/option[text()='HLA-A*01:01']").click()
            # time.sleep(1)
            # driver.find_element(By.XPATH, "//select[@name='length_list']/option[text()='10']").click()
            # time.sleep(1)

            driver.find_element(By.XPATH, "//select[@name='output_format']/option[text()='Text file']").click()
            driver.find_element(By.XPATH, "/html/body/div[3]/form/table/tbody/tr[14]/th/div/input[2]").click()

            wait = WebDriverWait(driver, 180)
            wait.until(ec.visibility_of_element_located((By.XPATH, "/html/body/pre")))

            text_body = driver.find_element(By.XPATH, '/html/body/pre').text.splitlines()
            driver.close()

            method_used = text_body[2].split(' ')[2]

            text_body_split_rows = []
            for row in text_body[3:-3]:
                split_row = row.split(" ")
                text_body_split_rows.append(split_row)

            df = pd.DataFrame(text_body_split_rows[1:], columns=text_body_split_rows[0])
            df["percentile_rank"] = pd.to_numeric(df["percentile_rank"])
            df = df.loc[df['percentile_rank'] <= 1]

            rows = [[i for i in row[1:]] for row in df.itertuples()]
            for i in rows:
                i = [key] + [conserved_sequence] + [method_used] + i
                mhci_results.append(i)

    return mhci_results


class Bepipred(Spider):
    name = 'bepipred'
    start_urls = ['http://tools.iedb.org/bcell']
    global bepipred

    def parse(self, response):
        yield FormRequest.from_response(
            response,
            url=self.start_urls[0],
            formdata={
                'swissprot': 'P0C6U8',
                'sequence_text': '',
                'method': 'Bepipred'
            },
            callback=self.get_results)

    def get_results(self, results_page):
        row = []
        for cell in results_page.xpath('/html/body/div[3]/table[2]/tbody/tr'):
            if int(cell.xpath('td[5]//text()').extract_first()) >= 7:
                row.append(cell.xpath('td[1]//text()').extract_first())
                row.append(cell.xpath('td[2]//text()').extract_first())
                row.append(cell.xpath('td[3]//text()').extract_first())
                row.append(cell.xpath('td[4]//text()').extract_first())
                row.append(cell.xpath('td[5]//text()').extract_first())
                bepipred.append(row)
                row = []
class Emini(Spider):
    name = 'Emini'
    start_urls = ['http://tools.iedb.org/bcell']
    global emini

    def parse(self, response):
        yield FormRequest.from_response(
            response,
            url=self.start_urls[0],
            formdata={
                'swissprot': 'P0C6U8',
                'sequence_text': '',
                'method': 'Emini'
            },
            callback=self.get_results)

    def get_results(self, results_page):
        row = []
        for cell in results_page.xpath('/html/body/div[3]/table[2]/tbody/tr'):
            if int(cell.xpath('td[5]//text()').extract_first()) >= 7:
                row.append(cell.xpath('td[1]//text()').extract_first())
                row.append(cell.xpath('td[2]//text()').extract_first())
                row.append(cell.xpath('td[3]//text()').extract_first())
                row.append(cell.xpath('td[4]//text()').extract_first())
                row.append(cell.xpath('td[5]//text()').extract_first())
                emini.append(row)
                row = []
class Kolaskar(Spider):
    name = 'Kolaskar-Tongaonkar'
    start_urls = ['http://tools.iedb.org/bcell']
    global kolaskar_tongaonkar

    def parse(self, response):
        yield FormRequest.from_response(
            response,
            url=self.start_urls[0],
            formdata={
                'swissprot': 'P0C6U8',
                'sequence_text': '',
                'method': 'Kolaskar-Tongaonkar'
            },
            callback=self.get_results)

    def get_results(self, results_page):
        row = []
        for cell in results_page.xpath('/html/body/div[3]/table[2]/tbody/tr'):
            if int(cell.xpath('td[5]//text()').extract_first()) >= 7:
                row.append(cell.xpath('td[1]//text()').extract_first())
                row.append(cell.xpath('td[2]//text()').extract_first())
                row.append(cell.xpath('td[3]//text()').extract_first())
                row.append(cell.xpath('td[4]//text()').extract_first())
                row.append(cell.xpath('td[5]//text()').extract_first())
                kolaskar_tongaonkar.append(row)
                row = []
class Bepipred2(Spider):
    name = 'bepipred2'
    start_urls = ['http://tools.iedb.org/bcell']
    global bepipred2

    def parse(self, response):
        yield FormRequest.from_response(
            response,
            url=self.start_urls[0],
            formdata={
                'swissprot': 'P02185',
                'sequence_text': '',
                'method': 'Bepipred2'
            },
            callback=self.get_results)

    def get_results(self, results_page):
        row = []
        for cell in results_page.xpath('/html/body/div[3]/table[2]/tbody/tr'):
            if int(cell.xpath('td[5]//text()').extract_first()) >= 7:
                row.append(cell.xpath('td[1]//text()').extract_first())
                row.append(cell.xpath('td[2]//text()').extract_first())
                row.append(cell.xpath('td[3]//text()').extract_first())
                row.append(cell.xpath('td[4]//text()').extract_first())
                row.append(cell.xpath('td[5]//text()').extract_first())
                bepipred2.append(row)
                row = []
class Chou_Fasman(Spider):
    name = 'Chou-Fasman'
    start_urls = ['http://tools.iedb.org/bcell']
    global choufasman

    def parse(self, response):
        yield FormRequest.from_response(
            response,
            url=self.start_urls[0],
            formdata={
                'swissprot': 'P0C6U8',
                'sequence_text': '',
                'method': 'Chou-Fasman'
            },
            callback=self.get_results)

    def get_results(self, results_page):
        row = []
        for cell in results_page.xpath('/html/body/div[3]/table[2]/tbody/tr'):
            if cell.xpath('td[5]//text()').extract_first():
                row.append(cell.xpath('td[1]//text()').extract_first())
                row.append(cell.xpath('td[2]//text()').extract_first())
                row.append(cell.xpath('td[3]//text()').extract_first())
                row.append(cell.xpath('td[4]//text()').extract_first())
                row.append(cell.xpath('td[5]//text()').extract_first())
                row.append(cell.xpath('td[6]//text()').extract_first())
                choufasman.append(row)
                row = []
class Karplus_Schulz(Spider):
    name = 'Karplus-Schulz'
    start_urls = ['http://tools.iedb.org/bcell']
    global karplusschulz

    def parse(self, response):
        yield FormRequest.from_response(
            response,
            url=self.start_urls[0],
            formdata={
                'swissprot': 'P0C6U8',
                'sequence_text': '',
                'method': 'Karplus-Schulz'
            },
            callback=self.get_results)

    def get_results(self, results_page):
        row = []
        for cell in results_page.xpath('/html/body/div[3]/table[2]/tbody/tr'):
            if cell.xpath('td[5]//text()').extract_first():
                row.append(cell.xpath('td[1]//text()').extract_first())
                row.append(cell.xpath('td[2]//text()').extract_first())
                row.append(cell.xpath('td[3]//text()').extract_first())
                row.append(cell.xpath('td[4]//text()').extract_first())
                row.append(cell.xpath('td[5]//text()').extract_first())
                row.append(cell.xpath('td[6]//text()').extract_first())
                karplusschulz.append(row)
                row = []
class Parker(Spider):
    name = 'Parker'
    start_urls = ['http://tools.iedb.org/bcell']
    global parker

    def parse(self, response):
        yield FormRequest.from_response(
            response,
            url=self.start_urls[0],
            formdata={
                'swissprot': 'P0C6U8',
                'sequence_text': '',
                'method': 'Parker'
            },
            callback=self.get_results)

    def get_results(self, results_page):
        row = []
        for cell in results_page.xpath('/html/body/div[3]/table[2]/tbody/tr'):
            if cell.xpath('td[5]//text()').extract_first():
                row.append(cell.xpath('td[1]//text()').extract_first())
                row.append(cell.xpath('td[2]//text()').extract_first())
                row.append(cell.xpath('td[3]//text()').extract_first())
                row.append(cell.xpath('td[4]//text()').extract_first())
                row.append(cell.xpath('td[5]//text()').extract_first())
                row.append(cell.xpath('td[6]//text()').extract_first())
                parker.append(row)
                row = []
class Discotope(Spider):
    name = 'Discotope'
    start_urls = ['http://tools.iedb.org/discotope/']
    global discotope

    def parse(self, response):
        yield FormRequest.from_response(
            response,
            url=self.start_urls[0],
            formdata={
                'pdb': '6vxx',
                'chain': 'A',
                'version': '1.1'
            },
            callback=self.results_page)

    def results_page(self, response):
        url = 'http://tools.iedb.org/discotope/table/'
        yield scrapy.Request(url=url, callback=self.get_results, dont_filter=True)

    def get_results(self, results):
        row = []
        for cell in results.xpath('/html/body/div[3]/table/tbody/tr'):
            row.append(cell.xpath('td[1]//text()').extract_first())
            row.append(cell.xpath('td[2]//text()').extract_first())
            row.append(cell.xpath('td[3]//text()').extract_first())
            row.append(cell.xpath('td[4]//text()').extract_first())
            row.append(cell.xpath('td[5]//text()').extract_first())
            row.append(cell.xpath('td[6]//text()').extract_first())
            discotope.append(row)
            row = []
class Ellipro(Spider):
    name = 'Ellipro'
    start_urls = ['http://tools.iedb.org/ellipro/']
    global ellipro

    def parse(self, response):
        yield FormRequest.from_response(
            response,
            url=self.start_urls[0],
            formdata={
                'pred_tool': 'ellipro',
                'source': 'html',
                'pdb_id': '6vxx',
                'pdb_file': '(binary)',
                'min_score': '0.5',
                'max_distance': '6',
                'submit': 'Submit'
            },
            callback=self.second_page)

    def second_page(self, sec_page):
        yield FormRequest.from_response(
            sec_page,
            formdata={
                'chain': ['A', 'B', 'C']
            },
            callback=self.get_results)

    def get_results(self, results_page):
        row = []
        for cell in results_page.xpath('/html/body/div[3]/table[2]/tbody/tr'):
            row.append(cell.xpath('td[1]//text()').extract_first())
            row.append(cell.xpath('td[2]//text()').extract_first())
            row.append(cell.xpath('td[3]//text()').extract_first())
            row.append(cell.xpath('td[4]//text()').extract_first())
            row.append(cell.xpath('td[5]//text()').extract_first())
            row.append(cell.xpath('td[6]//text()').extract_first())
            row.append(cell.xpath('td[7]//text()').extract_first())
            row.append(cell.xpath('td[8]//text()').extract_first())
            ellipro.append(row)
            row = []
class MhcII(Spider):
    name = 'MhcII'
    start_urls = ['http://tools.iedb.org/mhcii/']
    global mhcii
    global fasta
    global alleles_mhcii
    global length_mhcii

    def parse(self, response):
        yield FormRequest.from_response(
            response,
            url=self.start_urls[0],
            formdata={
                'sequence_text': fasta,
                'method': '1222',
                'locus_list': 'DR',
                'refset': 'on',
                'allele': alleles_mhcii,
                'allele_list': '',
                'allele_list_a': '',
                'allele_file': '(binary)',
                'len-select-options': 'mul_rec',
                'length': length_mhcii,
                'sort_output': 'adjusted_rank',
                'output_format': 'ascii'
            },
            callback=self.get_results)

    def get_results(self, results_page):
        text = results_page.xpath('/html/body/pre//text()').extract_first()
        rows = text.splitlines()
        mhcii = []
        for row in rows:
            split_row = row.split("\t")
            mhcii.append(split_row)
            print(split_row)
class Netctlpan(Spider):
    name = 'Netctlpan'
    start_urls = ['http://tools.iedb.org/netchop/']
    global netctlpan
    global fasta

    def parse(self, response):
        yield FormRequest.from_response(
            response,
            url=self.start_urls[0],
            formdata={
                'pred_tool': 'netchop',
                'pred_method': 'netctlpan',
                'sequence_text': fasta,
                'sequence_file': '(binary)',
                'method': '0',
                'netchop_threshold': '0.5',
                'netctl_cleavage': '0.15',
                'netctl_tap': '0.05',
                'supertype': 'A1',
                'netctl_threshold': '0.75',
                'species_list': 'human',
                'freq': 'freq',
                'allele_list': 'HLA-A01:01',
                'length_list': '9',
                'netctlpan_threshold': '-99.9',
                'netctlpan_cleavage': '0.225',
                'netctlpan_tap': '0.025',
                'epitope_threshold': '1.0'
            },
            callback=self.results_page)

    def results_page(self, response):
        url = 'http://tools.iedb.org/netchop/table/'
        yield scrapy.Request(url=url, callback=self.get_results, dont_filter=True)

    def get_results(self, results):
        row = []
        for cell in results.xpath('/html/body/div[3]/table/tbody/tr'):
            row.append(cell.xpath('td[1]//text()').extract_first())
            row.append(cell.xpath('td[2]//text()').extract_first())
            row.append(cell.xpath('td[3]//text()').extract_first())
            row.append(cell.xpath('td[4]//text()').extract_first())
            row.append(cell.xpath('td[5]//text()').extract_first())
            row.append(cell.xpath('td[6]//text()').extract_first())
            netctlpan.append(row)
            row = []
        print(netctlpan)


#configure_logging()
#settings = get_project_settings()
#runner = CrawlerRunner(settings)


@defer.inlineCallbacks
def crawl():
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
    reactor.stop()


#crawl()
#reactor.run()  # the script will block here until the last crawl call is finished


