from unittest import expectedFailure
import requests
import re
import time
import urllib
import os
import matplotlib.pyplot as plt
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as ec
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.keys import Keys

swissprot_example_ids = ['P59594', 'P0DTC2', 'K9N5Q8', 'P36334', 'Q0ZME7', 'P15423', 'Q6Q1S2', 'Q5MQD0', 'Q14EB0']

options = Options()
options.headless = True

# for id in swissprot_example_ids:
#         baseUrl="http://www.uniprot.org/uniprot/"
#         currentUrl=baseUrl+id+".fasta"
#         response = requests.post(currentUrl)
#         cData=''.join(response.text)
#         print(len(cData))

# pdb_ids = []
# for id in swissprot_example_ids:
#     try:
#         uniprot_url = 'https://www.uniprot.org/uniprot/' + id
#         uniprot = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
#         uniprot.get(uniprot_url)
#         wait_uniprot = WebDriverWait(uniprot, 30)
#         wait_uniprot.until(ec.visibility_of_element_located((By.XPATH, '//div[@id="structure"]/protvista-uniprot-structure/div/div[2]/protvista-datatable/table/tbody')))
#         uniprot_pdb_ids = uniprot.find_element(By.XPATH, '//div[@id="structure"]/protvista-uniprot-structure/div/div[2]/protvista-datatable/table/tbody').text.splitlines()
#         for row in uniprot_pdb_ids:
#             newrow = row.split()
#             if (newrow[2] == 'EM' or newrow[2] == 'X-ray') and (newrow[6].startswith('1-') == True):
#                 print(newrow[2], newrow[1],newrow[6])
#                 print(id)
#                 pdb_ids.append(newrow[1])
#                 break
#             elif (newrow[2] == 'EM' or newrow[2] == 'X-ray') and re.search("^1.-", newrow[6]):
#                 print(newrow[2], newrow[1],newrow[6])
#                 #pdb_ids.append(newrow[1])
#             else:
#                 continue
#         uniprot.close()
#     except:
#         print("No PDB format for: " + id)
#         uniprot.close()

# print(pdb_ids)

my_pdb_ids = ['5WRG','6VSB','5W9H','5I08','6U7H','7KIP']

for id in my_pdb_ids:

    # try:
    #     url = "https://files.rcsb.org/download/" + id + ".pdb"
    #     urllib.request.urlretrieve(url, os.getcwd() + '/' + id + ".pdb")
    # except:
    #     print("No PDB format for: " + id)

    try:
        dssp_url = 'https://www3.cmbi.umcn.nl/xssp/'
        dssp = webdriver.Firefox(options=options, executable_path = '../ScrapyEpitope/geckodriver')
        dssp.get(dssp_url)
        wait_dssp = WebDriverWait(dssp, 60)

        
        # dssp.find_element(By.XPATH, "/html/body/div/form/div[1]/div[2]/div/select/option[3]").click()
        # dssp.find_element(By.ID, "file_").send_keys(os.getcwd() + "/" + id + ".pdb")

        dssp.find_element(By.ID, "pdb_id").send_keys(id)
        time.sleep(3)
        dssp.find_element(By.XPATH, "/html/body/div/form/div[3]/button[2]").click()
        # time.sleep(90)
        time.sleep(5)
        dssp_results = dssp.find_element(By.XPATH, '/html/body/div/div/div/textarea').text.splitlines()
        dssp.close()
        dssp_only_table = dssp_results[28:]
        # os.remove(os.getcwd() + "/" + id + ".pdb")
        
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
                                                print("Failed")
                                                print(newrow)
        # print(acc_results)
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
        #plt.show()

    except:
        print("Failed to fetch DSSP output for: " + id)
        dssp.close()


# # API:
# for id in my_pdb_ids[:1]:
#     data = {'pdb_id' : id,
#             'seq' : "",
#             'file_' : '(binary)'}

#     rest_url = 'https://www3.cmbi.umcn.nl/xssp/'
#     url_create = '{}api/create/pdb_id/dssp/'.format(rest_url)

#     r = requests.post(url_create, data=data)
#     r.raise_for_status()
#     job_id = json.loads(r.text)['id']
#     print("Job submitted successfully. Id is: '{}'".format(job_id))

#     ready = False
#     while not ready:
#         url_status = '{}api/status/pdb_id/dssp{}/'.format(rest_url, job_id)
#         r = requests.get(url_status)
#         r.raise_for_status()

#         status = json.loads(r.text)['status']
#         print("Job status is: '{}'".format(status))

#         if status == 'SUCCESS':
#             ready = True
#         elif status in ['FAILURE', 'REVOKED']:
#             raise Exception(json.loads(r.text)['message'])
#         else:
#             time.sleep(5)
#     else:

#         url_result = '{}api/result/pdb_id/dssp/{}/'.format(rest_url, job_id)
#         r = requests.get(url_result)
#         r.raise_for_status()
#         result = json.loads(r.text)['result']

#         print(result)