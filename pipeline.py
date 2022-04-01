from selenium import webdriver
from selenium.webdriver.common.by import By
import scrapy
from scrapy import Spider
from scrapy.http import FormRequest
from scrapy.utils.project import get_project_settings
from scrapy.utils.log import configure_logging
from scrapy.crawler import CrawlerRunner
from twisted.internet import reactor, defer

"""
swissprot_id = 'P0DTC2'
bepipred = 'http://tools.iedb.org/bcell/'
other = 'http://tools.iedb.org/discotope/'


driver = webdriver.Chrome(executable_path="chromedriver")
driver.maximize_window()
driver.get(bepipred)

otherd = webdriver.Chrome(executable_path="chromedriver")
otherd.maximize_window()
otherd.get(other)

submit_protein = driver.find_element(By.NAME, "swissprot")
submit_protein.send_keys(swissprot_id)
submit = driver.find_element(By.NAME, "submit")
submit.click()

columns = []
for i in range(1, 6):
    columns.append(driver.find_element(By.XPATH, "/html/body/div[3]/table[2]/thead/tr/th["+ str(i)+ "]").text)
print(columns)

cells = []
#body = driver.find_element(By.XPATH, "/html/body/div[3]/table[2]/tbody")
for tr in driver.find_elements(By.XPATH, "/html/body/div[3]/table[2]/tbody"):
    tds = tr.find_elements(By.TAG_NAME, 'td')
    if tds:
        cells.append([td.text for td in tds])
print(cells)
"""

bepipred = []
emini = []
kolaskar_tongaonkar = []
bepipred2 = []
choufasman = []
karplusschulz = []
parker = []
discotope = []
ellipro = []
mhcii = []
netctlpan = []

fasta = '''>West Nile virus envelope glycoprotein
        FNCLGMSNRDFLEGVSGATWVDLVLEGDSCVTIMSKDKPTIDVKMMNMEAANLAEVRSYCYLATVSDLST
        KAACPTMGEAHNDKRADPAFVCRQGVVDRGWGNGCGLFGKGSIDTCAKFACSTKAIGRTILKENIKYEVA
        IFVHGPTTVESHGNYSTQVGATQAGRFSITPAAPSYTLKLGEYGEVTVDCEPRSGIDTNAYYVMTVGTKT
        FLVHREWFMDLNLPWSSAGSTVWRNRETLMEFEEPHATKQSVIALGSQEGALHQALAGAIPVEFSSNTVK
        LTSGHLKCRVKMEKLQLKGTTYGVCSKAFKFLGTPADTGHGTVVLELQYTGTDGPCKVPISSVASLNDLT
        PVGRLVTVNPFVSVATANAKVLIELEPPFGDSYIVVGRGEQQINHHWHKSGSSIGKAFTTTLKGAQRLAA
        LGDTAWDFGSVGGVFTSVGKAVHQVFGGAFRSLFGGMSWITQGLLGALLLWMGINARDRSIALTFLAVGG
        VLLFLSVNVHA'''
alleles_mhcii = ['HLA-DRB1*01:01','HLA-DRB1*03:01', 'HLA-DRB1*04:01']
length_mhcii = ['12','13']
alleles_mhci = ['','','']
length_mhci = ['','']

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


configure_logging()
settings = get_project_settings()
runner = CrawlerRunner(settings)


@defer.inlineCallbacks
def crawl():
    #yield runner.crawl(Ellipro)
    #yield runner.crawl(Discotope)
    #yield runner.crawl(Bepipred2)
    #yield runner.crawl(Bepipred)
    #yield runner.crawl(Emini)
    #yield runner.crawl(Kolaskar)
    #yield runner.crawl(Chou_Fasman)
    #yield runner.crawl(Karplus_Schulz)
    #yield runner.crawl(Parker)
    #yield runner.crawl(MhcII)
    yield runner.crawl(Netctlpan)
    reactor.stop()


crawl()
reactor.run()  # the script will block here until the last crawl call is finished
