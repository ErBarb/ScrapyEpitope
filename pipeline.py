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
discotope = []
ellipro = []


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
            callback=self.results_page)

    def results_page(self, response):
        url = 'http://tools.iedb.org/bcell/result/'
        yield scrapy.Request(url=url, callback=self.get_results, dont_filter=True)

    def get_results(self, results):
        row = []
        for cell in results.xpath('/html/body/div[3]/table[2]/tbody/tr'):
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
            callback=self.results_page)

    def results_page(self, response):
        url = 'http://tools.iedb.org/bcell/result/'
        yield scrapy.Request(url=url, callback=self.get_results, dont_filter=True)

    def get_results(self, results):
        row = []
        for cell in results.xpath('/html/body/div[3]/table[2]/tbody/tr'):
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
                'swissprot': 'P0C6U8',
                'sequence_text': '',
                'method': 'Bepipred2'
            },
            callback=self.results_page)

    def results_page(self, response):
        url = 'http://tools.iedb.org/bcell/result/'
        yield scrapy.Request(url=url, callback=self.get_results, dont_filter=True)

    def get_results(self, results):
        row = []
        for cell in results.xpath('/html/body/div[3]/table[2]/tbody/tr'):
            if int(cell.xpath('td[5]//text()').extract_first()) >= 7:
                row.append(cell.xpath('td[1]//text()').extract_first())
                row.append(cell.xpath('td[2]//text()').extract_first())
                row.append(cell.xpath('td[3]//text()').extract_first())
                row.append(cell.xpath('td[4]//text()').extract_first())
                row.append(cell.xpath('td[5]//text()').extract_first())
                bepipred2.append(row)
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

    #def results_page(self, response):
    #    sec_page_url = 'http://tools.iedb.org/ellipro/result/chains/'
    #    yield scrapy.Request(url=sec_page_url, callback=self.second_page, dont_filter=True)

    def second_page(self, sec_page):
        print('works')
        print(sec_page)
        yield FormRequest.from_response(
            sec_page,
            #url=sec_page,
            formdata={
                'chain': ['A','B','C']
                #'chain': 'B',
                #'chain': 'C'
            },
            callback=self.results_page)

    def results_page(self, response):
        url = 'http://tools.iedb.org/ellipro/result/predict/'
        yield scrapy.Request(url=url, callback=self.get_results, dont_filter=True)

    def get_results(self, results):
        row = []
        for cell in results.xpath('/html/body/div[3]/table[2]/tbody/tr'):
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


configure_logging()
settings = get_project_settings()
runner = CrawlerRunner(settings)


@defer.inlineCallbacks
def crawl():
    #yield runner.crawl(Ellipro)
    #yield runner.crawl(Discotope)
    #yield runner.crawl(Bepipred2)
    yield runner.crawl(Bepipred)
    #yield runner.crawl(Emini)
    #yield runner.crawl(Kolaskar)
    reactor.stop()


crawl()
reactor.run()  # the script will block here until the last crawl call is finished

print('bepipred')
for s in bepipred:
    print(*s)
print('emini')
for s in emini:
    print(*s)
print('kolaskar_tongaonkar')
for s in kolaskar_tongaonkar:
    print(*s)
print('bepipred2')
for s in bepipred2:
    print(*s)
print('discotope')
for s in discotope:
    print(*s)
print('ellipro')
for s in ellipro:
    print(*s)
