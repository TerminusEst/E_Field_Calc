#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
from bs4 import BeautifulSoup, SoupStrainer
from spacepy import pycdf
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import codecs
from urllib.request import urlopen
##########################################################################################
output_folder = "/home/blake/Desktop/Iono_output_temp/"

switch = 1
if switch == True:
    page = "https://ccmc.gsfc.nasa.gov/RoR_WWW/output_files/SWPC_SWMF_030311_1b/IONO-2D_CDF/"
    browser = webdriver.Chrome(executable_path=r'/home/blake/Drive/NASA_Shite/CDF_TEST/chromedriver')
    browser.get(page)

    html_source = browser.page_source

# download CDF files
soup = BeautifulSoup(html_source, "lxml")
links = []
for link in soup.findAll('a'):
    newlink = link.get('href')
    
    if ".swmf.it" in newlink:
        print(newlink)
        browser.get(page + newlink)
        time.sleep(0.1)

# test file
#cdf = pycdf.CDF("/home/blake/Downloads/null.swmf.i_e20120723-162900-000.cdf")

"""
# download Mag Data
count = 0
lines = html_source.split("\n")
for l in lines:
    if "mag_grid" in l:
        file_thing = l.split('a href="')[1].split('"')[0]
        print(file_thing)
        count += 1
        
        html = urlopen(page + file_thing).read().decode('utf-8')

        output_file = output_folder + file_thing[:-3] + ".txt"

        f = open(output_file, 'w')
        f.write(html)
        f.close()
"""        














