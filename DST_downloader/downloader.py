"""Crude script to scrape DST values from http://wdc.kugi.kyoto-u.ac.jp/

Have to download from dst_final/ and dst_provisional/ and then stitch together

Output is: 1 column of time floats, 24 rows of DST values (one per hour), one row per day
"""


import time
from bs4 import BeautifulSoup, SoupStrainer
from spacepy import pycdf
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import codecs
from urllib.request import urlopen
from draco import EField as ef
################################################################################

datez, dst_values = [], []

for year_int in range(2015, 2017, 1):
    for month_int in range(1, 13, 1):
    
        print(year_int, month_int)
        
        year = str(year_int)
        month = "%02d" % month_int

        page = "http://wdc.kugi.kyoto-u.ac.jp/dst_provisional/" + year + month + "/index.html"
        html = urlopen(page).read().decode('utf-8').split("\n")

        for index, line in enumerate(html):
            if "DAY" in line:
                start_index = int(index+1)
                break

        for i in range(start_index+1, 100):
            if "<" in html[i]:
                end_index = int(i-1)
                break
                

        for i in range(start_index, end_index):
            day_values = []
            line = html[i].replace("-", " -").split(" ")
            line = [int(x) for x in line if x!= ""]
            
            if (len(line) == 1) or (len(line) == 0):
                continue
            else:
                datez.append(datetime.datetime(int(year), int(month), int(line[0]), 0, 0))
                day_values.extend(line[1:])
                dst_values.append(day_values)

        time.sleep(0.5)

################################################################################
# Stitch together the data and save it to a file
"""
dst_values = np.array(dst_values).astype(float)
datez_float = np.array([ef.Time2Float(datez)]).T
output = np.concatenate((datez_float, dst_values), axis = 1)
savetxt("/home/blake/Drive/NASA_Shite/DST_downloader/DST_VALUES_2015_2016.txt", output)

data1 = np.loadtxt("/home/blake/Drive/NASA_Shite/DST_downloader/DST_VALUES_1991_2014.txt")
data2 = np.loadtxt("/home/blake/Drive/NASA_Shite/DST_downloader/DST_VALUES_2015_2016.txt")
data_final = np.concatenate((data1, data2), axis = 0)
savetxt("/home/blake/Drive/NASA_Shite/DST_downloader/DST_VALUES_1991_2016.txt", data_final)
"""

# Make a purty plot of minimum DST values per day
clf()

plot(ef.Float2Time(data_final[:,0]), np.min(data_final[:,1:], axis = 1))

xlabel("TIME", fontsize = 24)
ylabel("Minimum Daily DST (nT)", fontsize = 24)
show()
























