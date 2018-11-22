# SELENIUM web scraper 
# used to download a shit ton of INTERMAGNET data files

import time
import sys
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import zipfile
import os
from sh import gunzip
####################################################################

browser = webdriver.Chrome(executable_path='/home/blake/Drive/NASA_Shite/DST_downloader/chromedriver_244')
browser.get("http://intermagnet.org/data-donnee/download-eng.php")

switch = input("\nDo you want to continue? (y/n)\n")
if switch != "y":
    sys.exit()
####################################################################

input_csv = "/home/blake/Drive/NASA_Shite/DST_downloader/Storm_Info.csv"
input_data = np.loadtxt(input_csv, skiprows = 1, dtype = str, delimiter=",")

for line in input_data[21:]:
    input_date = line[0]
    input_days = int(line[1])
    input_switch = int(line[3])
    
    if input_switch == 1:
        continue
    
    input_datetime = datetime.datetime(int(input_date[:4]), int(input_date[5:7]), int(input_date[8:]), 0, 0)

    start_date = input_datetime + datetime.timedelta(days = -1)
    end_date = input_datetime + datetime.timedelta(days = input_days-1)

    start_year = start_date.year
    start_month = start_date.month
    start_day = start_date.day

    end_year = end_date.year
    end_month = end_date.month
    end_day = end_date.day

    print("Trying to Download:")
    print(input_date)
    print("From: ", start_year, start_month, start_day)
    print("To:   ", end_year, end_month, end_day)
    print("\n")

    browser.get("http://intermagnet.org/data-donnee/download-eng.php#view")
    ####################################################################
    
    # first, change the dates:
    browser.find_element_by_name("from_year").click()
    time.sleep(1)
    browser.find_element_by_name("from_year").send_keys(str(start_year))
    time.sleep(1)
    browser.find_element_by_name("to_year").click()
    time.sleep(1)
    browser.find_element_by_name("to_year").send_keys(str(end_year))

    time.sleep(1)
    browser.find_element_by_name("from_month").click()
    time.sleep(1)
    browser.find_element_by_name("from_month").send_keys("%02d" % start_month)
    time.sleep(1)
    browser.find_element_by_name("to_month").click()
    time.sleep(1)
    browser.find_element_by_name("to_month").send_keys("%02d" % end_month)

    time.sleep(1)
    browser.find_element_by_name("from_day").click()
    time.sleep(1)
    browser.find_element_by_name("from_day").send_keys("%02d" % start_day)
    time.sleep(1)
    browser.find_element_by_name("to_day").click()
    time.sleep(1)
    browser.find_element_by_name("to_day").send_keys("%02d" % end_day)
    time.sleep(1)

    # Scroll to bottom
    browser.execute_script("window.scrollTo(0, document.body.scrollHeight*2);")

    # Press the "Search for data" button
    browser.find_element_by_id("search").click()

    # Select all Data
    time.sleep(5)
    browser.find_element_by_id("select_all").click()

    # Download all data
    time.sleep(5)
    browser.find_element_by_id("download_data").click()

    # Input email address
    time.sleep(5)
    browser.find_element_by_id("email").send_keys("sean.blake@nasa.gov")

    # Accept terms and conditions
    time.sleep(1)
    browser.find_element_by_id("accept").click()

    # now click the download button
    time.sleep(5)
    browser.execute_script("window.scrollTo(0, document.body.scrollHeight*100);")
    time.sleep(1)

    xpath_string = '//*[@id="wb-main-in"]/div[1]/form/input['
    for i in range(100, 500):
        new_xpath_string = xpath_string + str(i) + "]"
        try:
            aaa = browser.find_element_by_xpath(new_xpath_string)
            aaa.click()
            break
        except:
            continue
    time.sleep(1)

    # Now download the zipped folder
    time.sleep(10)
    html_lines = browser.page_source.split("\n")
    for i in html_lines:
        if "Download ZIP file" in i:
            break
    download_string = i.split(">")[2].split("<")[0]
    browser.find_element_by_link_text(download_string).click()

    time.sleep(30)

    ####################################################################
    print(input_date, "unzipping")
    # Now to unzip 
    downloadfolder = "/home/blake/Downloads/"
    newfolder = "/home/blake/Drive/NASA_Shite/DST_downloader/INTERMAGNET_DATA/" + input_date.replace("-", "_")

    switch = True
    while switch == True:
        filenames = os.listdir(downloadfolder)
        fname = downloadfolder + filenames[0]
        if "crdownload" in fname:
            time.sleep(1)
            continue
        else:
            switch = False
            
    zip_ref = zipfile.ZipFile(fname)
    zip_ref.extractall(newfolder)
    zip_ref.close()

    os.remove(newfolder + "/conditions_of_use.txt")
    os.remove(newfolder + "/log.txt")

    for iii in os.listdir(newfolder):
        if ".gz" in iii:
            gunzip(newfolder + "/" + iii)

    os.remove(downloadfolder + filenames[0])

    browser.get("http://intermagnet.org/data-donnee/download-eng.php#view")





























