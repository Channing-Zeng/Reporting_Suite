#!/usr/bin/python

from selenium import webdriver
# initiate the browser. It will open the url, 
# and we can access all its content, and make actions on it. 
browser = webdriver.Chrome()
url = 'file:///Users/vlad/Dropbox/az/analysis/TGEN/final/2014-07-16_exomes/qc/targQC/targQC.html'
# the page test.html is changing constantly its content by receiving sockets, etc. 
#So we need to save its "status" when we decide for further retrieval)
browser.get(url)
# wait until we want to save the content (this could be a buttonUI action, etc.):
#raw_input("Press to print web page")  
# save the html rendered content in that moment: 
html_source = browser.page_source
# display to check: 
with open('STATIC.html', 'w') as f:
  f.write(html_source)
#print html_source
