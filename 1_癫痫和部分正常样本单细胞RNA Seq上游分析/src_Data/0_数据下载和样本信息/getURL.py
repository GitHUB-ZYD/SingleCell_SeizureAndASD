# coding=utf-8
import requests

from xml.dom.minidom import parse
SRR = []
downloadUrls = []
with open("SRR_Acc_List.txt","r") as f:
    x = f.readlines()
    for i in x:
        SRR.append(i.rstrip())
count=1

with open("downloadUrl.csv","w") as f1:
    f1.write("Accession,downloadUrl\n")
    for i in SRR:
        url = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=%s"%(i)
        #print(url)
        response = requests.get(url)
        with open("temp.xml","w",encoding="utf-8") as f2:
            f2.write(response.text)
        domTree = parse("temp.xml")
        #获得根节点
        rootNode = domTree.documentElement
        for j in rootNode.getElementsByTagName('SRAFile'):
            url = j.getAttribute('url')
            if "sra-pub-run-odp.s3.amazonaws.com" in url:
                #print(url)
                downloadUrls.append(url)
                f1.write("%s,%s\n"%(i,url))
                break
        print(count)
        count=count+1