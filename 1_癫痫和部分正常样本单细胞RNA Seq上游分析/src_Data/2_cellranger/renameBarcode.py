samplename = ["5607_BA9","5981_BA9_Nova","5546_BA9_Nova","5839_BA9","5893_BA9","5881_BA9","5873_BA9","5609_BA9_Nova","5787_BA9_Nova","5892_BA9","5850_BA9","5929_BA9","5876_BA9"]
barcodes = []
with open("merge_barcodes.tsv","r") as f:
    tmp = f.readlines()
    for i in tmp:
        barcodes.append(i)
renameBarcodes = []



with open("rename_barcode.tsv","w") as f:
    for i in range(13):
        endNum = "-"+str(i+1)+"\n"
        tmp = filter(lambda x:endNum in x,barcodes)
        for j in list(tmp):
            j = j.rstrip()+"_"+samplename[i]
            renameBarcodes.append(j)
            f.write(j+"\n")

    


