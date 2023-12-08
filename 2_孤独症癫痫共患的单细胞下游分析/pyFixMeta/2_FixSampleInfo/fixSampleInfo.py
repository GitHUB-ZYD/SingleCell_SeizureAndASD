import pandas as pd
import numpy as np
meta = pd.read_csv("ASD_Seizure_meta_fixCellType.csv")
sampleInfo = pd.read_csv("sampleInfo.csv")

sampleID = list(sampleInfo["Sample name"])
sampleInfoCol = list(sampleInfo.columns)


for i in sampleID:
    info = np.array(sampleInfo.loc[sampleInfo["Sample name"]==i,])[0]
    condition = meta["cellID"].str.contains(i)
    for j in range(len(sampleInfoCol)):
        meta.loc[condition,sampleInfoCol[j]]=info[j]
meta.to_csv("ASD_Seizure_meta_fixSampleInfo.csv")
