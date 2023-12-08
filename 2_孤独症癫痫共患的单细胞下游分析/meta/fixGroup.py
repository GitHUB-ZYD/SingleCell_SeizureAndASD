import pandas as pd

cellType = ["Neu-NRGN","L5/6","Oligodendrocytes","OPC","AST","Endothelial","Microglia","IN-VIP","L5/6-CC","IN-SV2C","L2/3","IN-PV","L4","IN-SST","Neu-mat"]


meta = pd.read_csv("ASD_Seizure_meta_fix.csv")

condition = (meta["Seizure"]=="no") & (meta["ASD"]=="no")
meta.loc[condition,"group"]="Control"

condition = (meta["Seizure"]=="no") & (meta["ASD"]=="yes")
meta.loc[condition,"group"]="ASD"

condition = (meta["Seizure"]=="yes") & (meta["ASD"]=="no")
meta.loc[condition,"group"]="Seizure"

condition = (meta["Seizure"]=="yes") & (meta["ASD"]=="yes")
meta.loc[condition,"group"]="Comorbidity"

meta.to_csv("ASD_Seizure_meta_fix_fixGroup.csv")