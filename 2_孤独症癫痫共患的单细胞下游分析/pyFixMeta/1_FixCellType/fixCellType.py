import pandas as pd

meta = pd.read_csv("ASD_Seizure_meta.csv")
meta["eyeCluster"] = meta["cluster"]
celltype = pd.read_csv("eye.csv")
celltype = list(celltype["celltype"])

clusters = set(meta["seurat_clusters"])

for i in range(len(celltype)):
    condition = (meta["cluster"].isna())&(meta["seurat_clusters"]==i)
    meta.loc[condition,"eyeCluster"]=celltype[i]
meta.to_csv("ASD_Seizure_meta_fixCellType.csv")