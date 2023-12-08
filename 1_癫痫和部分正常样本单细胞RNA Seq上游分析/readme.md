# 2023-11-08

1. 由于只有癫痫样本没有转录组表达数据，需要自己做单细胞上游，因此下载并且sra转化为fastq

[爬虫脚本(下载SRA的地址获取)，样本的信息](./src_data/0_数据下载和样本信息)

只下载了seizures的所有单细胞样本，因为有其他的ASD样本的基因表达矩阵，没必要重复跑

```bash
#多线程下载
#首先将urls_eplicy_all中的15个url平均分为3份，分别时url1,url2,url3
nohup wget -i url1 1>url1.log 2>&1 &
nohup wget -i url2 1>url2.log 2>&1 &
nohup wget -i url3 1>url3.log 2>&1 &
# 3-4个线程即可，不要太多，因为我不知道会不会速度更慢
```

[sra转化为fastq的脚本和运行时间统计](./src_Data/1_sra2fastq)

```bash
time fasterq-dump --split-3 -e 24 SRR9264381 -p
# syh说可以直接从ebi下载fastq文件，但是我去尝试了一下，只有200kb的速度，不如自己转！！
# -e指定24个线程，一个sra大概运行150-250 分钟
# 较大的线程数目在某一步会消耗大量内存，当我指定64个线程的时候，运行期间会内存不足死掉；指定20个线程的时候，峰值消耗了64g内存，fyp老师机器的50%；后面我都是指定24个线程运行
# 有多少个进程，在tmp文件中就会有多少组数据，如果是20，那么就会有00-19这些文件
# -p 可以输出fasterq的进度
## fasterq-dump没有--gzip参数，不能直接输出为fastq.gz文件，造成存储空间的浪费；后期可以通过pigz对其进行多线程压缩，速度嘎嘎快，比gzip快很多，但是消耗cpu运算资源



```

```bash
#遇到一个小问题
time ls
#输出是这样子的
real	0m0.003s
user	0m0.000s
sys	0m0.003s

nohup time ls & 
#这样子运行的话，就会失去对应的输出格式，变成下面这种，格式变了不好看
0.00user 0.00system 3:35:15elapsed 16%CPU (0avgtext+0avgdata 794064maxresident)k
155999312inputs+502334488outputs (0major+12816759minor)pagefaults 0swaps
# 然而在fasterq.log这个文件中，输出格式却还是很清楚
```



# 2023-11-09

### 开始cell ranger count

**需要通过软链接将fasqt.gz改名为cell ranger能够识别的格式**

比如

```bash
ln -s ../SRR9264381_1.fastq SRR9264381_S1_L001_R1_001.fastq
ln -s ../SRR9264381_2.fastq SRR9264381_S1_L001_R2_001.fastq
```

15个脚本的批量coun见[cellCount.sh](./src_Data\2_cellranger/0_cellCount.sh)

### 开始cell ranger  aggr

count完之后要合并，使用aggr先合并seizures和control的，然后再将合并后的和论文提供给的合并，使用seurat5

```bash
# 如果我是用绝对路径指定软链接的target，比如
ln -s /source/SRR_1.fastq /target/SRR_1.fastq
# 那么 当我修改了source的名字为source1，软链接失效

# 如果我使用了相对路径，然后做软链接，然后又将该软链接移动到其他目录，由于是相对路径，导致相对路径无法找到对应的文件，然后就链接失效
```

```bash
# cellranger aggr
 cellranger aggr --id=seizure_merge --csv=info_run.csv --localmem 64G
# 指定使用64G内存，避免使用过多内存影响别人
# 合并以后，编写python脚本，给barcode.tsv改名字，让它的命名格式和别人给的barcode文件是一样的
```

[info_run.csv](.\src_Data\2_cellranger\info_run.csv)

[改名脚本](.\src_Data\2_cellranger\renameBarcode.py)

### 论文给的矩阵数据和自己aggr得到的矩阵数据见[目录](.\src_Data\3_mtx)



**注意**

论文给的矩阵数据只包括了ASD和control的41个样本，没有包括癫痫的control和癫痫样本，并且，癫痫中使用的control样本有两个和论文给的矩阵重复了(5958_BA9,5577_BA9)

```python
# 统计有多少了细胞
import pandas as pd

asd = []
with open('exprMatrix_50.tsv', 'r', encoding='utf-8') as f:
    for line in f:
    # 读取一行后，末尾一般会有一个\n，所以用strip函数去掉
        line = line.strip('\n').split('\t')
        for i in line:
            i = i.split('-1_')[-1]
            asd.append(i)
        print(len(set(asd)))
        
        break
```

因此只需要将**SRR9264381**到**SRR9264393**跑一下单细胞上游，然后合并，然后再合并到总表即可。





# 2023-11-13

| BA24 | Ventral Anterior cingulate cortex | 腹前扣带皮层     | emotional and cognitive processing                   | ACC  |
| ---- | --------------------------------- | ---------------- | ---------------------------------------------------- | ---- |
| BA9  | Dorsolateral prefrontal cortex    | 背外侧前额叶皮层 | prefrontal associational integration                 | PFC  |
| BA46 | Dorsolateral prefrontal cortex    | 背外侧前额叶     | participates in prefrontal associational integration | PFC  |



## Seurat5中找高表达基因和差异表达基因的函数

1. **FindVariableFeatures**：
   - **功能**：`FindVariableFeatures` 用于识别在数据集中具有变异性的特征（通常是基因）。这些特征在细胞之间的表达水平变化较大，可能包括差异表达的基因。
   - **用途**：通常，`FindVariableFeatures` 的输出结果用于减少数据维度，选择最具信息量的特征进行后续分析，如主成分分析（PCA）或 t-分布随机近邻嵌入（t-SNE）等。它有助于减少计算复杂性并提高可视化的效果。
2. **FindMarkers**：
   - **功能**：`FindMarkers` 用于执行差异表达分析，比较不同组（通常是细胞群组或簇）之间的基因表达，以识别在不同组之间具有显著差异表达的基因。
   - **用途**：`FindMarkers` 的输出结果包括每个群组（或组别）中的差异表达基因列表，以及相关的统计信息（如p值和折叠变化）。这些结果有助于识别在不同细胞群组之间的生物学差异。
