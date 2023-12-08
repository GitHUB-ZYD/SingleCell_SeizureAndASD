library(VennDiagram)

x = read.csv("kegg.csv",header = 1)

ASD = x[x$group=="ASD",1]
Seizure = x[x$group=="Seizure",1]
Comorbidity = x[x$group=="Comorbidity",1]


venn.diagram(x=list(ASD,Seizure,Comorbidity),
             
             scaled = F, # 根据比例显示大小
             
             alpha= 1, #透明度
             
             lwd=1,lty=1, #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             
             cex = 2, # 数字大小
             
             fontface = "bold",  # 字体粗细；加粗bold
             
             fill=c('#FFFFCC','#CCFFFF',"#FFCCCC"), # 填充色 配色https://www.58pic.com/
             
             category.names = c("ASD", "Seizure","Comorbidity") , #标签名
             
             cat.dist = 0, # 标签距离圆圈的远近
             
             cat.pos = c(-120, -240, -180), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             
             cat.cex = 2, #标签字体大小
             
             cat.fontface = "bold",  # 标签字体加粗
             
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             
             output=TRUE,
             
             filename='./result.png',# 文件保存
             
             imagetype="png",  # 类型（tiff png svg）
             
             resolution = 300
             
)

# 获取交集
inter <- get.venn.partitions(list(ASD,Seizure,Comorbidity))
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.csv(inter[-c(5, 6)], 'result.txt', row.names = FALSE, sep = '\t', quote = FALSE)
