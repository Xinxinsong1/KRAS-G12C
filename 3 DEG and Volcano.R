res = read.csv("res4_log_order_R_vs_C.csv")
res2 <-subset (res,pvalue >10e-300)
head(res)
head(res2)

#提取差异基因分析的结果
diff_gene_deseq2<-subset(res,padj<0.05&abs(log2FoldChange) > 0)
dim(diff_gene_deseq2)
head(diff_gene_deseq2)
write.csv(diff_gene_deseq2,file = "DEG_CM_vs_C.csv") #输出为一个文件

#火山图
library(ggplot2)
library(EnhancedVolcano)
library(airway)
library(magrittr)
library(ggrepel)
EnhancedVolcano(res2,
                lab =res2$external_gene_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                title = 'R vs C',
                pointSize = 2.0,
                labSize = 3.0,
                colAlpha = 0.75,
                legendLabSize = 8)
#火山图run到这里就可以了


#备用
EnhancedVolcano(res2,
                lab = rownames(res2),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c('ENSG00000102794') ,#IRG1
                drawConnectors = TRUE,
                widthConnectors = 0.4,
                colConnectors = 'black',
                xlim = c(-15, 15),
                title = 'M155 0h vs 6h volcano',
                pCutoff = 10e-2,
                FCcutoff = 2,
                col=c('black', 'blue', 'green', 'red1'),
                colAlpha = 1,
                legend=c('NS','Log2 FC','P value',
                         'P value & Log2 FC'),
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 5.0)