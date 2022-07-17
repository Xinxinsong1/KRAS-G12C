mergedcounts<- read.table("Xinxin_merged_gene_counts.txt", header = TRUE, quote = '\t',skip =0)
Untreat<- mergedcounts[,c(1,2,3,4,5,9,10,11)]
sampleNames <- c("Geneid","Symbol","C1","C2","C3","R1","R2","R3")
names(Untreat)[1:8] <- sampleNames
head(Untreat)
write.table(Untreat,"Untreat.txt", row.names = FALSE, quote = FALSE)
write.csv(Untreat,"Untreat.csv", row.names = FALSE, quote = FALSE)



