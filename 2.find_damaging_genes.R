setwd("/Users/shenjingting/Desktop/face_database/we_WGS/1.analysis/") 
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

#读取case
case_syn <- read_excel("1.zoutput/case_output.xlsx", sheet=1) 
case_LGD <- read_excel("1.zoutput/case_output.xlsx", sheet=2) 
case_Dmis <- read_excel("1.zoutput/case_output.xlsx", sheet=3) 
case_Dmis_HC <- read_excel("1.zoutput/case_output.xlsx", sheet=4)
case_LGD["type"] <- "LGD"
case_Dmis["type"] <- "Dmis"
case_Dmis_HC["type"] <- "Dmis_HC"
case_syn["type"] <- "synonymous"
case_damaging <- rbind(case_LGD, case_Dmis, case_Dmis_HC)
case_damaging <- case_damaging[!duplicated(case_damaging[,1:(ncol(case_damaging)-1)]), ]
#write.table(case_damaging, file="3.zoutput/case_damaging.txt", sep="\t", col.names = F, quote = F)
#write.table(case_syn, file="3.zoutput/case_synonymous.txt", sep="\t", col.names = F, quote = F)

#读取control
control_syn <- read_excel("1.zoutput/control_output.xlsx", sheet=1) 
control_LGD <- read_excel("1.zoutput/control_output.xlsx", sheet=2) 
control_Dmis <- read_excel("1.zoutput/control_output.xlsx", sheet=3) 
control_Dmis_HC <- read_excel("1.zoutput/control_output.xlsx", sheet=4) 
control_LGD["type"] <- "LGD"
control_Dmis["type"] <- "Dmis"
control_Dmis_HC["type"] <- "Dmis_HC"
control_syn["type"] <- "synonymous"
control_damaging <- rbind(control_LGD, control_Dmis,control_Dmis_HC)
control_damaging <- control_damaging[!duplicated(control_damaging[,1:(ncol(control_damaging)-1)]), ]
write.table(control_damaging, file="3.zoutput/control_damaging.txt", sep="\t", col.names = F, quote = F)
write.table(control_syn, file="3.zoutput/control_synonymous.txt", sep="\t", col.names = F, quote = F)

##########################对特定基因的变异统计###########################################
#筛选CNCC formation 相关基因变异
SOX9 <- subset(case_damaging, SYMBOL=="SOX9")
PAX3 <- subset(case_damaging, SYMBOL=="PAX3")
TFAP2A <- subset(case_damaging, SYMBOL=="TFAP2A")
TCOF1 <- subset(case_damaging, SYMBOL=="TCOF1")
POLR1C <- subset(case_damaging, SYMBOL=="POLR1C")
POLR1D <- subset(case_damaging, SYMBOL=="POLR1D")
RPS19BP1 <- subset(case_damaging, SYMBOL=="RPS19BP1")
DHODH <- subset(case_damaging, SYMBOL=="DHODH")
SF3B4 <- subset(case_damaging, SYMBOL=="SF3B4")

#筛选CNCC formation 相关基因变异
SNAIL2 <- subset(case_damaging, SYMBOL=="SNAIL2")
EDNRB <- subset(case_damaging, SYMBOL=="EDNRB")
EDN3 <- subset(case_damaging, SYMBOL=="EDN3")
MID1 <- subset(case_damaging, SYMBOL=="MID1")
ELP1 <- subset(case_damaging, SYMBOL=="ELP1")

#筛选CNCC differentiation 相关基因变异
PHOX2B <- subset(case_damaging, SYMBOL=="PHOX2B")
MITF <- subset(case_damaging, SYMBOL=="MITF")
SOX10 <- subset(case_damaging, SYMBOL=="SOX10")
KIT <- subset(case_damaging, SYMBOL=="KIT")

#筛选半侧颜面短小致病基因变异
OTX3 <- subset(case_damaging, SYMBOL=="OTX3")

###############################case和control中携带有害突变基因的人数统计###################################
Damaging_gene <- read.table("3.zoutput/case_damaging_gene.list")
case_damaging_dedup <- distinct(case_damaging[c("file", "SYMBOL")])
control_damaging_dedup <- distinct(control_damaging[c("file", "SYMBOL")])
Damaging_gene_case_carrier <- case_damaging_dedup %>% group_by(SYMBOL) %>% summarise(carriers = n())
Damaging_gene_control_carrier <- control_damaging_dedup %>% group_by(SYMBOL) %>% summarise(carriers = n())
write.table(Damaging_gene_case_carrier, file="3.zoutput/Damaging_gene_case_carrier.txt", sep="\t", col.names = F, row.names = T, quote = F)
write.table(Damaging_gene_control_carrier, file="3.zoutput/Damaging_gene_control_carrier.txt", sep="\t", col.names = F, row.names = T, quote = F)



################对damaging mutation基因的功能富集分析########################
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(forcats)
library(dplyr)
library(stringr)
library(tidyr)
library(org.Rn.eg.db)
library(tidyverse)
Damaging_gene <- read.table("3.zoutput/case_damaging_gene.list")
#【教训：这一步非常重要】
#Damaging_gene <- subset(Damaging_gene, Damaging_gene$V2 >= 3)
Damaging_gene_list <- Damaging_gene$V1
print(Damaging_gene_list)
#转化为ENTREZID
Damaging_gene_list <- bitr(Damaging_gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb = "org.Hs.eg.db")

#GO 富集
ego_Damaging <- enrichGO(gene = Damaging_gene_list$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.05,
                    readable = TRUE)
res_go <- data.frame(ego_Damaging)
write.csv(res_go, "3.zoutput/enrich_result/Damaging_gene_go_enrich.csv", quote = F)

# 计算Fold Enrichment的函数
calculate_fold_enrichment <- function(gene_ratio, bg_ratio) {
  # 解析GeneRatio
  gene_parts <- as.numeric(str_split(gene_ratio, "/", simplify = TRUE))
  gene_ratio_val <- gene_parts[1] / gene_parts[2]
  
  # 解析BgRatio
  bg_parts <- as.numeric(str_split(bg_ratio, "/", simplify = TRUE))
  bg_ratio_val <- bg_parts[1] / bg_parts[2]
  
  # 计算富集倍数
  fold_enrichment <- gene_ratio_val / bg_ratio_val
  
  return(fold_enrichment)
}

# 应用计算函数
res_go <- res_go %>%
  mutate(
    FoldEnrichment = map2_dbl(
      GeneRatio, 
      BgRatio, 
      ~ calculate_fold_enrichment(.x, .y)
    ),
    # 可选：添加-log10(padjust)用于可视化
    logPvalue = -log10(p.adjust)  
  )

table(ego_Damaging@result$ONTOLOGY)
dotplot(ego_Damaging, showCategory = 5, split='ONTOLOGY', title="GO") + facet_grid(ONTOLOGY~., scale="free") 

#提取BP
BP <- subset(res_go, res_go$ONTOLOGY=="BP")
rownames(BP) <- 1:nrow(BP)
BP$order <- factor(rev(as.integer(rownames(BP))), labels = rev(BP$Description))
BP <- BP[1:15, ]
ID_list <- BP$ID
###############################20260218新加入代码：加入钙离子信号通路###########
ID_list <- ID_list[-15]
ID_list <- ID_list[-8]
ID_list
ID_list[15] <- "GO:0097553"
ID_list
BP <- subset(res_go, res_go$ID %in% ID_list)
rownames(BP) <- 1:nrow(BP)
BP$order <- factor(rev(as.integer(rownames(BP))), labels = rev(BP$Description))
################################################################################
#作图
ggplot(BP,aes(y=order,x=FoldEnrichment))+
  
  geom_point(aes(size=Count,color=logPvalue))+
  
  scale_color_gradient(low = "purple",high ="yellow")+
  
  labs(color = "-log10(P.adj)",size="Count",
       
       x="Fold Enrichment",y="GO term",title="GO Enrichment") +
  
  theme_bw( base_size = 16 )

#########################对LGD和Dmis的分别进行GO富集#############################
LGD_gene_list <- unique(case_LGD$SYMBOL)
Dmis_gene_list <- unique(case_Dmis$SYMBOL)
LGD_gene_list <- bitr(LGD_gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb = "org.Hs.eg.db")
Dmis_gene_list <- bitr(Dmis_gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb = "org.Hs.eg.db")

ego_LGD <- enrichGO(gene = LGD_gene_list$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.05,
                    readable = TRUE)
res_LGD <- data.frame(ego_LGD)

ego_Dmis <- enrichGO(gene = Dmis_gene_list$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.05,
                    readable = TRUE)
res_Dmis <- data.frame(ego_Dmis)

#输出结果
write.csv(res_LGD, "3.zoutput/enrich_result/LGD_gene_go_enrich.csv", quote=F)
write.csv(res_Dmis, "3.zoutput/enrich_result/Dmis_gene_go_enrich.csv", quote=F)

#各自提取top20
res_LGD <- read.csv("3.zoutput/enrich_result/LGD_gene_go_enrich.csv")
res_Dmis <- read.csv("3.zoutput/enrich_result/Dmis_gene_go_enrich.csv")
table(res_LGD$ONTOLOGY) #BP CC MF
table(res_Dmis$ONTOLOGY) #BP CC MF
#lgd - top20-
res_LGD <- res_LGD[1:20, ]
res_LGD$type <- "LGD"
#dmis -top20 -
#加上GO:0006816
GO_item = res_Dmis[res_Dmis$ID =="GO:0006816", ]
colnames(GO_item) <- colnames(res_Dmis)
res_Dmis <- res_Dmis[1:20, ]
res_Dmis <- rbind(res_Dmis, GO_item)
res_Dmis$type <- "Dmis+Dmis-HC"

res_all <- rbind(res_LGD, res_Dmis)
res_all["logP"] <- -log10(res_all["p.adjust"])
res_all$Count <- as.numeric(res_all$Count)

#作图
ggplot(res_all,aes(y=Description,x=type))+
  geom_point(aes(size=Count,color=logP))+
  
  scale_color_gradient(low = "purple",high ="yellow")+
  
  labs(
       
       x="Fold Enrichment",y="GO term",color="-log10 (P-value)")+
  
  theme_bw(base_size = 16)


#优化作图
res_LGD_Dmis <- read_excel("3.zoutput/enrich_result/merge_LGD_Dmis_gene_go_enrich.xlsx", sheet = 1)
colnames(res_LGD_Dmis)
res_LGD_Dmis["p.adjust"]
res_LGD_Dmis["p.adjust"] <- as.numeric(res_LGD_Dmis$p.adjust)
res_LGD_Dmis["logP"] <- -log10(res_LGD_Dmis["p.adjust"])

ggplot(res_LGD_Dmis,aes(x=type, y=Description))+
  geom_point(aes(size=Count,color=logP))+
  
  scale_color_gradient(low = "purple",high ="yellow")+
  
  labs(x="Type",y="GO term", size = "Count", color = "-log10(P.adj)", title="GO Enrichment")+
  
  theme_bw(base_size = 14)



################对damaging mutation基因的kegg通路富集分析########################
ego_kegg <- enrichKEGG(gene = Damaging_gene_list$ENTREZID, 
                       organism = "hsa",
                       keyType = "kegg",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.01,
                       qvalueCutoff = 0.05)
res_kegg <- as.data.frame(ego_kegg)

#输出结果
write.csv(res_kegg, "3.zoutput/enrich_result/Damaging_gene_kegg_enrich.csv", quote=F)

res_kegg <- res_kegg[1:20, ]
# 绘制柱状图
ggplot(res_kegg, aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "white",high ="purple")+
  labs(title = "KEGG Pathway Enrichment", x = "Pathways", y = "Gene Count") +
  theme_minimal()


##########对significant damaging mutation (q<0.3) 基因的go通路富集 #############
Damaging_gene <- read.table("./5.zoutput/case_damaging_gene_pvalue0.3.txt")
#【教训：这一步非常重要】
#Damaging_gene <- subset(Damaging_gene, Damaging_gene$V2 >= 3)
Damaging_gene_list <- Damaging_gene$V1
print(Damaging_gene_list)
#转化为ENTREZID
Damaging_gene_list <- bitr(Damaging_gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb = "org.Hs.eg.db")

#GO 富集
ego_Damaging <- enrichGO(gene = Damaging_gene_list$ENTREZID,
                         OrgDb = org.Hs.eg.db,
                         ont = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.05,
                         readable = TRUE)
res_go <- data.frame(ego_Damaging)
